/**
 *	Generate a bunch of manifold arcs from a quasi-halo and correct them to be natural in the BC4BP
 *
 *	To compile: g++ --std=c++11 -ltpat -Wall -pedantic LPF_PropAllManifolds.cpp -o a.out
 */

#include "AllIncludes.hpp"

#include <cmath>
#include <cstdio>

using namespace astrohelion;

int main(void){
	int numHaloNodes = 8;
	int numManNodes = 10;
	int numManifolds = 3;

	CorrectionEngine corrector;
	SysData_cr3bp seSys("sun", "earth");
	SysData_cr3bp emSys("earth", "moon");
	SysData_bc4bp bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");

	// Determine geometry at t0 so that we can start at T = 0;
	DynamicsModel_bc4bp::orientAtEpoch(leaveHaloEpoch, &bcSys);
	SimEngine bcEngine;

	// std::vector<double> quasihaloIC {0.992062701730402, -0.00193297804496478, -0.000855312860994067,
	// 	-0.00019220600470211, -0.0110053275009781, -0.000462973103726882};
	// double quasihaloPeriod = 3.07407568042871;

	// State from Quasi-halo family #19, orbit number 38, beginning at initRev = 5
	// See LPF_getHaloSegments script to generate the quasi-halo segment
	std::vector<double> quasihaloIC {0.992077032646842, -0.000673288630772272, -0.00076892548988272,
		-1.80048025175285e-05, -0.011893178120389, -0.000323997259648554};
	double quasihaloPeriod = 3.06658003613767;

	// Create a simulation engine and propagate the quasihalo for one rev
	SimEngine engine;
	Traj_cr3bp quasihalo(&seSys);
	engine.runSim(quasihaloIC, quasihaloPeriod, &quasihalo);
	
	// Compute manifolds from the halo for a short amount of time to get ICs
	std::vector<Traj_cr3bp> manICs = getManifolds(Manifold_tp::MAN_U_P, &quasihalo, numManifolds, 1e-4);

	// Create two events used to propagate manifolds to XZ plane
	double earthX = 1 - seSys.getMu();
	Event cross_EarthPlane(&seSys, Event_tp::YZ_PLANE, 1, true, &earthX);
	Event cross_Y0(&seSys, Event_tp::XZ_PLANE, 1, true);

	double bcMaxDist_data[] = {1, 2};
	Event bcMaxDist(&bcSys, Event_tp::DIST, 0, true, bcMaxDist_data);
	bcEngine.addEvent(bcMaxDist);

	for(size_t n = 0; n < manICs.size(); n++){
		printf("Manifold %03zu:\n", n);
		// First, propagate to X = X_earth to leave the quasihalo
		engine.clearEvents();
		engine.addEvent(cross_EarthPlane);
		Traj_cr3bp o(&seSys);
		engine.runSim(manICs[n].getStateByIx(0), 2*PI, &o);
		std::vector<Event> endEvents = engine.getEndEvents(&o);

		if(std::find(endEvents.begin(), endEvents.end(), cross_EarthPlane) == endEvents.end()){
			astrohelion::printErr("  Manifold %zu did not encounter the Earth-X plane...\n", n);
		}else{
			// If the manifold reached that plane, propagate to the XZ plane
			engine.clearEvents();
			engine.addEvent(cross_Y0);
			Traj_cr3bp o2(&seSys);
			engine.runSim(o.getStateByIx(-1), PI, &o2);
			endEvents = engine.getEndEvents(&o2);
			if(std::find(endEvents.begin(), endEvents.end(), cross_Y0) == endEvents.end()){
				astrohelion::printErr("  Manifold %zu did not reach the XZ plane...\n", n);
			}else{
				// If the XZ plane has been reached, concatenate the two arcs
				o += o2;
				
				// Turn the manifold into a nodeset
				Nodeset_cr3bp manData(o, numManNodes);

				// Create some nodes to represent the quasihalo
				Nodeset_cr3bp haloData(o.getStateByIx(0), &seSys, -PI, numHaloNodes);
				double haloTOF = std::abs(haloData.getTotalTOF());

				// Reverse order, tack the manifold on to the end
				// haloData.reverseOrder();
				// haloData.deleteNode(-1);
				// haloData += manData;
				haloData.appendSetAtNode(&manData, 0, 0, 0);
				haloData.print();
				haloData.printInChrono();

				// Transform to BCR4BP system
				double t0 = 0;// t = 0 occurs at the beginning of the manifold
				Nodeset_bc4bp bcNodes = bcr4bpr_SE2SEM(haloData, &bcSys, 0, t0);

				// Constrain the final node to be on the XZ_Plane
				double fixY0_Data[] = {NAN, 0, NAN, NAN, NAN, NAN};
				Constraint fixY0(Constraint_tp::STATE, bcNodes.getNumNodes()-1, fixY0_Data, 6);
				bcNodes.addConstraint(fixY0);

				bcNodes.print();
				bcNodes.printInChrono();
				waitForUser();

				// Correct to be continuous
				try{
					Nodeset_bc4bp bcCorrected(&bcSys);
					corrector.setTol(9e-12);
					corrector.multShoot(&bcNodes, &bcCorrected);

					bcCorrected.print();
					bcCorrected.printInChrono();

					// Now propagate forward for a while from the final state
					Traj_bc4bp natArc(&bcSys);
					bcEngine.runSim(bcCorrected.getStateByIx(-1), bcCorrected.getEpochByIx(-1), 360*24*3600/bcSys.getCharT(), &natArc);
					Traj_bc4bp fullTraj = Traj_bc4bp::fromNodeset(bcCorrected);
					fullTraj += natArc;

					// Save the trajectory to file
					char name[32];
					sprintf(name, "data/LPF_QH_4B_NaturalManifolds/Traj%03zu_SEM.mat", n);
					fullTraj.saveToMat(name);

					Traj_cr3bp fullTraj_SE = bcr4bpr_SEM2SE(fullTraj, &seSys);
					sprintf(name, "data/LPF_QH_4B_NaturalManifolds/Traj%03zu_SE.mat", n);
					fullTraj_SE.saveToMat(name);

					Traj_cr3bp fullTraj_EM = cr3bp_SE2EM(fullTraj_SE, &emSys, bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
					sprintf(name, "data/LPF_QH_4B_NaturalManifolds/Traj%03zu_EM.mat", n);
					fullTraj_EM.saveToMat(name);

					Traj_cr3bp fullTraj_ECI = cr3bp_rot2inert(fullTraj_EM, 0);
					sprintf(name, "data/LPF_QH_4B_NaturalManifolds/Traj%03zu_ECI.mat", n);
					fullTraj_ECI.saveToMat(name);

					Traj_cr3bp fullTraj_MCI = cr3bp_rot2inert(fullTraj_EM, 1);
					sprintf(name, "data/LPF_QH_4B_NaturalManifolds/Traj%03zu_MCI.mat", n);
					fullTraj_MCI.saveToMat(name);

				}catch(DivergeException &e){
					astrohelion::printErr("  Unable to correct manifold %03zu to be continuous in BC4BP... darn!\n", n);
				}
			}
		}
	}
}