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
	int numManifolds = 100;

	MultShootEngine corrector;
	SysData_cr3bp seSys("sun", "earth");
	SysData_cr3bp emSys("earth", "moon");
	SysData_bc4bp bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEphemerisTime("2016/04/18");

	// Define synodic month, seconds
	double P_syn = 2*PI*bcSys.getCharT()/(sqrt(bcSys.getMu()/pow(bcSys.getCharLRatio(), 3)) - bcSys.getK());

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
	Arcset_cr3bp quasihalo(&seSys);
	engine.runSim(quasihaloIC, quasihaloPeriod, &quasihalo);
	
	// Compute manifolds from the halo for a short amount of time to get ICs
	std::vector<Arcset_cr3bp> manICs = getManifolds(Manifold_tp::MAN_U_RIGHT, &quasihalo, numManifolds, 1e-4);

	// Create event that trigers when trajectory passes within lunar radius
	double evtData[] = {1, (1 - emSys.getMu())*emSys.getCharL()/seSys.getCharL()};
	Event enterLunarRad(&seSys, Event_tp::DIST, -1, true, evtData);

	// Create two events used to propagate manifolds to XZ plane
	double earthX = 1 - seSys.getMu();
	Event cross_EarthPlane(&seSys, Event_tp::YZ_PLANE, 1, true, &earthX);
	Event cross_Y0(&seSys, Event_tp::XZ_PLANE, 1, true);

	double bcMaxDist_data[] = {1, 2};
	Event bcMaxDist(&bcSys, Event_tp::DIST, 0, true, bcMaxDist_data);
	bcEngine.addEvent(bcMaxDist);

	for(unsigned int n = 0; n < manICs.size(); n++){
	// for(unsigned int n = 43; n < 44; n++){
		printf("Manifold %03zu:\n", n);

		// First propagate to lunar radius and determine epoch shift required to minimize lunar encounter distance
		Arcset_cr3bp traj(&seSys);
		engine.clearEvents();
		engine.addEvent(enterLunarRad);
		engine.runSim(manICs[n].getStateByIx(0), 4*PI, &traj);
		
		// Save the events (includes automatically-generated crash events)
		std::vector<Event> events = engine.getEvents();

		// Determine which events killed the simulation
		std::vector<Event> endEvents = engine.getEndEvents(&traj);
		
		double dT = 0;
		if(std::find(endEvents.begin(), endEvents.end(), enterLunarRad) == endEvents.end()){
			printErr("  Manifold %03zu did not reach lunar radius...\n", n);
			continue;
		}else{
			std::vector<SimEventRecord> allEvents = engine.getEventRecords();
			for(unsigned int i = 0; i < allEvents.size(); i++){
				if(events[allEvents[i].eventIx] == enterLunarRad){
					std::vector<double> state = traj.getStateByIx(allEvents[i].stepIx);
					double t = traj.getEpochByIx(allEvents[i].stepIx);

					// Convert to EM
					std::vector<double> stateEM = cr3bp_SE2EM_state(state, t,
						bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma(), emSys.getCharL(),
						emSys.getCharT(), seSys.getCharL(), seSys.getCharT(), seSys.getMu());

					double r[] = {stateEM[0], stateEM[1], stateEM[2]};		// Position in EM rot. coord.
					double rMag2 = sqrt(r[0]*r[0] + r[1]*r[1]);				// Magnitude of planar position components
					double theta = acos(r[0]/rMag2)*astrohelion::sign(r[1]);	// Angle between Moon and manifold
					dT = P_syn*theta/2/PI;									// Epoch shift to put encounter EXACTLY at the moon
					if(dT < 0)
						dT += P_syn;										// Can't leave earlier than t = 0, so add a month
				}
			}
		}

		// First, propagate to X = X_earth to leave the quasihalo
		Arcset_cr3bp o(&seSys);
		engine.clearEvents();
		engine.addEvent(cross_EarthPlane);
		engine.runSim(manICs[n].getStateByIx(0), 2*PI, &o);
		endEvents = engine.getEndEvents(&o);

		if(std::find(endEvents.begin(), endEvents.end(), cross_EarthPlane) == endEvents.end()){
			printErr("  Manifold %zu did not encounter the Earth-X plane...\n", n);
		}else{
			// If the manifold reached that plane, propagate to the XZ plane
			Arcset_cr3bp o2(&seSys);
			engine.clearEvents();
			engine.addEvent(cross_Y0);
			engine.runSim(o.getStateByIx(-1), PI, &o2);
			endEvents = engine.getEndEvents(&o2);
			if(std::find(endEvents.begin(), endEvents.end(), cross_Y0) == endEvents.end()){
				printErr("  Manifold %zu did not reach the XZ plane...\n", n);
			}else{
				// If the XZ plane has been reached, concatenate the two arcs and save
				o += o2;
				
				// Turn the manifold into a nodeset
				Arcset_cr3bp manData(o, numManNodes);

				// Create some nodes to represent the quasihalo
				Arcset_cr3bp haloData(o.getStateByIx(0), &seSys, -PI, numHaloNodes);
				double haloTOF = std::abs(haloData.getTotalTOF());

				// Tack the manifold on to the end
				haloData.appendSetAtNode(&manData, 0, 0, 0);

				// Transform to BCR4BP system
				double t0 = 0;	// t = 0 occurs at the beginning of the manifold
				Arcset_bc4bp bcNodes = bcr4bpr_SE2SEM(haloData, &bcSys, 0, t0);

				// Constrain the final node to be on the XZ_Plane
				double fixY0_Data[] = {NAN, 0, NAN, NAN, NAN, NAN};
				Constraint fixY0(Constraint_tp::STATE, bcNodes.getNumNodes()-1, fixY0_Data, 6);
				bcNodes.addConstraint(fixY0);

				// Correct to be continuous
				try{
					Arcset_bc4bp bcCorrected(&bcSys);
					corrector.setTol(9e-12);
					corrector.multShoot(&bcNodes, &bcCorrected);

					// Now propagate forward for a while from the final state
					Arcset_bc4bp natArc(&bcSys);
					// bcEngine.runSim(bcCorrected.getStatebyIx(-1), bcCorrected.getEpochByIx(-1), 360*24*3600/bcSys.getCharT(), &natArc);
					bcEngine.runSim(bcCorrected.getStateByIx(-1), bcCorrected.getEpochByIx(-1), 540*24*3600/bcSys.getCharT(), &natArc);
					Arcset_bc4bp fullTraj = Arcset_bc4bp::fromNodeset(bcCorrected);
					fullTraj += natArc;

					// Save the trajectory to file
					char name[32];
					snprintf(name, 32, "data/LPF_QH_4B_NaturalManifolds/Arcset%03zu_SEM.mat", n);
					fullTraj.saveToMat(name);

					Arcset_cr3bp fullTraj_SE = bcr4bpr_SEM2SE(fullTraj, &seSys);
					// snprintf(name, 32, "data/Arcset%03zu_SE.mat", n);
					// fullTraj_SE.saveToMat(name);

					Arcset_cr3bp fullTraj_EM = cr3bp_SE2EM(fullTraj_SE, &emSys, bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
					snprintf(name, 32, "data/LPF_QH_4B_NaturalManifolds/Arcset%03zu_EM.mat", n);
					fullTraj_EM.saveToMat(name);

					Arcset_cr3bp fullTraj_ECI = cr3bp_rot2inert(fullTraj_EM, 0);
					snprintf(name, 32, "data/LPF_QH_4B_NaturalManifolds/Arcset%03zu_ECI.mat", n);
					fullTraj_ECI.saveToMat(name);

					Arcset_cr3bp fullTraj_MCI = cr3bp_rot2inert(fullTraj_EM, 1);
					snprintf(name, 32, "data/LPF_QH_4B_NaturalManifolds/Arcset%03zu_MCI.mat", n);
					fullTraj_MCI.saveToMat(name);

				}catch(DivergeException &e){
					printErr("  Unable to correct manifold %03zu to be continuous in BC4BP... darn!\n", n);
				}
			}
		}
	}
}