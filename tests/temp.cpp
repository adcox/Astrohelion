#include "AllIncludes.hpp"

using namespace astrohelion;

void computeMap(const std::vector<double> ICs, const int numReturns, const double C, const SysData_cr3bp_lt* pSys,
	ControlLaw_cr3bp_lt *pLaw, char manType){

	unsigned int stateSize = pSys->getDynamicsModel()->getCoreStateSize();
	int numICs = ICs.size()/stateSize;

	printf("Beginning map creation for C = %.2f, Manifold type = %c, ", C, manType);
	printf("Law = %u, Isp = %.0f, F = %.2f mN, M0 = %.2f kg\n", pLaw->getLawType(), pLaw->getIsp(), pLaw->getThrust_dim()*1000, pSys->getRefMass());
	printf("  Will integrate %d orbits, %d total integrations\n", numICs, 
		numICs*numReturns);

	// Create the simulation engine
	SimEngine engine;
	engine.setVerbosity(Verbosity_tp::NO_MSG);
	engine.setMakeDefaultEvents(false);
	// engine.setVerbosity(Verbosity_tp::SOME_MSG);
	engine.setMaxCompTime(30);	// No more than 30 seconds per initial condition

	// Create an event that fires at the map crossing
	Event mapCross(Event_tp::XZ_PLANE, 1, true);
	mapCross.setStopCount(numReturns);
	engine.addEvent(mapCross);

	if(manType == 'U'){
		// Unstable manifolds propagate forward in time
		std::vector<double> minMassData {0.75};
		Event minMass(Event_tp::MASS, 0, true, minMassData);
		engine.addEvent(minMass);
	}else{
		// Unstable manifolds propagate backward in time
		std::vector<double> maxMassData {1.0};
		Event maxMass(Event_tp::MASS, 0, true, maxMassData);
		engine.addEvent(maxMass);
		engine.setRevTime(true);
	}

	// Create another event that fires if the arc strays too far from P1
	std::vector<double> maxDistData {0, 4.0};
	Event maxDist(Event_tp::DIST, 1, true, maxDistData);
	engine.addEvent(maxDist);

	// Create vector to store all map returns
	std::vector<double> returns, times;
	returns.reserve(numReturns*numICs*stateSize);		// Allocate space assuming all returns are reached
	times.reserve(numReturns*numICs);

	std::cout << "  Beginning computation loop" << std::endl;

	double startTime = omp_get_wtime();
	unsigned int count = 0;

	#pragma omp parallel for firstprivate(engine, mapCross, returns, times) schedule(dynamic)
	for(int n = 0; n < numICs; n++){
		char filename[128];
		sprintf(filename, "data_%c_law%u/JC-%08d_M0-%.0fkg_F-%.0fmN_Isp-%.0fsec_Orb%08d.mat",
			manType, pLaw->getLawType(), (int)floor(C*1e7), pSys->getRefMass(), pLaw->getThrust_dim()*1000, pLaw->getIsp(), n);

		if(!fileExists(filename)){
			std::vector<double> IC(ICs.begin()+n*stateSize, ICs.begin()+(n+1)*stateSize);

			Arcset_cr3bp_lt traj(pSys);

			try{
				engine.runSim(IC, numReturns*6*PI, &traj, pLaw);	// Run for up to 6*pi for each map return
			}catch(DivergeException &e){
				// printErr("Integrator failed (GSL): saving what data we have...\n");
			}catch(const Exception &e){
				// printErr("Error:\n%s\n", e.what());
				
				returns.clear();	// Reset some things
				times.clear();
				continue; // skip to next loop iteration
			}
			
			for(unsigned int i = 0; i < traj.getNumNodes(); i++){
				if(traj.getNodeRefByIx_const(i).getTriggerEvent() == mapCross.getType()){
					std::vector<double> state = traj.getStateByIx(i);
					returns.insert(returns.end(), state.begin(), state.end());
					times.push_back(traj.getTimeByIx(i));
				}
			}

			// printf("Manifold %d ended with %s event\n", n, Event::getEventTpStr(traj.getNodeByIx(-1).getTriggerEvent()));
			
			// Save data to file
			mat_t *matfp = Mat_CreateVer(filename, nullptr, MAT_FT_DEFAULT);
			saveMatrixToFile(matfp, "data", returns, returns.size()/stateSize, stateSize);
			saveMatrixToFile(matfp, "times", times, times.size(), 1);
			saveDoubleToFile(matfp, "IC_ix", static_cast<double>(n));
			saveDoubleToFile(matfp, "Thrust", pLaw->getThrust());
			saveDoubleToFile(matfp, "Isp", pLaw->getIsp());
			saveDoubleToFile(matfp, "EndEventType", static_cast<double>(traj.getNodeRefByIx_const(-1).getTriggerEvent()));
			pSys->saveToMat(matfp);	// save all relevant parameters
			Mat_Close(matfp);

			// Reset some things
			returns.clear();
			times.clear();

			#pragma omp atomic
			++count;

			if(count % 100 == 0){
				#pragma omp critical
				printf("  Progress: %d / %d (%.2f%%)", count, numICs, 100*(static_cast<double>(count)/static_cast<double>(numICs)));
				printf(" - Using %d / %d threads\n", omp_get_num_threads(), omp_get_max_threads());
				fflush(stdout);
			}
		}else{
			printf("%s already exists... skipping\n", filename);
		}// end of fileExists() check
	}// end of loop through ICs

	double endTime = omp_get_wtime();

	if(endTime - startTime > 60)
		printf("Done. Total time = %.2f min\n", (endTime - startTime)/60);
	else
		printf("Done. Total time = %.2f sec\n", endTime - startTime);

	printf("--------------------------------------\n\n");
}//============================================================

void propManifolds(Arcset_cr3bp_lt *pOrb, ControlLaw_cr3bp_lt *pLaw, int numArcs, int numReturns){
	const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt *>(pOrb->getSysData());
	Arcset_cr3bp_lt discreteOrb(pSys);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(pOrb->getStateByIx(0), pOrb->getTotalTOF(), numArcs + 1, &discreteOrb, pLaw);
	discreteOrb.saveToMat("data/discrete.mat");

	ManifoldEngine many;
	std::vector<Arcset_cr3bp_lt> traj_u = many.computeSetFromLTPeriodic(Manifold_tp::MAN_U, &discreteOrb, pLaw, numArcs, 0);
	std::vector<Arcset_cr3bp_lt> traj_s = many.computeSetFromLTPeriodic(Manifold_tp::MAN_S, &discreteOrb, pLaw, numArcs, 0);

	waitForUser();

	std::vector<double> ics_u, ics_s;
	for(unsigned int n = 0; n < traj_u.size(); n++){
		std::vector<double> state = traj_u[n].getStateByIx(0);
		ics_u.insert(ics_u.end(), state.begin(), state.end());
		ics_u.back() = 1.0;			// Unstable arcs begin with 100% mass

		// printf("Unstable IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n", state[0], state[1],
		// 	state[2], state[3], state[4], state[5], 1.0);

		state = traj_s[n].getStateByIx(0);
		ics_s.insert(ics_s.end(), state.begin(), state.end());
		ics_s.back() = 0.75;		// Stable arcs end with 75% mass (propagate in reverse time)
	}

	
	computeMap(ics_u, numReturns, pOrb->getJacobiByIx(0), pSys, pLaw, 'U');
	computeMap(ics_s, numReturns, pOrb->getJacobiByIx(0), pSys, pLaw, 'S');
}//====================================================

int main(){
	double C = 3.15;
	double M0 = 1000;
	double F = 1e-16;		// set F to a very small number to allow nontrivial conversion between control laws
	double f = 1e-2;
	double Isp = 1500;

	SysData_cr3bp_lt sys("earth", "moon", M0);
	Fam_cr3bp fam("../data/families/cr3bp_earth-moon/L1_Lyap.mat");

	std::vector<FamMember_cr3bp> matches = fam.getMemberByJacobi(C);

	if(matches.size() == 0){
		printErr("Could not find/correct any orbits!\n");
		return 0;
	}

	ControlLaw_cr3bp_lt lawNat(ControlLaw::NO_CTRL, 0, 1);
	ControlLaw_cr3bp_lt lawLeft(ControlLaw_cr3bp_lt::CONST_C_2D_LEFT, F, Isp);
	ControlLaw_cr3bp_lt lawRight(ControlLaw_cr3bp_lt::CONST_C_2D_RIGHT, F, Isp);
	ControlLaw_cr3bp_lt lawGeneral(ControlLaw_cr3bp_lt::GENERAL_CONST_F, F, Isp);
	
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	Arcset_cr3bp_lt natPO(&sys), leftLTPO(&sys), rightLTPO(&sys);
	std::vector<double> ic = matches[0].getIC();
	ic.push_back(1);

	// First run halfway around
	sim.runSim_manyNodes(ic, matches[0].getTOF()/2, 4, &natPO, &lawRight);
	ic = natPO.getStateByIx(-1);
	ic[6] = 1;	// reset mass

	// propagate from other side
	natPO = Arcset_cr3bp_lt(&sys);
	sim.runSim_manyNodes(ic, matches[0].getTOF(), 7, &natPO, &lawNat);
	sim.runSim_manyNodes(ic, matches[0].getTOF(), 7, &rightLTPO, &lawRight);
	sim.runSim_manyNodes(ic, matches[0].getTOF(), 7, &leftLTPO, &lawLeft);

	double n = static_cast<double>(natPO.getNodeByIx(-1).getID());
	std::vector<double> periodicityData {n, n, NAN, n, n, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, 0, periodicityData);

	Constraint fixJC(Constraint_tp::JC, 0, &C, 1);

	std::vector<double> initState {NAN, 0, NAN, NAN, NAN, NAN, 1};
	Constraint fixInitState(Constraint_tp::STATE, 0, initState);

	natPO.addConstraint(periodicityCon);
	natPO.addConstraint(fixJC);
	natPO.addConstraint(fixInitState);

	rightLTPO.addConstraint(periodicityCon);
	rightLTPO.addConstraint(fixJC);
	rightLTPO.addConstraint(fixInitState);

	leftLTPO.addConstraint(periodicityCon);
	leftLTPO.addConstraint(fixJC);
	leftLTPO.addConstraint(fixInitState);

	lawLeft.setThrust_nondim(f, &sys);
	lawRight.setThrust_nondim(f, &sys);

	Arcset_cr3bp_lt correctedNat(&sys), correctedLeft(&sys), correctedRight(&sys);

	MultShootEngine shooter;
	// shooter.setDoLineSearch(true);
	// shooter.setMaxIts(200);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);
	shooter.multShoot(&natPO, &correctedNat);
	shooter.multShoot(&leftLTPO, &correctedLeft);
	shooter.multShoot(&rightLTPO, &correctedRight);

	correctedNat.saveToMat("data/EM_L1_Lyap_nat.mat");
	correctedLeft.saveToMat("data/EM_L1_Lyap_left.mat");
	correctedRight.saveToMat("data/EM_L1_Lyap_right.mat");
	
	propManifolds(&correctedNat, &lawNat, 50, 5);
	propManifolds(&correctedLeft, &lawLeft, 50, 5);
	propManifolds(&correctedRight, &lawRight, 50, 5);
}//====================================================





