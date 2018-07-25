/*	This file builds a trajectory between a low-thrust L1 equilibrium point and 
 * 	a low-thrust L5 equilibrium point
 * 	
 * 	Initial guess begins at L1 with alpha = -60
 * 	Initial guess ends at L5 with alpha = -120
 *
 */
#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = ControlLaw_cr3bp_lt;

int main(int argc, char** argv){
	(void) argc;
	(void) argv;

	SysData_cr3bp natSys("earth", "moon");
	SysData_cr3bp_lt ltSys("earth", "moon", 1000);

	unsigned int lawID = ltlaw::GENERAL | ltlaw::CONST_F | ltlaw::CONST_M;
	double f = 7e-2;

	double alpha_unstable = -120*PI/180;
	double tof_unstable = 5.2605;
	std::vector<double> q0_unstable {1.154020, -0.026896, 0,
									-0.007312, 0.006925, 0, 1.0};

	double alpha_stable = -60*PI/180;
	// double tof_stable = 17.0326;
	double tof_stable = 9;
	std::vector<double> q0_stable {	0.232666, 0.531154, 0,
									-0.714302, 0.238309, 0, 1.0};

	double tof_spo = 6.5817;
	std::vector<double> q0_spo { 0.499343, -0.820068, 0,
								0.065382, 0.019300, 0, 1.0};

	ControlLaw_cr3bp_lt noLaw(ControlLaw::NO_CTRL, std::vector<double> {}),
		thrustLaw(lawID, std::vector<double>{f});

	SimEngine sim;

	Arcset_cr3bp_lt spo(&ltSys), stable(&ltSys), unstable(&ltSys);

	sim.runSim_manyNodes(q0_spo, std::vector<double>{}, 0, tof_spo, 4, 
		&spo, &noLaw);

	sim.runSim_manyNodes(q0_unstable, std::vector<double>{alpha_unstable, 0}, 
		0, tof_unstable, std::floor(tof_unstable), &unstable, &thrustLaw);

	sim.runSim_manyNodes(q0_stable, std::vector<double>{alpha_stable, 0}, 
		0, tof_stable, 5, &stable, &thrustLaw);

	// Concatenate Arcs into a single arcset
	unstable.appendSetAtNode(&stable, unstable.getNodeRefByIx(-1).getID(), 
		stable.getNodeRefByIx(0).getID(), 0);
	unstable.appendSetAtNode(&spo, unstable.getNodeRefByIx(-1).getID(),
		spo.getNodeRefByIx(0).getID(), 0);

	unstable.sortChrono();
	unstable.print();
	unstable.saveToMat("temp_transfer.mat");

	// Create constraints on the initial position and energy
	std::vector<double> q0 = unstable.getStateByIx(0);
	std::vector<double> initState {q0[0], q0[1], NAN, NAN, NAN, NAN, 1};
	Constraint conInitPos(Constraint_tp::STATE, unstable.getNodeRefByIx(0).getID(),
		initState);
	Constraint conInitJacobi(Constraint_tp::JC, unstable.getNodeRefByIx(0).getID(),
		std::vector<double> {unstable.getJacobiByIx(0)});

	// Constrain the nearly-final SPO state
	std::vector<double> finalState = unstable.getStateByIx(-3);
	finalState.erase(finalState.begin()+6, finalState.end());	// get rid of everything after six states
	Constraint conFinalState(Constraint_tp::STATE, unstable.getNodeRefByIx(-3).getID(),
		finalState);

	unstable.addConstraint(conInitPos);
	unstable.addConstraint(conInitJacobi);
	unstable.addConstraint(conFinalState);

	MultShootEngine shooter;
	// shooter.setDoLineSearch(true);
	// shooter.setMaxIts(50);
	// shooter.setIgnoreDiverge(true);
	Arcset_cr3bp_lt converged(&ltSys);
	try{
		shooter.multShoot(&unstable, &converged);
	}catch(const Exception &e){
		printErr("%s\n", e.what());
	}
	converged.saveToMat("converged_constM.mat");

	// Change law to variable mass
	lawID = ltlaw::GENERAL | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	thrustLaw.setParams(std::vector<double> {f, 3000});
	thrustLaw.setType(lawID);

	// Add a constraint to fix the initial mass
	// Constraint conInitM(Constraint_tp::STATE, converged.getNodeRefByIx(0).getID(),
	// 	std::vector<double> {NAN, NAN, NAN, NAN, NAN, NAN, 1});
	// converged.addConstraint(conInitM);

	converged.print();
	waitForUser();
	
	Arcset_cr3bp_lt converged2(&ltSys);
	try{
		shooter.multShoot(&converged, &converged2);
	}catch(const Exception &e){
		printErr("%s\n", e.what());
	}
	converged2.saveToMat("converged_varM.mat");

	return EXIT_SUCCESS;
}//====================================================