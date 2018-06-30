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

	double alpha_stable = -60*PI/180;
	// double tof_stable = 17.6933;
	// std::vector<double> q0_stable {	-0.035488, -0.265677, 0,
	// 								2.057210, -0.418170, 0, 1.0};
	double tof_stable = 15.9956;
	std::vector<double> q0_stable {	-0.621271, -0.003300, 0, -0.095895, -0.683197, 0, 1.0};

	double alpha_unstable = 96*PI/180;
	// double tof_unstable = 5.0333;
	// std::vector<double> q0_unstable {0.835832, -0.012246, 0,
	// 								-0.014001, 0.005538, 0, 1.0};
	double tof_unstable = 3.3246;
	std::vector<double> q0_unstable {0.833916, 0.019025, 0, -0.012808, 0.007409, 0, 1.0};

	double tof_spo = 6.5817;
	std::vector<double> q0_spo { 0.499343, -0.820068, 0.000000, 0.065382, 0.019300, 0, 1.0};

	ControlLaw_cr3bp_lt noLaw(ControlLaw::NO_CTRL, std::vector<double> {}),
		ltLaw(lawID, std::vector<double>{f});

	SimEngine sim;

	Arcset_cr3bp_lt spo(&ltSys), stable(&ltSys), unstable(&ltSys);

	sim.runSim_manyNodes(q0_spo, std::vector<double>{}, 0, tof_spo, 4, 
		&spo, &noLaw);

	sim.runSim_manyNodes(q0_stable, std::vector<double>{alpha_stable, 0}, 
		0, std::floor(tof_stable), tof_stable, &stable, &ltLaw);

	sim.runSim_manyNodes(q0_unstable, std::vector<double>{alpha_unstable, 0}, 
		0, std::floor(tof_unstable), tof_unstable, &unstable, &ltLaw);

	// Concatenate Arcs into a single arcset
	unstable.appendSetAtNode(&stable, unstable.getNodeRefByIx(-1).getID(), 
		stable.getNodeRefByIx(0).getID(), 0);
	unstable.appendSetAtNode(&spo, unstable.getNodeRefByIx(-1).getID(),
		spo.getNodeRefByIx(0).getID(), 0);

	unstable.sortChrono();
	unstable.print();
	unstable.saveToMat("temp_transfer.mat");

	// Create constraints on the initial position and energy
	std::vector<double> initState {q0_unstable[0], q0_unstable[1], NAN,
									NAN, NAN, NAN, NAN};
	Constraint conInitPos(Constraint_tp::STATE, unstable.getNodeRefByIx(0).getID(),
		initState);
	Constraint conInitJacobi(Constraint_tp::JC, unstable.getNodeRefByIx(0).getID(),
		std::vector<double> {unstable.getJacobiByIx(0)});

	// Constrain the nearly-final SPO state
	std::vector<double> finalState = spo.getStateByIx(-1);
	finalState.erase(finalState.begin()+6, finalState.end());	// get rid of everything after six states
	Constraint conFinalState(Constraint_tp::STATE, unstable.getNodeRefByIx(-3).getID(),
		finalState);

	unstable.addConstraint(conInitPos);
	unstable.addConstraint(conInitJacobi);
	unstable.addConstraint(conFinalState);

	MultShootEngine shooter;
	Arcset_cr3bp_lt converged(&ltSys);
	try{
		shooter.multShoot(&unstable, &converged);
	}catch(const Exception &e){
		printErr("%s\n", e.what());
	}
	converged.saveToMat("converged_constM.mat");

	// Change law to variable mass
	lawID = ltlaw::GENERAL | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	ltLaw.setParams(std::vector<double> {f, 10000});
	ltLaw.setType(lawID);

	// Add a constraint to fix the initial mass
	Constraint conInitM(Constraint_tp::STATE, converged.getNodeRefByIx(0).getID(),
		std::vector<double> {NAN, NAN, NAN, NAN, NAN, NAN, 1});
	converged.addConstraint(conInitM);

	converged.print();
	waitForUser();
	
	Arcset_cr3bp_lt converged2(&ltSys);
	shooter.setDoLineSearch(true);
	shooter.setMaxIts(200);
	try{
		shooter.multShoot(&converged, &converged2);
	}catch(const Exception &e){
		printErr("%s\n", e.what());
	}
	converged2.saveToMat("converged_varM.mat");

	return EXIT_SUCCESS;
}//====================================================