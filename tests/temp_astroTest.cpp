#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char** argv){
	(void) argc;
	(void) argv;

	double f = 1e-2;
	SysData_cr3bp_lt sys("earth", "moon", 100);
	ControlLaw_cr3bp_lt law(ControlLaw::NO_CTRL, f, 3000);

	// Compute Equilibrium solutions beginning with analytic or semi-analytic expressions for 
	// the solutions at alpha = 0
	std::vector<double> L1, L2, L3;
	DynamicsModel_cr3bp_lt::getEquilibPt(&sys, 1, law.getParam(0), 1e-6, &L1);
	DynamicsModel_cr3bp_lt::getEquilibPt(&sys, 2, law.getParam(0), 1e-6, &L2);
	DynamicsModel_cr3bp_lt::getEquilibPt(&sys, 3, law.getParam(0), 1e-6, &L3);

	// Indices in center-center and center-saddle regions at f = 1e-2
	int L3_center_pts[] = {433, 593, 672};
	int L3_saddle_pts[] = {40, 386, 755};
	int L2_saddle_pts[] = {0, 74, 149, 234};
	int L1_saddle_pts[] = {0, 74, 149, 230};

	// -------------------------------
	// Construct Linear Approximation
	// -------------------------------
	unsigned int i = 0;				// Which of the L3_center_pts to use
	unsigned int numNodes = 5;
	double x0[3] = {0.005, 0, 0};	// Offset from equilibria
	double eqPt[3] = {0};
	std::copy(L3.begin() + 3*L3_center_pts[i]+1, L3.begin() + 3*L3_center_pts[i]+3, eqPt);
	double alpha = L3[3*L3_center_pts[i] + 0];

	LinMotionEngine_cr3bp_lt linEngine;
	Arcset_cr3bp_lt linArc(&sys);

	linEngine.getLinear(eqPt, f, alpha, x0, LinMotion_tp::OSC, &linArc, &law, numNodes);

	std::vector<double> perConData {0, 0, NAN, 0, 0, NAN, NAN};
	Constraint perCon(Constraint_tp::MATCH_CUST, linArc.getNodeRefByIx(-1).getID(), perConData);

	std::vector<double> stateConData {NAN, eqPt[1], 0, NAN, NAN, 0, 1};
	Constraint stateCon(Constraint_tp::STATE, linArc.getNodeRefByIx(0).getID(), stateConData);

	linArc.addConstraint(perCon);
	linArc.addConstraint(stateCon);

	// -------------------------------
	// Natural Parameter Continuation
	// -------------------------------

	NatParamEngine npe;
	npe.setTol(9e-12);
	// npe.setNumOrbits(25);

	Arcset_cr3bp_lt a1(&sys), a2(&sys);
	std::vector<Arcset> allArcs {};
	std::vector<double> alwaysFixed {};
	std::vector<unsigned int> indVarIx {0, 4};
	std::vector<unsigned int> depVarIx {4};
	try{
		npe.continuePO(&linArc, &a1, &a2, allArcs, alwaysFixed, indVarIx, depVarIx);
	}catch(Exception &e){
		printErr("Exception: %s\n", e.what());
	}

	Family_PO fam(&sys);
	for(unsigned a = 0; a < allArcs.size(); a++){
		Arcset_cr3bp_lt temp(allArcs[a]);
		fam.addMember(temp);
	}

	// fam.sortEigs();
	fam.setName("Low-Thrust Periodic");
	fam.saveToMat("../../data/families_toBeChecked/familyTest_LT.mat");
	
	return EXIT_SUCCESS;
}