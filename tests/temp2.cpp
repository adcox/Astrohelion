#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(){

	// Parameters to isolate the equilibria
	double f = 7e-2;
	double alpha = 55.0*PI/180.0;
	SysData_cr3bp_lt sys("earth", "moon", 1);
	std::vector<double> params {};
	ControlLaw_cr3bp_lt law(ControlLaw::NO_CTRL, params);	// will be overwritten

	// For f = 7e-2, L3, L4, and L5 are all the same ZAC
	std::vector<double> LPts;
	DynamicsModel_cr3bp_lt::getEquilibPt(&sys, 3, f, 1e-6, &LPts);

	// Find Eq pt closest to the desired angle
	double minDiff = 4, dS = 0, dC = 0;
	double eqPt[3] = {0};

	// Find the equilibrium solution that is closes to the desired angle
	// and also above the x-axis (Looking for L4-like points)
	for(unsigned int i = 0; i < LPts.size()/3; i++){
		dS = std::abs(sin(LPts[3*i]) - sin(alpha));
		dC = std::abs(cos(LPts[3*i]) - cos(alpha));
		if(dS + dC < minDiff && LPts[3*i+2] > 0.6){
			minDiff = dS + dC;
			eqPt[0] = LPts[3*i+1];
			eqPt[1] = LPts[3*i+2];
		}
	}

	if(minDiff == 4){
		throw Exception("Did not find eq pt near desired angle");
	}

	LinMotionEngine_cr3bp_lt linEngine;
	linEngine.setVerbosity(Verbosity_tp::ALL_MSG);
	Arcset_cr3bp_lt linArc(&sys);

	double x0[] = {0.0075, 0, 0};	// Step off of equilibria
	linEngine.getLinear(eqPt, sqrt(f), alpha, x0, LinMotion_tp::OSC,
		&linArc, &law, 5);

	law.setType(ControlLaw_cr3bp_lt::CONST_MF_GENERAL);

	linArc.saveToMat("linArc.mat");

	MultShootEngine msEngine;
	msEngine.setVerbosity(Verbosity_tp::SOME_MSG);
	msEngine.setDoLineSearch(true);
	msEngine.setMaxIts(200);
	msEngine.setTOFType(MSTOF_tp::VAR_EQUALARC);

	// Periodicity Constraint
	double nf = linArc.getNumNodes() - 1;
	std::vector<double> perConData {nf, nf, NAN, NAN, nf, NAN};
	Constraint perCon(Constraint_tp::MATCH_CUST, 0, perConData);
	linArc.addConstraint(perCon);

	// Fix initial Mass
	std::vector<double> stateConData {NAN, NAN, NAN, 0, NAN, NAN, 1};
	Constraint stateCon(Constraint_tp::STATE, 0, stateConData);
	linArc.addConstraint(stateCon);

	Arcset_cr3bp_lt nonlinArc(&sys), nonlinArc2(&sys), nonlinArc3(&sys);

	try{
		msEngine.multShoot(&linArc, &nonlinArc);

		// nonlinArc.print();
		// char filename[128];
		// sprintf(filename, "guess_L4_a%06.2f.mat", alpha*180/PI);
		// nonlinArc.saveToMat(filename);
	}catch(DivergeException &e){
		printErr("Failed to converge: %s\n", e.what());
		return EXIT_SUCCESS;
	}

	// Propagate to xdot = 0
	std::vector<double> stateEvtData {3, 0};
	Event stateEvt(Event_tp::STATE_PLANE, 0, true, stateEvtData);
	SimEngine sim;
	// sim.setVerbosity(Verbosity_tp::ALL_MSG);
	sim.addEvent(stateEvt);
	Arcset_cr3bp_lt arc(&sys);
	bool foundPt = false;
	for(unsigned int n = 0; n < nonlinArc.getNumNodes(); n++){
		arc.reset();
		std::vector<double> q0 = nonlinArc.getStateByIx(n);
		std::vector<double> ctrl0 = nonlinArc.getNodeRefByIx(n).getExtraParamVec(PARAMKEY_CTRL);
		double tof = nonlinArc.getTOFByIx(n);
		sim.runSim(q0, ctrl0, 0, tof, &arc, &law);

		if(arc.getNodeByIx(-1).getTriggerEvent() == Event_tp::STATE_PLANE){
			// Great! Use this state as a new initial state for an arcset
			foundPt = true;
			q0 = arc.getStateByIx(-1);
			nonlinArc.reset();
			sim.clearEvents();
			sim.runSim_manyNodes(q0, ctrl0, 0, linArc.getTotalTOF(), 5, 
				&nonlinArc, &law);
			break;
		}
	}

	if(!foundPt){
		printErr("Did not find dx/dt = 0\n");
		return EXIT_SUCCESS;
	}

	// Update periodicity constraint
	nf = nonlinArc.getNodeByIx(-1).getID();
	perConData = std::vector<double> {nf, nf, NAN, nf, nf, NAN};
	perCon.setData(perConData);

	// Update state consraint
	stateConData = std::vector<double> {NAN, NAN, NAN, NAN, NAN, NAN, 1};
	stateCon.setData(stateConData);

	nonlinArc.addConstraint(perCon);
	nonlinArc.addConstraint(stateCon);

	// Correct again for periodicity
	try{
		msEngine.multShoot(&nonlinArc, &nonlinArc2);
		nonlinArc2.saveToMat("nonlin2.mat");
	}catch(DivergeException &e){
		printErr("Failed to converge: %s\n", e.what());
		return EXIT_SUCCESS;
	}
	// nonlinArc2.print();

	// Remove control from the free variables; it remains fixed
	std::vector<double> ctrl0 {alpha, 0}, ctrlCont {1, 1};
	Constraint fixCtrlCon(Constraint_tp::CTRL, 
		nonlinArc2.getNodeByIx(0).getID(), ctrl0);
	nonlinArc2.addConstraint(fixCtrlCon);
	for(unsigned int i = 0; i < nonlinArc2.getNumSegs(); i++){
		Constraint contCon(Constraint_tp::CONT_CTRL, 
			nonlinArc2.getSegByIx(i).getID(), ctrlCont);
		nonlinArc2.addConstraint(contCon);
	}

	// nonlinArc2.print();

	// Correct Again
	// msEngine.setVerbosity(Verbosity_tp::ALL_MSG);
	try{
		msEngine.multShoot(&nonlinArc2, &nonlinArc3);

		char filename[128];
		sprintf(filename, "guess_L4_a%06.2f.mat", alpha*180/PI);
		nonlinArc3.saveToMat(filename);
	}catch(DivergeException &e){
		printErr("Failed to converge: %s\n", e.what());
		return EXIT_SUCCESS;
	}

	return EXIT_SUCCESS;
}