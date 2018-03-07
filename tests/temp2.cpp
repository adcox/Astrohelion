#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(){
	// Load the L1 Lyapunov family
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp_earth-moon/L1_Lyap.mat");
	printf("Opening %s\n", filename);
	SysData_cr3bp crSys(filename);

	Family_PO lyapFam(&crSys);
	std::vector<ControlLaw*> loadedLaws {};
	lyapFam.readFromMat(filename, loadedLaws);

	unsigned int ix0 = 10;

	Arcset_periodic lyap = lyapFam.getMember(ix0);

	// Define low-thrust parameters
	double f = 7e-2;
	double Isp = 1500;

	// Create parameters and other data storage objects
	SysData_cr3bp_lt ltSys("earth", "moon", 1);	// Let M0 = 1
	std::vector<double> ltParams {sqrt(f), Isp},
		ctrl0 {0, 0},
		n0conData {NAN, NAN, 0, NAN, NAN, 0, 1},
		ctrlConData {1, 1};
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_MF_GENERAL, ltParams);

	// Create sim engine and storage arcset
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt ltGuess(&ltSys), ltConverged(&ltSys);

	// Create Multiple shooting engine
	MultShootEngine shooter;
	shooter.setVerbosity(Verbosity_tp::ALL_MSG);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);
	Family_PO fam(&ltSys);

	double alpha = 0;
	
	printf("Converging at alpha = %.2f deg\n", alpha*180/PI);

	ctrl0[0] = alpha;
	ltGuess.reset();
	ltConverged.reset();

	// Fill the storage arcset with arcs propagated from the natural Lyapunov nodes
	// with the specified low-thrust parameters
	double m = 1;
	std::vector<double> q0 {};
	for(unsigned int i = 0; i < lyap.getNumSegs(); i++){
		q0 = lyap.getStateByIx(i);
		q0.push_back(m);
		if(i == 0){
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &ltGuess, &law);
			
		}else{
			Arcset_cr3bp_lt temp(&ltSys);
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &temp, &law);
			ltGuess.appendSetAtNode(&temp, ltGuess.getNodeByIx(-1).getID(), 0, 0);
		}

		m = ltGuess.getStateByIx(-1).back();
	}

	// Periodicity Constraint
	int id0 = ltGuess.getNodeByIx(0).getID();
	int idf = ltGuess.getNodeByIx(-1).getID();

	std::vector<double> conData{ idf, idf, NAN, idf, NAN, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, id0, conData);

	ltGuess.addConstraint(periodicityCon);

	// Fix initial state components
	Constraint n0con(Constraint_tp::STATE, ltGuess.getNodeByIx(0).getID(), n0conData);
	ltGuess.addConstraint(n0con);

	// Remove control from design vector; always fixed
	for(unsigned int i = 0; i < ltGuess.getNumNodes(); i++){
		Constraint con(Constraint_tp::RM_CTRL, ltGuess.getNodeByIx(i).getID(), ctrlConData);
		ltGuess.addConstraint(con);
	}

	try{
		shooter.multShoot(&ltGuess, &ltConverged);
	}catch(DivergeException &e){
		printErr("Corrector diverged\n");
	}

	ltConverged.putInChronoOrder();

	// ltGuess.print();
	// ltConverged.print();

	// Now, pull out the constraints and adjust them so that we can converge
	// orbits over a range of alpha but at the same H_lt value
	std::vector<Constraint> allCons = ltConverged.getAllConstraints();
	ltConverged.clearAllConstraints();
	for(unsigned int i = 0; i < allCons.size(); i++){
		if(allCons[i].getType() == Constraint_tp::RM_CTRL){
			if(allCons[i].getID() != ltConverged.getNodeByIx(0).getID()){
				ltConverged.addConstraint(allCons[i]);
			}
		}else{
			ltConverged.addConstraint(allCons[i]);
		}
	}

	// Get the low-thrust Hamiltonian value; overwrite q0
	q0 = ltConverged.getSegRefByIx(0).getStateByRow(0);
	double H_lt = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q0[0]), &ltSys,
		&law);


	// Create a constraint to fix the low-thrust Hamiltonian
	Constraint HltCon(Constraint_tp::HLT, ltConverged.getNodeByIx(0).getID(),
		&H_lt, 1);
	ltConverged.addConstraint(HltCon);

	// Overwrite allCons with the constraints we want to apply to every orbit
	allCons = ltConverged.getAllConstraints();

	// Create a new constraint to fix the thrust angle
	Constraint ctrl0Con(Constraint_tp::CTRL, ltConverged.getNodeByIx(0).getID(),
		ctrl0);
	ltConverged.addConstraint(ctrl0Con);	// only add so it shows up in file
	
	// Add the converged member to the family
	fam.addMember(ltConverged);

	printColor(BOLDBLUE, "LT_CONVERGED\n");
	ltConverged.print();

	// Loop through alpha and converge orbits for all values, at same H_lt value
	double alphaStep = PI/180.0;
	Arcset_cr3bp_lt ltDiffAlpha(&ltSys), guess = ltConverged;
	for(alpha = alphaStep; std::abs(alpha) < PI; alpha += alphaStep){
		ctrl0[0] = alpha;
		ctrl0Con.setData(ctrl0);

		guess.clearAllConstraints();
		for(unsigned int i = 0; i < allCons.size(); i++){
			guess.addConstraint(allCons[i]);
		}
		guess.addConstraint(ctrl0Con);

		ltDiffAlpha.reset();
		try{
			printf("Converging at alpha = %.2f deg\n", alpha*180/PI);
			shooter.multShoot(&guess, &ltDiffAlpha);

			q0 = ltDiffAlpha.getSegRefByIx(0).getStateByRow(0);
			double H = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q0[0]),
				&ltSys, &law);
			printf("Err in H_lt = %.2e\n", std::abs(H - H_lt));

			guess = ltDiffAlpha;	// Update initiial guess
			
			printColor(BOLDBLUE, "GUESS\n");
			guess.print();
			printColor(BOLDBLUE, "LT_DIFF_ALPHA\n");
			ltDiffAlpha.print();

			waitForUser();
		}catch(DivergeException &e){
			printErr("Corrector diverged for alpha = %.2f deg\n", alpha*180/PI);
			if(alphaStep > 0){
				alphaStep *= -1;
				alpha = 0;
				guess = ltConverged;
				continue;
			}else{
				break;
			}
		}
		fam.addMember(ltDiffAlpha);

		// Go back to zero and reverse direction if you make it to PI
		if(alphaStep > 0 && std::abs(alpha + alphaStep) >= PI){
			alphaStep *= -1;
			alpha = 0;
			guess = ltConverged;
		}
	}

	fam.sortEigs();

	sprintf(filename, "L1_Lyap_f%.1e_Hlt_%.3f_law%d.mat", f,
		H_lt, law.getType());
	fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}