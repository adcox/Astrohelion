#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

// Define low-thrust parameters
double f = 7e-2;		// nondim
double Isp = 1500;		// sec
double alpha0 = 0;		// rad

/**
 * @brief Construct an initial guess for a LTPO with a fixed alpha = 0 value 
 * from a natural orbit (e.g., a Lyapunov)
 * @details [long description]
 * 
 * @param famFile File in which the natural family is stored
 * @param C Use a family member with this Jacobi value to seed the initial guess
 * @param pSys Pointer to CR3BP-LT system data object
 * @param pLaw Pointer to CR3BP-LT control law
 * 
 * @return A CR3BP-LT arcset with discontinuities.
 */
Arcset_cr3bp_lt initGuessFromNat(const char *famFile, double C, 
	SysData_cr3bp_lt *pSys, ControlLaw_cr3bp_lt *pLaw){

	printf("Opening %s\n", famFile);
	SysData_cr3bp crSys(famFile);

	Family_PO_cr3bp lyapFam(&crSys);
	std::vector<ControlLaw*> loadedLaws {};
	lyapFam.readFromMat(famFile, loadedLaws);

	std::vector<Arcset_periodic> matches = lyapFam.getMemberByJacobi(C);
	if(matches.empty()){
		printErr("Could not locate any orbits at C = %f\n", C);
		return EXIT_SUCCESS;
	}

	// Use the first match
	Arcset_periodic lyap = matches[0];

	*pSys = SysData_cr3bp_lt(crSys.getPrimary(0).c_str(),
		crSys.getPrimary(1).c_str(), 1);		// Let m0 = 1
	std::vector<double> ltParams {sqrt(f), Isp}, ctrl0 {0, 0};
	*pLaw = ControlLaw_cr3bp_lt(ControlLaw_cr3bp_lt::CONST_MF_GENERAL, ltParams);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	// Fill the storage arcset with arcs propagated from the natural Lyapunov 
	// nodes with the specified low-thrust parameters
	double m = 1;
	std::vector<double> q0 {};
	Arcset_cr3bp_lt guess(pSys);
	for(unsigned int i = 0; i < lyap.getNumSegs(); i++){
		q0 = lyap.getStateByIx(i);
		q0.push_back(m);
		if(i == 0){
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &guess, pLaw);
		}else{
			Arcset_cr3bp_lt temp(pSys);
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &temp, pLaw);
			guess.appendSetAtNode(&temp, guess.getNodeByIx(-1).getID(),0,0);
		}

		m = guess.getStateByIx(-1).back();
	}

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}

	return guess;
}//====================================================

/**
 * @brief Get an initial guess from a family of LTPOs with a fixed alpha value
 * @details [long description]
 * 
 * @param ltpoFile [description]
 * @param H [description]
 * @param pSys [description]
 * @param laws [description]
 * @return [description]
 */
Arcset_cr3bp_lt initGuessFromLTFam_fixedAlpha(const char *ltpoFile, double H,
	SysData_cr3bp_lt *pSys, std::vector<ControlLaw *> &laws){

	*pSys = SysData_cr3bp_lt(ltpoFile);
	Family_PO_cr3bp_lt fam(pSys);
	fam.readFromMat(ltpoFile, laws);

	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(H);
	if(matches.empty()){
		printErr("Could not locate any orbits at Hlt = %f\n", H);

		for(auto p : laws){
			delete p;
			p = nullptr;
		}

		return EXIT_SUCCESS;
	}

	// Get the first guess
	Arcset_cr3bp_lt guess = static_cast<Arcset_cr3bp_lt>(matches[0]);

	// Set alpha0 to the alpha that corresponds to the family
	std::vector<double> ctrl = guess.getNodeByIx(0).getExtraParamVec(PARAMKEY_CTRL);
	if(ctrl.size() == 2){
		alpha0 = ctrl[0];
	}else{
		printErr("Expected ctrl at node 0 of family member to have size 2."
			" Size is %zu\n", ctrl.size());

		for(auto p : laws){
			delete p;
			p = nullptr;
		}

		return EXIT_SUCCESS;
	}
}//====================================================

/**
 * @brief Compute a family of LTPO all with the same thrust angle but with a 
 * variety of energies (low-thrust Hamiltonian)
 * 
 */
int main(){
	// Define law and system objects; overwritten by data from file
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_MF_GENERAL, 
		std::vector<double> {sqrt(f), Isp});
	SysData_cr3bp_lt ltSys("earth", "moon", 1);
	std::vector<ControlLaw *> loadedLaws {};

	// Load the L1 Lyapunov family
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp_earth-moon/L1_Lyap.mat");

	/*
	 *	Construct an initial guess from a natural orbit in a family at
	 *	a specified Jacobi Constant value.
	 */
	// Arcset_cr3bp_lt ltGuess = initGuessFromNat(filename, 3.1, &ltSys, &law);

	/*
	 *	Constrcut an initial guess from a low-thrust periodic orbit in a family 
	 *	with fixed alpha at a specified low-thrust Hamiltonian value.
	 */
	Arcset_cr3bp_lt ltGuess = initGuessFromLTFam_fixedAlpha(filename, -1.5,
		&ltSys, loadedLaws);
	for(unsigned int s = 0; s < ltGuess.getNumSegs(); s++){
		ltGuess.getSegRefByIx(s).setCtrlLaw(&law);
	}

	// Update control law to use the correct type/parameters
	law.setType(ControlLaw_cr3bp_lt::CONST_MF_GENERAL);
	law.setParams(std::vector<double> {sqrt(f), Isp});

	// Create parameters and other data storage objects
	// SysData_cr3bp_lt ltSys("earth", "moon", 1);	// Let M0 = 1
	std::vector<double> ctrl0 {0, 0},
		n0conData {NAN, NAN, 0, NAN, NAN, 0, 1},
		ctrlConData {1, 1};
	// ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_MF_GENERAL, ltParams);

	// Create Multiple shooting engine
	MultShootEngine shooter;
	shooter.setVerbosity(Verbosity_tp::SOME_MSG);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);

	// Be extra careful when converging the first orbit
	shooter.setMaxIts(200);
	shooter.setDoLineSearch(true);
	
	printf("Converging at alpha = %.2f deg\n", alpha0*180/PI);

	ctrl0[0] = alpha0;
	Arcset_cr3bp_lt ltConverged(&ltSys);

	// Periodicity Constraint
	int id0 = ltGuess.getNodeByIx(0).getID();
	double idf = static_cast<double>(ltGuess.getNodeByIx(-1).getID());

	std::vector<double> conData{ idf, idf, NAN, idf, NAN, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, id0, conData);

	ltGuess.addConstraint(periodicityCon);

	// Fix initial state components
	Constraint n0con(Constraint_tp::STATE, id0, n0conData);
	ltGuess.addConstraint(n0con);

	// Remove control from design vector; always fixed
	for(unsigned int i = 0; i < ltGuess.getNumNodes(); i++){
		Constraint con(Constraint_tp::RM_CTRL, ltGuess.getNodeByIx(i).getID(), 
			ctrlConData);
		ltGuess.addConstraint(con);
	}

	try{
		shooter.multShoot(&ltGuess, &ltConverged);
	}catch(DivergeException &e){
		printErr("Corrector diverged\n");
		return EXIT_SUCCESS;
	}

	ltConverged.putInChronoOrder();

	// Now, pull out the constraints and adjust them so that we can converge
	// orbits over a range of alpha but at the same H_lt value
	std::vector<Constraint> allCons = ltConverged.getAllConstraints();
	ltConverged.clearAllConstraints();
	for(unsigned int i = 0; i < allCons.size(); i++){
		if(allCons[i].getType() != Constraint_tp::RM_CTRL){
			ltConverged.addConstraint(allCons[i]);
		}
	}
	// Constraint control continuity along all segments
	// Thus, constraining control at initial node should propagate to all
	// later nodes/segments
	for(unsigned int i = 0; i < ltConverged.getNumSegs(); i++){
		Constraint con(Constraint_tp::CONT_CTRL,
			ltConverged.getSegRefByIx(i).getID(), ctrlConData);
		ltConverged.addConstraint(con);
	}

	// Get the low-thrust Hamiltonian value; overwrite q0
	std::vector<double> q0 = ltConverged.getSegRefByIx(0).getStateByRow(0);
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
	
	printColor(BOLDBLUE, "LT_CONVERGED\n");
	ltConverged.print();
	printf("H_lt = %f\n", H_lt);
	waitForUser();

	// Reset Shooter to "normal" behavior
	shooter.setMaxIts(20);
	shooter.setDoLineSearch(false);

	// Loop through alpha and converge orbits for all values, at same H_lt value
	double alphaStep = PI/180.0;
	Arcset_cr3bp_lt ltDiffAlpha(&ltSys), guess = ltConverged;
	std::vector<Arcset_cr3bp_lt> allArcs {ltConverged};
	for(double alpha = alpha0 + alphaStep; std::abs(alpha) < PI; 
		alpha += alphaStep){

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

			if(std::abs(ltDiffAlpha.getTotalTOF()) < 1e-2)
				throw DivergeException("TOF is too small");

			guess = ltDiffAlpha;	// Update initiial guess
			
			// printColor(BOLDBLUE, "GUESS\n");
			// guess.print();
			// printColor(BOLDBLUE, "LT_DIFF_ALPHA\n");
			// ltDiffAlpha.print();

			// waitForUser();
		}catch(DivergeException &e){
			printErr("Corrector diverged for alpha = %.2f deg\n", alpha*180/PI);
			// Go back to alpha0 and try opposite direction
			if(alphaStep > 0){
				alphaStep *= -1;
				alpha = alpha0;
				guess = ltConverged;
				continue;
			}else{
				break;
			}
		}

		// Insert arcs so they are in order w.r.t. alpha
		if(alphaStep > 0)
			allArcs.push_back(ltDiffAlpha);
		else
			allArcs.insert(allArcs.begin(), ltDiffAlpha);

		// Go back to alpha0 and reverse direction if you make it to PI
		if(alphaStep > 0 && std::abs(alpha + alphaStep) >= PI){
			alphaStep *= -1;
			alpha = alpha0;
			guess = ltConverged;
		}
	}

	// Put Arcs in a family structure
	Family_PO fam(&ltSys);
	for(unsigned int i = 0; i < allArcs.size(); i++)
		fam.addMember(allArcs[i]);

	fam.sortEigs();	// Sort those eigenvalues and eigenvectors

	sprintf(filename, "L1_Lyap_f%.1e_Hlt%.3f_law%d.mat", f,
		H_lt, law.getType());
	fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}