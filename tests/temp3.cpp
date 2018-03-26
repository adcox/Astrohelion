#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

// Define low-thrust parameters
double f = 7e-2;		// nondim
double Isp = 1500;		// sec

// Extent of family: move up to this much from initial H_lt value
double H_range = 1.4;

// Which angle to converge orbits at
double alpha = 133*PI/180;	// rad

/**
 * @brief Construct an initial guess for a LTPO with a fixed alpha value from a 
 * natural orbit (e.g., a Lyapunov)
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
		crSys.getPrimary(1).c_str(), 1);
	std::vector<double> ltParams {sqrt(f), Isp},
		ctrl0 {alpha, 0};
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

Arcset_cr3bp_lt initGuessFromLTPO(const char *ltpoFile, SysData_cr3bp_lt *pSys, 
	std::vector<ControlLaw *> &laws){

	*pSys = SysData_cr3bp_lt(ltpoFile);
	Arcset_cr3bp_lt guess(pSys);
	
	guess.readFromMat(ltpoFile, laws);

	if(laws.size() != 1){
		char msg[256];
		sprintf(msg, "%zu laws were loaded from %s;\n Expected 1 law\n", 
			laws.size(), ltpoFile);

		// Free memory
		for(auto p : laws){
			delete p;
		}
		
		throw Exception(msg);
	}

	return guess;
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

	// Get the initial guess from a natural orbit
	char filename[128];
	
	// sprintf(filename, "../data/families/cr3bp_earth-moon/L1_Lyap.mat");
	// Arcset_cr3bp_lt ltGuess = initGuessFromNat(filename, 3.1, &ltSys, &law);
	// law.setType(ControlLaw_cr3bp_lt::CONST_MF_GENERAL);
	// law.setParams(std::vector<double> {sqrt(f), Isp});

	std::vector<ControlLaw *> laws {};
	sprintf(filename, 
		"../../LowThrust/MiscExplorations/RobertCollab/guess_L4_a133.00.mat");
	Arcset_cr3bp_lt ltGuess = initGuessFromLTPO(filename, &ltSys, laws);
	laws[0]->setType(ControlLaw_cr3bp_lt::CONST_MF_GENERAL);
	laws[0]->setParams(std::vector<double> {sqrt(f), Isp});
	law = *static_cast<ControlLaw_cr3bp_lt *>(laws[0]);	// Copy of the object

	Arcset_cr3bp_lt ltConverged(&ltSys);
	ltGuess.clearAllConstraints();

	// Create Multiple shooting engine
	MultShootEngine shooter;
	shooter.setVerbosity(Verbosity_tp::SOME_MSG);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);

	// Be extra careful when converging the first orbit
	shooter.setMaxIts(200);
	shooter.setDoLineSearch(true);

	printf("Converging at alpha = %.2f deg\n", alpha*180/PI);

	// Periodicity Constraint
	int id0 = ltGuess.getNodeByIx(0).getID();
	double idf = static_cast<double>(ltGuess.getNodeByIx(-1).getID());
	// std::vector<double> conData{ idf, idf, NAN, idf, NAN, NAN, NAN};
	std::vector<double> conData{ idf, idf, NAN, NAN, idf, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, id0, conData);
	ltGuess.addConstraint(periodicityCon);

	// Initial State Constraint
	std::vector<double> n0conData {NAN, NAN, 0, 0, NAN, 0, 1};
	Constraint n0con(Constraint_tp::STATE, ltGuess.getNodeByIx(0).getID(), 
		n0conData);
	ltGuess.addConstraint(n0con);

	// Remove control from design vector; fix at desired value
	std::vector<double> ctrl0 {alpha, 0};
	for(unsigned int i = 0; i < ltGuess.getNumNodes(); i++){
		Constraint con(Constraint_tp::RM_CTRL, ltGuess.getNodeByIx(i).getID(), 
			nullptr, 0);
		ltGuess.addConstraint(con);
		ltGuess.getNodeRefByIx(i).setExtraParamVec(PARAMKEY_CTRL, ctrl0);
	}

	// // Fix initial control and continuity on all segments
	// std::vector<double> ctrl0 {alpha, 0}, ctrlConData {1, 1};
	// Constraint fixCtrlCon(Constraint_tp::CTRL, 
	// 	ltGuess.getNodeByIx(0).getID(), ctrl0);
	// ltGuess.addConstraint(fixCtrlCon);
	// for(unsigned int i = 0; i < ltGuess.getNumSegs(); i++){
	// 	Constraint contCon(Constraint_tp::CONT_CTRL, 
	// 		ltGuess.getSegByIx(i).getID(), ctrlConData);
	// 	ltGuess.addConstraint(contCon);
	// }

	// Do corrections
	try{
		shooter.multShoot(&ltGuess, &ltConverged);
	}catch(DivergeException &e){
		printErr("Corrector diverged\n");
		return EXIT_SUCCESS;
	}catch(Exception &e){
		printErr("Error:\n%s\n", e.what());
		return EXIT_SUCCESS;
	}

	ltConverged.putInChronoOrder();

	// Now, pull out the constraints and adjust them so that we can converge
	// orbits over a range of H_lt but at the same alpha value
	std::vector<Constraint> allCons = ltConverged.getAllConstraints();

	// Get the low-thrust Hamiltonian value; overwrite q0
	std::vector<double> q0 = ltConverged.getSegRefByIx(0).getStateByRow(0);
	double H_lt = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q0[0]), &ltSys,
		&law);

	// Create a constraint to fix the low-thrust Hamiltonian
	Constraint HltCon(Constraint_tp::HLT, ltConverged.getNodeByIx(0).getID(),
		&H_lt, 1);
	ltConverged.addConstraint(HltCon);

	printColor(BOLDBLUE, "LT_CONVERGED\n");
	printColor(BOLDBLUE, "H_lt = %f\n", H_lt);
	printColor(BOLDBLUE, "alpha = %f deg\n", alpha*180/PI);

	ltConverged.print();
	ltConverged.saveToMat("ltConverged.mat");

	// Check Actual Periodicity
	q0 = ltConverged.getStateByIx(0);
	std::vector<double> qf = ltConverged.getStateByIx(-1);
	double diff = 0;
	assert(q0.size() == qf.size());
	for(unsigned int i = 0; i < q0.size(); i++){
		diff += pow(q0[i] - qf[i], 2);
	}
	if(sqrt(diff) > shooter.getTol()){
		printErr("Difference between q0 and qf is %e\n", sqrt(diff));
		// return EXIT_SUCCESS;
	}

	waitForUser();

	// Reset Shooter to "normal" behavior
	// shooter.setMaxIts(20);
	// shooter.setDoLineSearch(false);

	// Loop through alpha and converge orbits for all values, at same H_lt value
	double H_step = 0.001, H = 0;
	Arcset_cr3bp_lt ltDiffH(&ltSys), guess = ltConverged;
	std::vector<Arcset_cr3bp_lt> allArcs {ltConverged};
	for(H = H_lt; std::abs(H - H_lt) < H_range; H += H_step){
		HltCon.setData(&H, 1);

		guess.clearAllConstraints();
		for(unsigned int i = 0; i < allCons.size(); i++){
			guess.addConstraint(allCons[i]);
		}
		guess.addConstraint(HltCon);

		ltDiffH.reset();
		try{
			printf("Converging at H = %.4f\n", H);
			shooter.multShoot(&guess, &ltDiffH);

			if(std::abs(ltDiffH.getTotalTOF()) < 1e-2)
				throw DivergeException("TOF is too small");

			guess = ltDiffH;	// Update initiial guess
			
			// printColor(BOLDBLUE, "GUESS\n");
			// guess.print();
			// printColor(BOLDBLUE, "LT_DIFF_ALPHA\n");
			// ltDiffH.print();

			// waitForUser();
		}catch(DivergeException &e){
			printErr("Corrector diverged for H = %.4f\n", H);
			if(H_step > 0){
				H_step *= -1;
				H = H_lt;
				guess = ltConverged;
				continue;
			}else{
				break;
			}
		}catch(Exception &e){
			printErr("Error:\n%s\n", e.what());
			break;
		}

		// Insert arcs so they are in order w.r.t. alpha
		if(H_step > 0)
			allArcs.push_back(ltDiffH);
		else
			allArcs.insert(allArcs.begin(), ltDiffH);

		// Go back to zero and reverse direction if you make it to PI
		if(H_step > 0 && std::abs(H + H_step - H_lt) >= H_range){
			H_step *= -1;
			H = H_lt;
			guess = ltConverged;
		}
	}

	// Put Arcs in a family structure
	Family_PO fam(&ltSys);
	for(unsigned int i = 0; i < allArcs.size(); i++)
		fam.addMember(allArcs[i]);

	fam.sortEigs();	// Sort those eigenvalues and eigenvectors

	sprintf(filename, "L4_Lyap_f%.1e_a%06.2f_law%d.mat", f,
		alpha*180/PI, law.getType());
	fam.saveToMat(filename);

	for(auto p : laws){
		delete p;
		p = nullptr;
	}
}