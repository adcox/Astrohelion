/**
 * \file PseudoArcEngine.cpp
 * \brief
 * 
 * \author Andrew Cox
 * \version October 13, 2017
 * \copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PseudoArcEngine.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "AsciiOutput.hpp"
#include "Arcset.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "Calculations.hpp"
#include "Exceptions.hpp"
#include "Family_PO.hpp"
#include "SysData_cr3bp.hpp"
#include "MultShootData.hpp"
#include "MultShootEngine.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  \brief Construct a default pseudo-arclength engine object
 */
PseudoArcEngine::PseudoArcEngine() : ContinuationEngine(){}

/**
 *  \brief Copy a pseudo-arclength object into this one
 *  \param e reference to another pseudo-arclength engine
 */
PseudoArcEngine::PseudoArcEngine(const PseudoArcEngine &e) : ContinuationEngine(e){
	copyMe(e);
}//====================================================

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------
/**
 *  \brief Set this object equal to another
 * 
 *  \param e reference to another PAC engine
 *  \return a reference to this engine after the copy/assignment
 */
PseudoArcEngine& PseudoArcEngine::operator=(const PseudoArcEngine &e){
	copyMe(e);
	return *this;
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Tell the algorithm which nullspace vector to choose in the first
 *  iteration if the nullspace is multidimensional.
 *  \details While rare for properly formatted families, a multidimensional
 *  nullspace may occur if the initial guess is precisely a bifurcation between
 *  families. Multidimensional nullspaces may also occur if the constraints are 
 *  not independent or if the Jacobian matrix is poorly scaled.
 * 
 *  \param c which column (index begins at 0) of the nullspace to choose for
 *  the initial step
 */
void PseudoArcEngine::setNullspaceCol(unsigned int c){ nullspaceCol = c; }

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Use pseudo-arclength continuation to generate a family of symmetric periodic orbits in the CR3BP
 *  \details [long description]
 * 
 *  \param pFam pointer to a family object to store the orbits in
 *  \param pInitHalfPerGuess pointer to an arcset that represents the half-period arc of the first
 *  family member. Any constraints applied to this input guess are ignored and replaced with 
 *  mirror condition constraints
 *  \param mirrorType Describes how the orbit is symmetric
 *  \param initDir a vector that contains one non-zero element that indicates the sign of 
 *	the desired step for the family. The index corresponds to the state index. For example,
 *	if I wish the family to continue with a step in the negative z direction for the first node,
 *	I would input a vector of the form {0, 0, -1, 0, 0, 0, ...}. Technically, you can constrain
 * 	a step on any of the free variables, but only one step direction will be considered (the first one
 * 	as the algorithm reads through the vector) 
 */
void PseudoArcEngine::continueSymmetricPO_cr3bp(Family_PO *pFam, const Arcset_cr3bp *pInitHalfPerGuess, Mirror_tp mirrorType,
	std::vector<int> initDir){

	printf("Correcting Initial Guess...\n");

	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp *>(pInitHalfPerGuess->getSysData());
	Arcset_cr3bp temp_halfPer(*pInitHalfPerGuess);	// Copy the input initial guess
	temp_halfPer.clearAllConstraints();				// Clear constraints and add new ones

	// Apply mirror constraints to the half period arcset
	std::vector<unsigned int> fixedStates {};	// No fixed states
	cr3bp_addMirrorCons(&temp_halfPer, mirrorType, fixedStates);

	Arcset_cr3bp member(pSys), temp(pSys);	// create storage arcsets
	std::vector<Arcset> allArcs {};			// vector to store all half-period solutions in

	// Use pseudo-arclength continuation to evolve the half-period solution along the family
	pac(&temp_halfPer, &member, &temp, allArcs, initDir);

	// Convert the half-period solutions to full-period arcsets using symmetry properties
	Arcset_periodic perOrbit(pSys);
	for(unsigned int a = 0; a < allArcs.size(); a++){
		Arcset_cr3bp halfPerArc(allArcs[a]);
		perOrbit.reset();

		cr3bp_halfPO2fullPO(&halfPerArc, &perOrbit, mirrorType);
		pFam->addMember(perOrbit);	// Save the solution to the family or periodic orbits
	}
}//====================================================

/**
 *  \brief The core of the pseudo-arclength continuation algorithm
 *  \details This function is designed to be system- and model-agnostic
 * 
 *  \param pInitGuess pointer to an initial guess arcset. This arcset must include all
 *  the necessary constraints to converge on a family member
 *  \param pMember pointer to a storage arcset; can be empty but must be initialized and non-null
 *  \param pTemp pointer to a storage arcset; can be empty but must be initialized and non-null
 *  \param allArcs reference to a vector in which all the converged family members are stored
 *  \param initDir a vector that contains one non-zero element that indicates the sign of 
 *	the desired step for the family. The index corresponds to the state index. For example,
 *	if I wish the family to continue with a step in the negative z direction for the first node,
 *	I would input a vector of the form {0, 0, -1, 0, 0, 0, ...}. Technically, you can constrain
 * 	a step on any of the free variables, but only one step direction will be considered (the first one
 * 	as the algorithm reads through the vector) 
 * 	\param pEngine a pointer to a multiple shooting engine that will be used as a template for the 
 * 	pseudo-arclength continuation. Note that parameters required for this continuation, i.e. 
 * 	TOFType, will not be overwritten by the input engine. If left as a nullptr, no parameters
 * 	are copied.
 */
void PseudoArcEngine::pac(const Arcset *pInitGuess, Arcset *pMember, Arcset *pTemp,
	std::vector<Arcset>& allArcs, const std::vector<int> &initDir, const MultShootEngine *pEngine){

	if(pInitGuess == nullptr || pMember == nullptr || pTemp == nullptr)
		throw Exception("PseudoArcEngine::pac: input arcset pointers cannot be null");
	
	printf("Beginning Pseudo Arclength Continuation\n");
	cleanEngine();	// Reset all run-dependent parameters for the new continuation process
	double tof0 = pInitGuess->getTotalTOF();

	// Set up the multiple shooter
	MultShootEngine corrector;

	// Copy all settings from the input engine
	if(pEngine){
		corrector = *pEngine;
	}

	// Override these settings that are required for the pseudo-arclength continuation or set
	// by the PAC engine itself
	corrector.setTOFType(MSTOF_tp::VAR_EQUALARC);	// MUST use equal arc time to get propper # of constraints
	corrector.setTol(tol);
	corrector.setIgnoreCrash(true);					// Ignore crashes into primary
	corrector.setFullFinalProp(false);				// Accept minimal information from the corrector
	corrector.setVerbosity(verbosity);	// Print out convergence details

	MultShootData familyItData(pInitGuess);		// Initialize data structure to store corrections details

	// Correct the initial solution using whatever constraints apply
	// Store the converged solution in pMember, which serves as the new solution for future iterations
	// Also store the iteration data to construct the nullspace vector(s)
	try{
		familyItData = corrector.multShoot(pInitGuess, pMember);
		allArcs.push_back(*pMember);
	}catch(DivergeException &e){
		printErr("PseudoArcEngine::pac: Could not converge initial guess!\n");
		return;
	}catch(LinAlgException &e){
		printErr("PseudoArcEngine::pac: There was a linear algebra error when correcting the initial guess...\n");
		return;
	}

	printf("Applying continuation to compute family...\n");

	// Initialize storage containers
	Eigen::VectorXd convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);
	Eigen::VectorXd prevN(familyItData.totalFree, 1);

	while(orbitCount < numOrbits){

		SparseMatXCd sparseDF(familyItData.totalCons, familyItData.totalFree);
		sparseDF.setFromTriplets(familyItData.DF_elements.begin(), familyItData.DF_elements.end());
		MatrixXRd DF(sparseDF);

		/* 
		 *	The first iteration should have a DF matrix that is (n-1) x n, but all further iterations will
		 * 	have an extra row for the pseudo-arc-length constraint; we want to remove that row and take the
		 * 	nullspace of the submatrix
		 */
		if(familyItData.totalCons == familyItData.totalFree){
			DF = DF.block(0, 0, familyItData.totalCons-1, familyItData.totalFree);
		}

		// toCSV(DF, "data/DF.csv");

		// Compute nullptr space of previously computed member's Jacobian Matrix
		Eigen::FullPivLU<MatrixXRd> lu(DF);
		lu.setThreshold(1e-14);
		MatrixXRd N = lu.kernel();

		// toCSV(N, "data/N.csv");
		// waitForUser();
		
		printf("DF has dimensions %ld x %ld\n", DF.rows(), DF.cols());
		
		if(!chooseNullVec(N, initDir, prevN))
			break;

		prevN = N;	// Save the valid nullspace vector for the next iteration

		// Create new initial guess in pMember; pTemp is overwritten
		getNextPACGuess(pMember, pTemp, convergedFreeVarVec, N, familyItData);

		/*
		 *	Apply multiple shooting to converge the new guess to be a member of the family
		 */
		bool killLoop = false;
		try{
			while(stepSize >= minStepSize){
				try{
					pTemp->reset();	// Reset arcset to store converged solution

					// Attempt to correct the new member; pTemp stores the corrected result
					familyItData = corrector.multShoot(pMember, pTemp);

					// If we reach this point, the corrector CONVERGED
					// Check to see if the solution is valid (i.e., hasn't jumped to another family)
					if(!checkPACSoln(familyItData, convergedFreeVarVec, killLoop)){

						// Solution is not valid; decrease the step size if possible
						if(decreaseStepSize()){
							// Re-Create the initial guess using the new step size; pTemp is overwritten
							getNextPACGuess(pMember, pTemp, convergedFreeVarVec, N, familyItData);
						}else{
							killLoop = true;
							break;
						}
					}else{
						// ------ SOLUTION IS VALID -------

						// If convergence was fast, take bigger steps
						if(familyItData.count < static_cast<int>(stepCount_increase)){
							stepSize = stepSize*stepScaleFactor < maxStepSize ? stepSize*stepScaleFactor : maxStepSize;
							printColor(MAGENTA, "Increased Step Size to %.4e (max %.4e)!\n", stepSize, maxStepSize);
						}

						// Save new converged family vector
						convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);

						// Exit this loop; we converged!
						break;
					}
				}catch(DivergeException &e){	// The corrector did not converge
					printWarn("Corrector diverged...\n");

					if(decreaseStepSize()){
						// Re-Create the initial guess using the new step size; pTemp is overwritten
						getNextPACGuess(pMember, pTemp, convergedFreeVarVec, N, familyItData);
					}else{
						killLoop = true;
						break;
					}
				}
			}
		}catch(LinAlgException &e){
			printErr("PseudoArcEngine::pac: There was a linear algebra error...\n");
			printErr("%s\n", e.what());
			killLoop = true;
		}

		if(killLoop)
			break;

		printf("Orbit %03d converged!\n", orbitCount);

		if(pTemp->getTotalTOF()*tof0 < 0){
			printErr("Time-of-Flight changed sign; ending continuation\n");
			break;
		}

		orbitCount++;

		// Add orbit to storage vector
		allArcs.push_back(*pTemp);
	}// end of while loop
}//====================================================

/**
 *	\brief Create a arcset that contains an initial guess for a family member
 *	using pseudo arclength continuation
 *
 *	\param pNewMember a pointer to an arcset object in which to store the guess for the
 *	next family member
 *	\param pTemp a pointer to an arbitrary arcset. This arcset is cleared and overwritten; the 
 *	final contents are not important.
 *	\param convergedFreeVarVec a matrix containing the free variable vector of the previous
 *	(nearest) converged family member
 *	\param N a 1D nullspace vector that lies tangent to the family
 *	\param pFamilyItData pointer to a MultShootData object containing corrections information about the
 *	previous (nearest) converged family member. This data, combined with the converged freeVarVec, is leveraged
 *	to construct the new solution
 *	\param it reference to the iteration data that describes the converged solution
 */
void PseudoArcEngine::getNextPACGuess(Arcset *pNewMember, Arcset *pTemp, const Eigen::VectorXd &convergedFreeVarVec,
	const Eigen::VectorXd &N, MultShootData& it){

	/**
	 *	Step away from previously converged solution
	 */

	Eigen::VectorXd newFreeVarVec = convergedFreeVarVec + stepSize*N;
	double *X = newFreeVarVec.data();
	it.X = std::vector<double>(X, X + newFreeVarVec.rows());

	pTemp->reset();	// clean everything out
	it.nodesOut = pTemp;

	pNewMember->getSysData()->getDynamicsModel()->multShoot_createOutput(it);

	// Get rid of any pre-existing pseudo arclength constraints from previous corrections
	std::vector<Constraint> arcCons = pTemp->getArcConstraints();
	pTemp->clearArcConstraints();
	for(Constraint &con : arcCons){
		if(con.getType() != Constraint_tp::PSEUDOARC)
			pTemp->addConstraint(con);
	}

	/* 
	 *	Form the Pseudo-Arclength Continuation constraint
	 */
	std::vector<double> pacCon_data = it.X;
	// Append the nullptr vector (i.e. step direction)
	pacCon_data.insert(pacCon_data.end(), N.data(), N.data()+N.rows());
	// Append the step size
	pacCon_data.insert(pacCon_data.end(), stepSize);
	// Create the actual constraint
	Constraint pacCon(Constraint_tp::PSEUDOARC, it.numNodes-1, pacCon_data);
	pTemp->addConstraint(pacCon);

	// Outputs for debugging and sanity checks
	printf("New IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...] tof = %.4f\n",
		X[0], X[1], X[2], X[3], X[4], X[5], X[newFreeVarVec.rows()-1]);

	*pNewMember = *pTemp;
}//====================================================

/**
 *  \brief Decrease the step size, if possible
 *  
 *  \return whether or not the step size was decreased. A return value
 *  of FALSE indicates that the step size cannot be decreased any further
 */
bool PseudoArcEngine::decreaseStepSize(){
	if(stepSize > minStepSize){
		printWarn("trying smaller step size\n");

		// Decrease step size and try again
		stepSize = stepSize/stepScaleFactor > minStepSize ? stepSize/stepScaleFactor : minStepSize;
		printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

		return true;	// success!
	}else{
		printColor(MAGENTA, "Reached minimum step size!\n");
		return false;	// Can't go any smaller...
	}
}//====================================================

/**
 *  \brief Check to see if the converged solution is valid
 * 
 *  \param it Iteration data from the most recent multiple shooting process
 *  \param convergedFreeVarVec Free variable vector describing the previous 
 *  converged and accepted solution
 *  \param killLoop Whether or not to end the continuation. For example,
 *  if the converged solutions have ceased changing then killLoop is set 
 *  to true.
 *  \return whether or not the converged solution is valid
 */
bool PseudoArcEngine::checkPACSoln(const MultShootData &it, const Eigen::VectorXd &convergedFreeVarVec,
	bool& killLoop){

	killLoop = false;

	// Check to see if the converged family vector is significantly different from previously computed family member
	// Note that only the core states (the IC for the trajectory) is checked; differences in other nodes 
	// are assumed to be insignificant (if IC is the same, only possible change is change in TOF)
	double sumSquared = 0;
	unsigned int core_dim = it.nodesIn->getSysData()->getDynamicsModel()->getCoreStateSize();
	for(unsigned int i = 0; i < core_dim; i++){
		sumSquared += (it.X[i] - convergedFreeVarVec[i])*(it.X[i] - convergedFreeVarVec[i]);
	}
	double diffX = sqrt(sumSquared);
	printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDYELLOW, "||diff in X(1:%u)|| = %.4e\n", core_dim, diffX);

	if(diffX < tol){
		printErr("Solutions from pseudo-arc-length have ceased changing; ending continuation\n");
		killLoop = true;	// End continuation
		return true;		// The solution is fine
	}

	if(diffX > 10*stepSize){
		printErr("Solution changed by %.4e > 10*(stepSize = %.4e).\n", diffX, stepSize);
		return false;	// Solution is no good
	}

	return true;
}//====================================================

/**
 *  \brief Check to see if the converged solution is periodic
 * 
 *  \param pSol Converged solution from the most recent multiple shooting process
 *  \param killLoop Whether or not to end the continuation. For example,
 *  if the converged solutions have ceased changing then killLoop is set 
 *  to true.
 *  \return whether or not the converged solution is periodic
 */
bool PseudoArcEngine::checkPACSoln_periodic(const Arcset_periodic *pSol, bool& killLoop){

	killLoop = false;
	
	// Check periodicity
	std::vector<double> state0 = pSol->getStateByIx(0);
	std::vector<double> statef = pSol->getStateByIx(-1);
	double sumSquared = 0;
	for(unsigned int i = 0; i < state0.size(); i++){
		sumSquared += (state0[i] - statef[i])*(state0[i] - statef[i]);
	}

	if(sqrt(sumSquared) > 10*tol){
		printColor(RED, "NOT PERIODIC! Error = %.4e\n", sqrt(sumSquared));
		return false;
	}

	return true;	// Default: Solution is ok!
}//====================================================

/**
 *  \brief Choose the nullspace vector as well as the correct sign
 *  \details Use the orientation of the nullspace vector in n-space to 
 *  determine the correct sign (i.e., choose the sign that results in the
 *  closest orientation to the previous converged member). Multi-dimensional
 *  nullspaces are not currently handled, though could be via the same logic.
 * 
 *  \param N A matrix that stores each nullspace vector as a column
 *  \return whether or not the selection process was successful. A return value of FALSe
 *  indicates that the continuation should end
 */
bool PseudoArcEngine::chooseNullVec(MatrixXRd &N, std::vector<int> initDir, const MatrixXRd &prevN){
	// Check to make sure the IS a nullspace
	if(N.rows() == 1){
		printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: nullspace is zero-dimensional; cannot proceed...\n");
		return false;
	}		

	// // For debugging, save nullspace vectors to file
	// char filename[32];
	// sprintf(filename, "N%02d.csv", orbitCount);
	// toCSV(N, filename);

	/**
	 *	Choose the nullspace vector that is closest to the previous one (which converged)
	 */
	printf("Choosing nullspace Vector (%ldD, %ld elements)\n", N.cols(), N.rows());
	if(orbitCount == 0){
		if(N.cols() > 1){
			printErr("Nullspace is multidimensional on first iteration; choosing column %u\n", nullspaceCol);
			if(nullspaceCol < N.cols()){
				Eigen::VectorXd tempVec = N.col(nullspaceCol);
				N = tempVec;
			}else{
				printErr("Nullspace column = %u is out of bounds (N has %d cols)\n", nullspaceCol, N.cols());
				return false;
			}
		}

		bool sameDir = true;
		for(unsigned int i = 0; i < initDir.size(); i++){
			// Find a non-zero element
			if(initDir[i] != 0 && i < static_cast<unsigned int>(N.rows())){
				// If signs are different, assume direction is different
				if(N(i,0)*initDir[i] < 0){
					sameDir = false;
					break;
				}
			}
		}

		// Reverse direction of nullspace
		if(!sameDir)
			N *= -1;

	}else{
		/* Make sure nullspace direction stays consistent by choosing the 
		 * nullptr vector that is closest to the same direction as the previous one
		 */
		int best_ix = 0;	// index of the column of the best vector option
		int best_sign = 1;	// sign associated with best vector
		double best_angle = 181;	// the best (smallest) angle found
		for(int i = 0; i < N.cols(); i++){
			// Compute angle from dot product
			Eigen::VectorXd col_i = N.cols() > 1 ? N.col(i) : N;
			Eigen::VectorXd dotProd = prevN.transpose()*col_i / (prevN.norm()*N.norm());
			double angle = std::acos(dotProd(0));
			int sign = 1;

			// Flip the sign if the angle is greater than 90 and change the value to 180 - angle
			printf("dot product = %.4f\n", dotProd(0));
			printf("Angle = %.4f deg\n", angle*180/PI);
			if(angle > PI/2.0){
				printColor(CYAN, "Angle is %.4f > pi/2; changing sign and angle\n", angle);
				angle = PI - angle;
				sign = -1;
			}

			// Keep track of the best option
			if(angle < best_angle){
				best_angle = angle;
				best_sign = sign;
				best_ix = i;
			}
		}

		printf("best ix = %d, sign = %d\n", best_ix, best_sign);
		Eigen::VectorXd tempVec = N.cols() > 1 ? N.col(best_ix) : N;
		N = best_sign*tempVec;	// Apply sign change, if needed
	}

	// Turn N into a unit vector
	N = N/N.norm();

	printf("Chose N with first elements = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...]\n",
			N(0), N(1), N(2), N(3), N(4), N(5));

	return true;
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Make a copy of the engine
 *  \param e reference to another engine object to copy into this one
 */
void PseudoArcEngine::copyMe(const PseudoArcEngine &e){
	ContinuationEngine::copyMe(e);
	stepSize = e.stepSize;
	orbitCount = e.orbitCount;
}//====================================================

/**
 *  \brief Reset all parameters, including user-specified options and run-specific iterators
 */
void PseudoArcEngine::reset(){
	ContinuationEngine::reset();
	stepSize = 0.001;
	orbitCount = 0;
}//====================================================

/**
 *  \brief Reset only run-specific iterators
 *  \details User-specified options (e.g., step size limits) are not changed
 */
void PseudoArcEngine::cleanEngine(){
	ContinuationEngine::cleanEngine();

	orbitCount = 0;
}//====================================================

}// End of astrohelion namespace