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

PseudoArcEngine::PseudoArcEngine() : ContinuationEngine(){}

PseudoArcEngine::PseudoArcEngine(const PseudoArcEngine &e) : ContinuationEngine(e){
	copyMe(e);
}//====================================================

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

PseudoArcEngine& PseudoArcEngine::operator=(const PseudoArcEngine &e){
	copyMe(e);
	return *this;
}//====================================================

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

void PseudoArcEngine::continueSymmetricPO_cr3bp(Family_PO *pFam, const Arcset_cr3bp *pInitHalfPerGuess, Mirror_tp mirrorType, std::vector<int> initDir){
	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp *>(pInitHalfPerGuess->getSysData());

	cleanEngine();

	Arcset_cr3bp temp_halfPer(*pInitHalfPerGuess);	// Copy the input initial guess
	double tof0 = pInitHalfPerGuess->getTotalTOF();

	printf("Correcting Initial Guess...\n");
	
	// Clear constraints and add new ones
	temp_halfPer.clearAllConstraints();

	/*	Constraint Method 1:
	 *	* Apply a periodicity constraint that forces the first and final node to be collocated,
	 *	  	ignoring one to avoid numerical troubles. It is best to ignore one of the planar 
	 *		position components, especially for planar families - ignoring z or z-dot can shift the
	 *		nullptr vector to have non-zero elements in the z and z-dot spots, effectively stepping 
	 *		out of the plane, which is not desireable for a planar family...
	 *	* Constrain one extra state to be fixed; I used something that makes sense for a perpendicular
	 *		plane crossing here (avoid constraining z or z-dot to be zero for planar families)
	 *
	 *	These were the descriptions for three input ints for this method:
	 *
	 *	\param periodicityIgnoreIx (int) index of a state variable to ignore when enforcing periodicity. It is best
	 *	to ignore one of the planar components (i.e. x or y, corresponding to indices 0 and 1, respectively)
	 *	\param fixToVal_ix (int) index of a variable to fix to <tt>fixToVal_val</tt> at the first node. If the 
	 *	index is between 0 and 5, it will constraint one of the usual state variables. An index of 6 will 
	 *	constrain total TOF, and an index of 7 will constrain Jacobi at the first node
	 *	\param fixToVal_val (double) the value to constrain state <tt>fixToVal_ix</tt> to
	 */
	 /*
	// Create a periodicity constraint
	double periodicConData[] = {0,0,0,0,0,0};
	periodicConData[periodicityIgnoreIx] = NAN;
	Constraint periodicCon(Constraint_tp::MATCH_CUST, temp_halfPer.getNumNodes()-1, periodicConData, 6);

	Constraint extraCon;
	if(fixToVal_ix < 6){
		double extraCon_data[] = {NAN, NAN, NAN, NAN, NAN, NAN};
		extraCon_data[fixToVal_ix] = fixToVal_val;
		extraCon = Constraint(Constraint_tp::STATE, 0, extraCon_data, 6);
	}else if(fixToVal_ix == 6){
		double val = fixToVal_val;
		extraCon = Constraint(Constraint_tp::TOF_TOTAL, 0, &val, 1);
	}else if(fixToVal_ix == 7){
		double val = fixToVal_val;
		extraCon = Constraint(Constraint_tp::JC, 0, &val, 1);
	}

	temp_halfPer.addConstraint(periodicCon);
	temp_halfPer.addConstraint(extraCon);
	*/

	/* Constraint Method 2:
	 *
	 *	Apply two constraints that enforce perpendicular crossings at the initial and final nodes
	 *	This gives MUCH better performance for the Lyapunov families
	 */
	std::vector<unsigned int> fixedStates {};	// No fixed states
	cr3bp_addMirrorCons(&temp_halfPer, mirrorType, fixedStates);

	// Correct the arcset to retrieve a free-variable vector for a family member
	MultShootEngine corrector;
	corrector.setTOFType(MSTOF_tp::VAR_EQUALARC);	// MUST use equal arc time to get propper # of constraints
	corrector.setTol(tol);
	corrector.setIgnoreCrash(true);		// Ignore crashes into primary
	corrector.setFullFinalProp(false);	// Accept minimal information from the corrector
	corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	MultShootData familyItData(&temp_halfPer);
	Arcset_cr3bp halfPer_corrected(pSys);

	// Initialize this arcset outside the loop because the familyItData will end up with a pointer
	// to this arcset after the multiple shooting processs; if the declaration is in the loop,
	// the arcset is destroyed each iteration and the pointer ceases to be useful.
	Arcset_periodic perOrbit(pSys);
	Arcset_cr3bp newMember(pSys);

	try{
		familyItData = corrector.multShoot(&temp_halfPer, &halfPer_corrected);
	}catch(DivergeException &e){
		printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: Could not converge initial guess!\n");
	}catch(LinAlgException &e){
		printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: There was a linear algebra error...\n");
	}

	printf("Applying continuation to compute family...\n");
	
	// Initialize storage containers
	Eigen::VectorXd convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);
	Eigen::VectorXd prevN(familyItData.totalFree, 1);
	
	std::vector<Arcset_cr3bp> members;

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

		// Compute nullptr space of previously computed member's Jacobian Matrix
		Eigen::FullPivLU<MatrixXRd> lu(DF);
		lu.setThreshold(1e-14);
		MatrixXRd N = lu.kernel();

		printf("DF has dimensions %ld x %ld\n", DF.rows(), DF.cols());
		
		if(!chooseNullVec(N, initDir, prevN))
			break;

		prevN = N;	// Update memory

		newMember = getNextPACGuess_cr3bp(convergedFreeVarVec, N, stepSize, &familyItData);

		/*
		 *	Apply multiple shooting to converge the new guess to be a member of the family
		 */
		bool killLoop = false;
		try{
			while(stepSize >= minStepSize){
				try{

					// Reset perOrbit
					perOrbit = Arcset_periodic(pSys);
					Arcset_cr3bp temp_converged(pSys);

					familyItData = corrector.multShoot(&newMember, &temp_converged);
					cr3bp_halfPO2fullPO(&temp_converged, &perOrbit, mirrorType);

					if(!checkPACSolution_cr3bp(&perOrbit, &familyItData, convergedFreeVarVec, stepSize, &killLoop)){

						if(stepSize > minStepSize){
							stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
							printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

							// Re-Create the initial guess using the new step size
							newMember = getNextPACGuess_cr3bp(convergedFreeVarVec, N, stepSize, &familyItData);
							continue;	// Go to next iteration
						}else{
							printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: Could not converge new family member!\n");
							killLoop = true;
							break;
						}
					}else{

						// If convergence was fast, take bigger steps
						if(familyItData.count < 5){
							stepSize = stepSize*2 < maxStepSize ? stepSize*2 : maxStepSize;
							printColor(MAGENTA, "Increased Step Size to %.4e (max %.4e)!\n", stepSize, maxStepSize);
						}else if(familyItData.count > 10){
							// If convergence was slow, take smaller steps
							stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
							printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);
						}

						// Exit this loop; we converged!
						break;
					}
				}catch(DivergeException &e){
					if(stepSize > minStepSize){
						astrohelion::printWarn("PseudoArcEngine::continueSymmetricPO_cr3bp: Corrector diverged... trying smaller step size\n");

						// Decrease step size and try again
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

						// Re-Create the initial guess using the new step size
						newMember = getNextPACGuess_cr3bp(convergedFreeVarVec, N, stepSize, &familyItData);
					}else{
						printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: Could not converge new family member!\n");
						killLoop = true;
						break;
					}
				}
			}
		}catch(LinAlgException &e){
			printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: There was a linear algebra error...\n");
			killLoop = true;
		}

		if(killLoop)
			break;

		printf("Orbit %03d converged!\n", orbitCount);

		// Save new converged family vector
		convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);

		if(perOrbit.getTotalTOF()*tof0 < 0){
			printErr("Time-of-Flight changed sign; ending continuation\n");
			break;
		}

		orbitCount++;

		// Compute eigenvalues
		MatrixXRd mono = perOrbit.getMonodromy();
		
		double monoErr = std::abs(1.0 - mono.determinant());
		if(monoErr > 1e-5)
			printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);

		// Add orbit to family
		pFam->addMember(perOrbit);
	}// end of while loop
}//====================================================

/**
 *	\brief Create a arcset that contains an initial guess for a family member
 *	using pseudo arclength continuation
 *
 *	\param convergedFreeVarVec a matrix containing the free variable vector of the previous
 *	(nearest) converged family member
 *	\param N a 1D nullspace vector that lies tangent to the family
 *	\param stepSize scales the size of the step by scaling the nullspace vector
 *	\param pFamilyItData pointer to a MultShootData object containing corrections information about the
 *	previous (nearest) converged family member
 */
Arcset_cr3bp PseudoArcEngine::getNextPACGuess_cr3bp(const Eigen::VectorXd &convergedFreeVarVec,
	const Eigen::VectorXd &N, double stepSize, MultShootData *pFamilyItData){

	/**
	 *	Step away from previously converged solution
	 */

	Eigen::VectorXd newFreeVarVec = convergedFreeVarVec + stepSize*N;
	double *X = newFreeVarVec.data();
	pFamilyItData->X = std::vector<double>(X, X + newFreeVarVec.rows());

	// Convert into a new arcset (TODO: Make this more flexible by putting conversion code in a model?)
	const SysData_cr3bp *sys = static_cast<const SysData_cr3bp *>(pFamilyItData->nodesIn->getSysData());
	Arcset_cr3bp newMember(sys);
	pFamilyItData->nodesOut = &newMember;

	sys->getDynamicsModel()->multShoot_createOutput(pFamilyItData);

	// Get rid of any pre-existing pseudo arclength constraints from previous corrections
	std::vector<Constraint> arcCons = newMember.getArcConstraints();
	newMember.clearArcConstraints();
	for(Constraint &con : arcCons){
		if(con.getType() != Constraint_tp::PSEUDOARC)
			newMember.addConstraint(con);
	}

	/* 
	 *	Form the Pseudo-Arclength Continuation constraint
	 */
	std::vector<double> pacCon_data = pFamilyItData->X;
	// Append the nullptr vector (i.e. step direction)
	pacCon_data.insert(pacCon_data.end(), N.data(), N.data()+N.rows());
	// Append the step size
	pacCon_data.insert(pacCon_data.end(), stepSize);
	// Create the actual constraint
	Constraint pacCon(Constraint_tp::PSEUDOARC, pFamilyItData->numNodes-1, pacCon_data);
	newMember.addConstraint(pacCon);

	// Outputs for debugging and sanity checks
	printf("New IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...] tof = %.4f\n",
		X[0], X[1], X[2], X[3], X[4], X[5], X[newFreeVarVec.rows()-1]);

	return newMember;
}//====================================================

/**
 *  \brief Check to see if the converged solution is valid
 * 
 *  \param pSol Converged solution from the most recent multiple shooting process
 *  \param pIt Iteration data from the most recent multiple shooting process
 *  \param convergedFreeVarVec Free variable vector describing the previous 
 *  converged and accepted solution
 *  \param killLoop Whether or not to end the continuation. For example,
 *  if the converged solutions have ceased changing then killLoop is set 
 *  to true.
 *  \return whether or not the converged solution is valid
 */
bool PseudoArcEngine::checkPACSolution_cr3bp(const Arcset_periodic *pSol, const MultShootData *pIt,
	const Eigen::VectorXd &convergedFreeVarVec, double stepSize, bool *killLoop){

	*killLoop = false;
	
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

	if(pIt->totalFree < 6){
		char msg[128];
		sprintf(msg, "PseudoArcEngine::checkPACSolution_cr3bp: there are %d free variables; expect 6!", pIt->totalFree);
		throw Exception(msg);
	}

	// Check to see if the converged family vector is significantly different from previously computed family member
	// Note that only the first 6 states (the IC for the trajectory) is checked; differences in other nodes 
	// are assumed to be insignificant (if IC is the same, only possible change is change in TOF)
	sumSquared = 0;
	for(unsigned int i = 0; i < 6; i++){
		sumSquared += (pIt->X[i] - convergedFreeVarVec[i])*(pIt->X[i] - convergedFreeVarVec[i]);
	}
	double diffX = sqrt(sumSquared);
	astrohelion::printErr("||diff in X(1:6)|| = %.4e\n", diffX);

	if(diffX < tol){
		astrohelion::printErr("Solutions from pseudo-arc-length have ceased changing; ending continuation\n");
		*killLoop = true;	// End continuation
		return true;		// The solution is fine
	}

	if(diffX > 2*stepSize){
		astrohelion::printErr("Solution changed by %.4e > 2*(stepSize = %.4e).\n", diffX, stepSize);
		return false;	// Solution is no good
	}

	// Default: Solution is ok!
	return true;
}//====================================================

/**
 *  \brief Choose the nullspace vector as well as the correct sign
 *  \details [long description]
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
	// char filename[16];
	// sprintf(filename, "N%02d.csv", orbitCount);
	// N.astrohelion::toCSV(filename);

	/**
	 *	Choose the nullspace vector that is closest to the previous one (which converged)
	 */
	printf("Choosing nullspace Vector (%ldD, %ld elements)\n", N.cols(), N.rows());
	if(orbitCount == 0){
		if(N.cols() > 1){
			printErr("PseudoArcEngine::continueSymmetricPO_cr3bp: nullspace is multidimensional on first iteration; unsure how to proceed...\n");
			return false;
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
		Eigen::VectorXd temp = N.cols() > 1 ? N.col(best_ix) : N;
		N = best_sign*temp;	// Apply sign change, if needed
	}

	printf("Chose N with first elements = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...]\n",
			N(0), N(1), N(2), N(3), N(4), N(5));

	return true;
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

void PseudoArcEngine::copyMe(const PseudoArcEngine &e){
	ContinuationEngine::copyMe(e);
	stepSize = e.stepSize;
	orbitCount = e.orbitCount;
}//====================================================

void PseudoArcEngine::reset(){
	ContinuationEngine::reset();
	stepSize = 0.001;
}//====================================================

void PseudoArcEngine::cleanEngine(){
	ContinuationEngine::cleanEngine();

	orbitCount = 0;
}//====================================================

}// End of astrohelion namespace