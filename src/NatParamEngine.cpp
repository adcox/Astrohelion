/**
 * \file NatParamEngine.cpp
 * \brief
 * 
 * \author Andrew Cox
 * \version June 27, 2017
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

#include "NatParamEngine.hpp"

#include <Eigen/Dense>

#include "AsciiOutput.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "Calculations.hpp"
#include "Exceptions.hpp"
#include "Family_PO.hpp"
#include "MultShootData.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  \brief Default constructor
 */
NatParamEngine::NatParamEngine() : ContinuationEngine() {}

NatParamEngine::NatParamEngine(const NatParamEngine &e) : ContinuationEngine(e) {
	copyMe(e);
}//====================================================

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

NatParamEngine& NatParamEngine::operator=(const NatParamEngine &e){
	copyMe(e);
	return *this;
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Set the number of solutions to store in the "curve fit memory"
 *  \details The number of solutions saved informs the least-squares
 *  curve-fitting algorithm that improves the performance of the natural 
 *  parameter continuation.
 * 
 *  \param mem Number of solutions to use when fitting a curve to predict
 *  the continuation variables
 */
void NatParamEngine::setCurveFitMem(unsigned int mem){ curveFitMem = mem; }

/**
 *  \brief Set the number of simple continuation (i.e., no curve fitting)
 *  iterations to perform before switching to curve-fitted solutions
 * 
 *  \param num the number of simple continuation iterations to perform
 */
void NatParamEngine::setNumSimple(unsigned int num){ numSimple = num; }

/**
 *  \brief Set the slope threshold that is the limit between
 *  two step directions
 *  \details The natural parameter continuation typically leverages 
 *  several "independent" stepping variables to facilitate continuation 
 *  along a nonlinear contour. When a step in one variable results in
 *  large changes in the other, the stepping variable is changed to 
 *  limit the sensitivity.
 * 
 *  \param thresh slope threshold between stepping strategies
 */
void NatParamEngine::setSlopeThresh(double thresh){ slopeThresh = thresh; }

/**
 *  \brief Set the step size for simple (no curve-fitting) continuation
 *  \param step the step size for simple (no curve-fitting) continuation
 */
void NatParamEngine::setStep_simple(double step){ step_simple = step; }

/**
 *  \brief Set the step size for fitted continuation: variable 1
 *  \param step the step size for fitted continuation: variable 1
 */
void NatParamEngine::setStep_fitted_1(double step){ step_fitted_1 = step; }

/**
 *  \brief Set the step size for fitted continuation: variable 2
 *  \param step the step size for fitted continuation: variable 2
 */
void NatParamEngine::setStep_fitted_2(double step){ step_fitted_2 = step; }

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *	\brief Continue a family of periodic orbits via natural parameter continuation
 *	
 * 	\param fam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	\param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	\param mirrorTypes a vector of variables that describe how the family mirrors in the rotating 
 *	reference frame. Each entry corresponds to an independent variable in <tt>indVarIx</tt>
 *	\param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	\param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 *
 *	\throws Exception if <tt>indVarIx</tt> has fewer than two elements
 *	\throws Exception if <tt>mirrorTypes</tt> does not have the same size as <tt>indVarIx</tt>
 *	\throws Exception if the eigenvalues of the monodromy matrix cannot be computed
 *	\throws Exception if one of the indices stored in <tt>indVarIx</tt> or <tt>depVarIx</tt> is
 *	out of range
 */
void NatParamEngine::continueSymmetricPO_cr3bp(Family_PO *fam, const Arcset_cr3bp *initialGuess,
	std::vector<Mirror_tp> mirrorTypes, std::vector<unsigned int> indVarIx, std::vector<unsigned int> depVarIx){

	// Assume family is CR3BP
	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp*>(fam->getSysData());

	if(indVarIx.size() < 2)
		throw Exception("FamGenerator::cr3bp_natParamCont: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw Exception("FamGenerator::cr3bp_natParamCont: there must be an equal number of ind. vars and mirror types");

	int indVar1 = indVarIx[0];
	int indVar2 = indVarIx[1];
	Mirror_tp mirrorType = mirrorTypes[0];

	// Initially assume that we're fixing indVar1
	std::vector<unsigned int> fixStates;
	fixStates.push_back(indVar1);

	// Get info from the initial guess trajectory
	std::vector<double> IC = initialGuess->getStateByIx(0);
	double tof = initialGuess->getTotalTOF();
	double tof0 = tof;

	// Initialize counters and storage containers
	unsigned int orbitCount = 0;
	double indVarSlope = NAN, deltaVar1 = 1, deltaVar2 = 1;

	std::vector<Arcset_periodic> members;
	bool diverged = false;

	// Create a dummy arcset and create an iteration data object on the stack
	// The cr3bp_getPeriodic() function will only pass an iteration data pointer back
	// if the one passed in is not nullptr, hence we create a valid object and delete it
	// before exiting the function
	Arcset_cr3bp halfPerGuess(*initialGuess);
	Arcset_cr3bp tempCorrected(static_cast<const SysData_cr3bp *>(initialGuess->getSysData()));
	Arcset_cr3bp halfPerCorrected(static_cast<const SysData_cr3bp *>(initialGuess->getSysData()));
	MultShootData *pItData = new MultShootData(&halfPerGuess);

	while(orbitCount < numOrbits){
		Arcset_periodic perOrbit(pSys);
		try{
			printf("Guess for IC: [%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f] %.4f\n", IC[0], IC[1], IC[2], IC[3],
				IC[4], IC[5], tof);
			printf("Fix States: ");
			for(unsigned int i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			// Add constraints to the guess to enforce the mirror conition and fix the specified states
			cr3bp_addMirrorCons(&halfPerGuess, mirrorType, fixStates);
			// Correct using multiple shooting, save the corrected half-period arc for next time
			perOrbit = cr3bp_getSymPO(&halfPerGuess, &tempCorrected, mirrorType, tol, pItData);

			diverged = false;
			halfPerCorrected = tempCorrected;
			printf("Orbit %03d converged!\n", static_cast<int>(members.size()));
		}catch(DivergeException &e){
			diverged = true;
		}catch(LinAlgException &e){
			astrohelion::printErr("There was a linear algebra error during family continuation...\n");
			break;
		}

		// Check for large changes in period to detect leaving family
		if(!diverged && orbitCount > 2){
			// difference in TOF; use abs() because corrector may employ reverse time and switch to forward time
			double dTOF = std::abs(perOrbit.getTotalTOF()) - std::abs(members[members.size()-1].getTotalTOF());
			double percChange = std::abs(dTOF/perOrbit.getTotalTOF());
			if(percChange > 0.25){
				printWarn("percChange = %.4f\n", percChange);
				printWarn("Period jumped (now = %.5f)! Left the family! Trying smaller step size...\n", perOrbit.getTotalTOF());
				diverged = true;
			}
		}

		if(diverged && orbitCount == 0){
			astrohelion::printErr("Could not converge on the first family member; try a smaller step size\n");
			break;
		}else if(diverged && orbitCount > 0){
			perOrbit = members.back();	// Use previous solution as the converged solution to get correct next guess

			if(orbitCount <= numSimple){
				if(step_simple > minStepSize){
					step_simple = step_simple/2 > minStepSize ? step_simple/2 : minStepSize;
					printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", step_simple, minStepSize);
				}else{
					printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}else{
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq > minStepSize){
					dq = dq/2 > minStepSize ? dq/2 : minStepSize;
					printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", dq, minStepSize);

					if(std::abs(indVarSlope) > slopeThresh)
						step_fitted_1 = dq;
					else
						step_fitted_2 = dq;
				}else{
					printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}
		}else{	// Did not diverge

			// Save the computed orbit
			members.push_back(perOrbit);
			orbitCount++;

			// Check to see if we should update the step size
			if(pItData->count < 4 && orbitCount > numSimple){
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq < maxStepSize){
					dq = 2*dq < maxStepSize ? 2*dq : maxStepSize;
					printColor(MAGENTA, "  Increased step size to %0.4e (max = %.4e)\n", dq, maxStepSize);
					if(std::abs(indVarSlope) > slopeThresh)
						step_fitted_1 = dq;
					else
						step_fitted_2 = dq;
				}
			}

			// Compute eigenvalues
			MatrixXRd mono = perOrbit.getMonodromy();

			double monoErr = std::abs(1.0 - mono.determinant());
			if(monoErr > 1e-5)
				printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);

			// Add orbit to family
			fam->addMember(perOrbit);
		}

		// Create next initial guess
		tof = perOrbit.getTimeByIx(-1);

		if(tof*tof0 < 0){
			printErr("Time-of-Flight changed sign: ending continuation process\n");
			break;
		}

		if(orbitCount < numSimple){
			// Use simple continuation; copy the converged IC, step forward in the independent variable
			IC = perOrbit.getStateByIx(0);
			IC.at(indVar1) += step_simple;
		}else{

			// Compute the slope for the first time
			if(orbitCount == numSimple){
				deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVar1] - members[orbitCount-2].getStateByIx(0)[indVar1];
				deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVar2] - members[orbitCount-2].getStateByIx(0)[indVar2];
				indVarSlope = deltaVar1/deltaVar2;
			}

			// Use least-squares fitting to predict the value for the independent and dependent variables
			std::vector<double> prevStates;
			int first = static_cast<int>(members.size()) - static_cast<int>(curveFitMem) < 0 ? 0 : members.size() - curveFitMem;

			for(unsigned int n = first; n < members.size(); n++){
				std::vector<double> ic = members[n].getStateByIx(0);
				prevStates.insert(prevStates.end(), ic.begin(), ic.begin()+6);
				prevStates.push_back(members[n].getTimeByIx(-1));

				Arcset_cr3bp temp = static_cast<Arcset_cr3bp>(members[n]);
				prevStates.push_back(temp.getJacobiByIx(0));
			}

			// This will hold the input depVars plus the unused independent variable
			std::vector<unsigned int> allDepVars;
			std::vector<double> predictedIC;
			if(std::abs(indVarSlope) > slopeThresh){
				mirrorType = mirrorTypes[0];
				// Use continuation in indVar1
				IC.at(indVar1) = perOrbit.getStateByIx(0).at(indVar1) + astrohelion::sign(deltaVar1)*step_fitted_1;
				fixStates.clear();
				fixStates.push_back(indVar1);

				// Use Least Squares to predict dependent vars and unusued ind. vars
				allDepVars = depVarIx;
				if(std::find(depVarIx.begin(), depVarIx.end(), indVar2) == depVarIx.end()){
					allDepVars.push_back(indVar2);	// only add if it isn't already part of depVarIx
				}
				predictedIC = familyCont_LS(indVar1, IC.at(indVar1), allDepVars, prevStates);
			}else{
				mirrorType = mirrorTypes[1];
				// Use continuation in indVar2
				IC.at(indVar2) = perOrbit.getStateByIx(0).at(indVar2) + astrohelion::sign(deltaVar2)*step_fitted_2;
				fixStates.clear();
				fixStates.push_back(indVar2);

				// Use Least Squares to predict dependent vars and unusued ind. vars
				allDepVars = depVarIx;
				if(std::find(depVarIx.begin(), depVarIx.end(), indVar1) == depVarIx.end()){
					allDepVars.push_back(indVar1);	// only add if it isn't already part of depVarIx
				}
				predictedIC = familyCont_LS(indVar2, IC.at(indVar2), allDepVars, prevStates);
			}

			// Update IC with predicted variables
			for(unsigned int n = 0; n < allDepVars.size(); n++){
				unsigned int ix = allDepVars[n];
				if(ix < 6)
					IC[ix] = predictedIC[ix];
				else if(ix == 6)
					tof = predictedIC[ix];
				else{
					if(pItData){
						delete(pItData);
						pItData = nullptr;
					}
					throw Exception("FamGenerator::cr3bp_natParamCont: Cannot update independent variable; index out of range");
				}
			}

			// Update slope
			deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVar1] - members[orbitCount-2].getStateByIx(0)[indVar1];
			deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVar2] - members[orbitCount-2].getStateByIx(0)[indVar2];
			indVarSlope = deltaVar1/deltaVar2;
		}

		// Update the initial guess
		halfPerGuess = halfPerCorrected;
		halfPerGuess.setStateByIx(0, IC);	// Use the new guess for the initial state

		// Change the TOF on the final segment
		// WARNING - might cause issues of TOF goes negative
		double dt = tof/2.0 - halfPerGuess.getTotalTOF();
		halfPerGuess.getSegRefByIx(-1).setTOF(halfPerGuess.getTOFByIx(-1) + dt);
		
		if(halfPerGuess.getTOFByIx(-1)*halfPerCorrected.getTOFByIx(-1) < 0)
			printWarn("NatParamEngine::continueSymmetricPO_cr3bp: updated TOF has changed sign\n");
	}// end of while loop

	if(pItData){
		delete(pItData);
		pItData = nullptr;
	}
}//==================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Make a copy of the natural parameter engine
 * 
 *  \param engine the source of the copy
 */
void NatParamEngine::copyMe(const NatParamEngine &engine){
	ContinuationEngine::copyMe(engine);
	curveFitMem = engine.curveFitMem;
	numSimple = engine.curveFitMem;
	slopeThresh = engine.slopeThresh;
	step_simple = engine.step_simple;
	step_fitted_1 = engine.step_fitted_1;
	step_fitted_2 = engine.step_fitted_2;
}//====================================================

/**
 *  \brief Reset all continuation parameters
 */
void NatParamEngine::reset(){
	ContinuationEngine::reset();
	curveFitMem = 5;
	numSimple = 3;
	slopeThresh = 1;
	step_simple = 5e-4;
	step_fitted_1 = 5e-3;
	step_fitted_2 = 5e-3;
}//====================================================

/**
 *  \brief Reset any variables specific to an individual continuation process
 */
void NatParamEngine::cleanEngine(){
	ContinuationEngine::cleanEngine();
}//====================================================

}// End of astrohelion namespace


