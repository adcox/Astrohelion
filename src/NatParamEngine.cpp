/**
 * @file NatParamEngine.cpp
 * @brief
 * 
 * @author Andrew Cox
 * @version June 27, 2017
 * @copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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
#include "Arcset.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "Calculations.hpp"
#include "Exceptions.hpp"
#include "Family_PO.hpp"
#include "MultShootData.hpp"
#include "MultShootEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  @brief Default constructor
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
 *  @brief Get the number of solutions to store in the "curve fit memory"
 *  @details The number of solutions saved informs the least-squares
 *  curve-fitting algorithm that improves the performance of the natural 
 *  parameter continuation.
 * 
 *  @return the number of solutions to use when fitting a curve to predict
 *  the continuation variables
 */
unsigned int NatParamEngine::getCurveFitMem() const { return curveFitMem; }

/**
 *  @brief Get the number of simple continuation (i.e., no curve fitting)
 *  iterations to perform before switching to curve-fitted solutions
 * 
 *  @return the number of simple continuation iterations to perform
 */
unsigned int NatParamEngine::getNumSimple() const { return numSimple; }

/**
 *  @brief Get the number of iterations the corrector can take when a
 *  line search is NOT employed
 * 
 *  @return max number of iterations for the corrector
 */
unsigned int NatParamEngine::getNoSearchMaxIts() const { return noSearchMaxIts; }

/**
 *  @brief Get the number of iterations the corrector can take when a 
 *  line search IS employed
 * 
 *  @return max number of iterations for the corrector
 */
unsigned int NatParamEngine::getLineSearchMaxIts() const { return lineSearchMaxIts; }

/**
 *  @brief Get the slope threshold that is the limit between
 *  two step directions
 *  @details The natural parameter continuation typically leverages 
 *  several "independent" stepping variables to facilitate continuation 
 *  along a nonlinear contour. When a step in one variable results in
 *  large changes in the other, the stepping variable is changed to 
 *  limit the sensitivity.
 * 
 *  @return slope threshold between stepping strategies
 */
double NatParamEngine::getSlopeThresh() const { return slopeThresh; }

/**
 *  @brief Get the step size for simple (no curve-fitting) continuation
 *  @return the step size for simple (no curve-fitting) continuation
 */
double NatParamEngine::getStep_simple() const { return step_simple; }

/**
 *  @brief Get the step size for fitted continuation: variable 1
 *  @return the step size for fitted continuation: variable 1
 */
double NatParamEngine::getStep_fitted_1() const { return step_fitted_1; }

/**
 *  @brief Get the step size for fitted continuation: variable 2
 *  @return the step size for fitted continuation: variable 2
 */
double NatParamEngine::getStep_fitted_2() const { return step_fitted_2; }

/**
 *  @brief Set the number of solutions to store in the "curve fit memory"
 *  @details The number of solutions saved informs the least-squares
 *  curve-fitting algorithm that improves the performance of the natural 
 *  parameter continuation.
 * 
 *  @param mem Number of solutions to use when fitting a curve to predict
 *  the continuation variables
 */
void NatParamEngine::setCurveFitMem(unsigned int mem){ curveFitMem = mem; }

/**
 *  @brief Set the number of simple continuation (i.e., no curve fitting)
 *  iterations to perform before switching to curve-fitted solutions
 * 
 *  @param num the number of simple continuation iterations to perform
 */
void NatParamEngine::setNumSimple(unsigned int num){ numSimple = num; }

/**
 *  @brief Set the number of iterations the corrector can take when a
 *  line search is NOT employed
 * 
 *  @param num max number of iterations for the corrector
 */
void NatParamEngine::setNoSearchMaxIts(unsigned int num){ noSearchMaxIts = num; }

/**
 *  @brief Set the number of iterations the corrector can take when a 
 *  line search IS employed
 * 
 *  @param num max number of iterations for the corrector
 */
void NatParamEngine::setLineSearchMaxIts(unsigned int num){ lineSearchMaxIts = num; }

/**
 *  @brief Set the slope threshold that is the limit between
 *  two step directions
 *  @details The natural parameter continuation typically leverages 
 *  several "independent" stepping variables to facilitate continuation 
 *  along a nonlinear contour. When a step in one variable results in
 *  large changes in the other, the stepping variable is changed to 
 *  limit the sensitivity.
 * 
 *  @param thresh slope threshold between stepping strategies
 */
void NatParamEngine::setSlopeThresh(double thresh){ slopeThresh = thresh; }

/**
 *  @brief Set the step size for simple (no curve-fitting) continuation
 *  @param step the step size for simple (no curve-fitting) continuation
 */
void NatParamEngine::setStep_simple(double step){ step_simple = step; }

/**
 *  @brief Set the step size for fitted continuation: variable 1
 *  @param step the step size for fitted continuation: variable 1
 */
void NatParamEngine::setStep_fitted_1(double step){ step_fitted_1 = step; }

/**
 *  @brief Set the step size for fitted continuation: variable 2
 *  @param step the step size for fitted continuation: variable 2
 */
void NatParamEngine::setStep_fitted_2(double step){ step_fitted_2 = step; }

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Use natural parameter continuation to generate a family of periodic orbits in the CR3BP
 *  @details [long description]
 * 
 *  @param pFam pointer to a family in which the results are stored
 *  @param initGuess An initial guess for the first periodic orbit in the family
 *  @param alwaysFixStateVals A vector of initial state values that are always fixed
 *  @param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	@param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 */
void NatParamEngine::continuePO_cr3bp(Family_PO *pFam, const Arcset_cr3bp *initGuess, std::vector<double> alwaysFixStateVals,
	std::vector<unsigned int> indVarIx, std::vector<unsigned int> depVarIx){


	// Assume family is CR3BP
	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp*>(pFam->getSysData());

	Arcset_cr3bp guess(*initGuess), a1(pSys), a2(pSys);

	std::vector<double> perConData {0, 0, 0, 0, 0, NAN};
	Constraint perCon(Constraint_tp::MATCH_CUST, guess.getNodeRefByIx(-1).getID(), perConData);

	guess.addConstraint(perCon);
	std::vector<Arcset> allArcs {};
	continuePO(&guess, &a1, &a2, allArcs, alwaysFixStateVals, indVarIx, depVarIx);

	for(unsigned int i = 0; i < allArcs.size(); i++){
		pFam->addMember(allArcs[i]);
	}
}//====================================================

/**
 *	@brief Continue a family of periodic orbits via natural parameter continuation
 *	
 * 	@param pFam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	@param init_halfPerGuess a trajectory that is a good initial guess for half of 
 *	the first family member
 *	@param mirrorTypes a vector of variables that describe how the family mirrors in the rotating 
 *	reference frame. Each entry corresponds to an independent variable in `indVarIx`
 *	@param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	@param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 *
 *	@throws Exception if `indVarIx` has fewer than two elements
 *	@throws Exception if `mirrorTypes` does not have the same size as `indVarIx`
 *	@throws Exception if the eigenvalues of the monodromy matrix cannot be computed
 *	@throws Exception if one of the indices stored in `indVarIx` or `depVarIx` is
 *	out of range
 */
void NatParamEngine::continueSymmetricPO_cr3bp(Family_PO *pFam, const Arcset_cr3bp *init_halfPerGuess,
	std::vector<Mirror_tp> mirrorTypes, std::vector<unsigned int> indVarIx, std::vector<unsigned int> depVarIx){

	// This function makes some shortcuts for the symmetric periodic orbits and, thus, does not leverage
	// the more general continuePO() function

	if(pFam == nullptr)
		throw Exception("NatParamEngine::continueSymmetricPO_cr3bp: Cannot proceed with null family");

	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp*>(pFam->getSysData());

	if(indVarIx.size() < 2)
		throw Exception("NatParamEngine::cr3bp_natParamCont: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw Exception("NatParamEngine::cr3bp_natParamCont: there must be an equal number of ind. vars and mirror types");

	// Reset counters and storage containers
	printf("Beginning Natural Parameter Continuation\n");
	cleanEngine();
	bool diverged = false;

	Mirror_tp mirrorType = mirrorTypes[0];

	// Get info from the initial guess trajectory
	std::vector<double> IC = init_halfPerGuess->getStateByIx(0);
	double tof = init_halfPerGuess->getTotalTOF();
	double tof0 = tof;

	// Initially assume that we're fixing indVar1
	fixStates.push_back(indVarIx[0]);

	// Create a dummy arcset and create an iteration data object on the stack
	// The cr3bp_correctHalfPerSymPO() function will only pass an iteration data pointer back
	// if the one passed in is not nullptr, hence we create a valid object and delete it
	// before exiting the function
	Arcset_cr3bp halfPerGuess(*init_halfPerGuess), tempCorrected(pSys), halfPerCorrected(pSys);
	MultShootData *pItData = new MultShootData(&halfPerGuess);
	std::vector<Arcset> members {};

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
			perOrbit = cr3bp_correctHalfPerSymPO(&halfPerGuess, &tempCorrected, mirrorType, tol, pItData);

			diverged = false;
			halfPerCorrected = tempCorrected;
			printf("Orbit %03d converged!\n", static_cast<int>(members.size()));
		}catch(const DivergeException &e){
			diverged = true;
		}catch(const LinAlgException &e){
			printErr("There was a linear algebra error during family continuation...\n");
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
			printErr("Could not converge on the first family member; try a smaller step size\n");
			break;
		}else if(diverged && orbitCount > 0){
			perOrbit = members.back();	// Use previous solution as the converged solution to get correct next guess

			if(!decreaseStepSize())
				break;

		}else{	// Did not diverge

			// Save the computed orbit
			members.push_back(perOrbit);
			orbitCount++;

			increaseStepSize(pItData->count);

			// Compute eigenvalues
			MatrixXRd mono = perOrbit.getMonodromy();

			double monoErr = std::abs(1.0 - mono.determinant());
			if(monoErr > 1e-5)
				printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);

			// Add orbit to family
			pFam->addMember(perOrbit);
		}

		// Create next initial guess
		tof = perOrbit.getTimeByIx(-1);

		if(tof*tof0 < 0){
			printErr("Time-of-Flight changed sign: ending continuation process\n");
			break;
		}

		if(!updateIC(perOrbit, IC, tof, indVarIx, depVarIx, members))
			break;

		// Update the mirror type (not done in the general updateIC function)
		// *** Specific to symmetric continuation
		bool foundMatch = false;
		for(unsigned int v = 0; v < indVarIx.size(); v++){
			if(fixStates[0] == indVarIx[v]){
				mirrorType = mirrorTypes[v];
				foundMatch = true;
				break;
			}
		}
		if(!foundMatch)
			throw Exception("NatParamEngine::continueSymmetricPO_cr3bp: fixStates[0] is set to an unknown index");

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

/**
 *  @brief A general natural parameter continuation algorithm
 * 
 *  @param pInitGuess pointer to the initial guess. Any constraints applied to this initial arc will 
 *  be copied to all future family members and enforced during continuation
 *  @param pInput pointer to an arbitrary arcset object of the same model/system as the initial guess
 *  @param pConverged pointer to an arbitrary arcset object of the same model/system as the initial guess
 *  @param allArcs Storage vector for all converged family members
 *  @param alwaysFixStateVals a vector of initial state values that should always be fixed; a value of NAN indicates
 *  that the specified state should not be fixed.
 *  @param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	@param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 *	@param pEngineTemplate pointer to a multiple shooter engine that will be used, with
 *	some modifications, to correct orbits. You might want to control how times-of-flight
 *	are parameterized, whether or not to ignore crashes into the primaries, etc.
 */
void NatParamEngine::continuePO(const Arcset *pInitGuess, Arcset *pInput, Arcset *pConverged, 
	std::vector<Arcset>& allArcs, std::vector<double> alwaysFixStateVals,
	std::vector<unsigned int> indVarIx, std::vector<unsigned int> depVarIx,
	MultShootEngine *pEngineTemplate){

	if(pInitGuess == nullptr || pInput == nullptr || pConverged == nullptr)
		throw Exception("NatParamEngine::continuePO: cannot function with null input arcset pointers");

	const unsigned int core_dim = pInitGuess->getSysData()->getDynamicsModel()->getCoreStateSize();

	if(indVarIx.size() < 2)
		throw Exception("NatParamEngine::continuePO: Must specify two independent variables");

	if(alwaysFixStateVals.size() > core_dim)
		throw Exception("NatParamEngine::continuePO: alwaysFixStateVals vector cannot be longer than core state dimension");

	// Add extra entries to fill out the full core dimension, if not already accomplished
	std::vector<double> extraFill(core_dim - alwaysFixStateVals.size(), NAN);
	alwaysFixStateVals.insert(alwaysFixStateVals.end(), extraFill.begin(), extraFill.end());

	// Check to make sure the independent and dependent state indices are in bounds
	char msg[128];
	for(unsigned int i = 0; i < indVarIx.size(); i++){
		if(indVarIx[i] >= core_dim){
			sprintf(msg, "NatParamEngine::continuePO: indVarIx[%u] = %u is outside the core state dimension", i, indVarIx[i]);
			throw Exception(msg);
		}
	}
	for(unsigned int i = 0; i < depVarIx.size(); i++){
		if(depVarIx[i] >= core_dim){
			sprintf(msg, "NatParamEngine::continuePO: depVarIx[%u] = %u is outside the core state dimension", i, depVarIx[i]);
			throw Exception(msg);
		}
	}

	// Reset counters and storage containers
	printf("Beginning Natural Parameter Continuation\n");
	cleanEngine();
	bool diverged = false;

	// Initially assume that we're fixing indVar1
	fixStates.push_back(indVarIx[0]);

	// Get info from the initial guess trajectory
	std::vector<double> IC = pInitGuess->getStateByIx(0);
	double tof = pInitGuess->getTotalTOF();
	double tof0 = tof;

	MultShootData *pItData = new MultShootData(pInput);
	*pInput = *pInitGuess;	// Copy the initial guess

	// Create the multiple shooting engine object
	MultShootEngine msEngine;

	if(pEngineTemplate){
		// Use parameters from input engine
		msEngine = *pEngineTemplate;
	}else{
		// Set some default parameters
		msEngine.setIgnoreCrash(true);	// Ignore crashes into primary
		msEngine.setTOFType(MSTOF_tp::VAR_FIXSIGN);
	}

	// Always apply these 
	msEngine.setTol(tol);
	msEngine.setFullFinalProp(false);	// Accept minimal information from msEngine

	// Periodicity constraint(s)
	std::vector<Constraint> perCons = pInitGuess->getAllConstraints();

	while(orbitCount < numOrbits){
		pConverged->reset();	// Reset storage arcset for converged solution

		pInput->clearAllConstraints();		// Remove all constraints

		// Re-apply the periodicity constraints
		for(unsigned int c = 0; c < perCons.size(); c++){ pInput->addConstraint(perCons[c]); }

		// Fix specific initial states
		std::vector<double> initConData(core_dim, NAN);
		for(unsigned int i = 0; i < fixStates.size(); i++){
			if(std::isnan(alwaysFixStateVals[fixStates[i]]))
				initConData.at(fixStates[i]) = IC.at(fixStates[i]);
			else
				printWarn("Cannot fix state %u; it is permenantly constrained to be %f\n", fixStates[i], alwaysFixStateVals[fixStates[i]]);
		}

		for(unsigned int i = 0; i < alwaysFixStateVals.size(); i++){
			if(!std::isnan(alwaysFixStateVals[i]))
				initConData[i] = alwaysFixStateVals[i];
		}		

		Constraint initCon(Constraint_tp::STATE, pInput->getNodeRefByIx(0).getID(), initConData);
		pInput->addConstraint(initCon);

		try{
			printf("Guess for IC: [");
			for(unsigned int i = 0; i < core_dim; i++){ printf(" %7.4f", IC[i]); }
			printf("] tof = %.4f\n", tof);
			printf("Fix States: ");
			for(unsigned int i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			msEngine.setDoLineSearch(false);
			msEngine.setMaxIts(noSearchMaxIts);
			msEngine.multShoot(pInput, pConverged, pItData);
			
			diverged = false;
			printf("Orbit %03d converged!\n", static_cast<int>(allArcs.size()));
		}catch(DivergeException &de){
			try{
				printErr("Corrector diverged: %s\n", de.what());
				printf("Attempting to converge orbit %03d again with a line search for the step size\n", static_cast<int>(allArcs.size()));
				pConverged->reset();
				msEngine.setDoLineSearch(true);
				msEngine.setMaxIts(lineSearchMaxIts);

				msEngine.multShoot(pInput, pConverged, pItData);
				diverged = false;
				printf("Orbit %03d converged!\n", static_cast<int>(allArcs.size()));
			}catch(DivergeException &de2){
				printErr("Corrector diverged: %s\n", de2.what());
				diverged = true;	
			}catch(LinAlgException &lae2){
				printErr("There was a linear algebra error during family continuation: %s\n", lae2.what());
				break;
			}
		}catch(LinAlgException &lae){
			printErr("There was a linear algebra error during family continuation: %s\n", lae.what());
			break;
		}

		// Check for large changes in period to detect leaving family
		if(!diverged && orbitCount > 2){
			// difference in TOF; use abs() because corrector may employ reverse time and switch to forward time
			double dTOF = std::abs(pInput->getTotalTOF()) - std::abs(allArcs[allArcs.size()-1].getTotalTOF());
			double percChange = std::abs(dTOF/pInput->getTotalTOF());
			if(percChange > 0.25){
				printWarn("percChange = %.4f\n", percChange);
				printWarn("Period jumped (now = %.5f)! Left the family! Trying smaller step size...\n", pInput->getTotalTOF());
				diverged = true;
			}
		}

		if(diverged && orbitCount == 0){
			printErr("Could not converge on the first family member; try a smaller step size\n");
			break;
		}else if(diverged && orbitCount > 0){
			*pConverged = allArcs.back();	// Use previous solution as the converged solution to get correct next guess

			if(!decreaseStepSize())
				break;

		}else{	// Converged (did not diverge)
			allArcs.push_back(*pConverged);	// save the computed orbit
			orbitCount++;

			increaseStepSize(pItData->count);
		}

		// Create the next guess
		tof = pInput->getTotalTOF();

		if(tof*tof0 < 0){
			printErr("Time-of-Flight changed sign: ending continuation process\n");
			break;
		}

		if(!updateIC(*pConverged, IC, tof, indVarIx, depVarIx, allArcs))
			break;

		// Update the initial guess
		*pInput = *pConverged;
		pInput->setStateByIx(0, IC);

		// Change the TOF on the final segment
		// WARNING - might cause issues of TOF goes negative
		double dt = (tof - pInput->getTotalTOF())/(pInput->getNumSegs());
		double newTOF = 0;
		for(unsigned int s = 0; s < pInput->getNumSegs(); s++){
			newTOF = pInput->getTOFByIx(s) + dt;

			if(newTOF * pConverged->getTOFByIx(s) > 0){
				pInput->getSegRefByIx(s).setTOF(newTOF);
			}else{
				printErr("Did not update TOF on segment %u because it would change sign\n", s);
			}
		}
	}

	if(pItData){
		delete(pItData);
		pItData = nullptr;
	}
}//====================================================

/**
 *  @brief Decrease the step size
 *  @details The step size is decreased by a factor of <code>stepScaleFactor</code> if
 *  possible (without going lower than the minimum step size)
 *  
 *  @return Whether or not the step size was decreased. A return value of FALSE indicates that the step
 *  size is already at the minimum value.
 */
bool NatParamEngine::decreaseStepSize(){
	if(orbitCount <= numSimple){
		if(step_simple > minStepSize){
			step_simple = step_simple/stepScaleFactor > minStepSize ? step_simple/stepScaleFactor : minStepSize;
			printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", step_simple, minStepSize);
		}else{
			printErr("Minimum step size reached, could not converge... exiting\n");
			return false;
		}
	}else{
		double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

		if(dq > minStepSize){
			dq = dq/stepScaleFactor > minStepSize ? dq/stepScaleFactor : minStepSize;
			printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", dq, minStepSize);

			if(std::abs(indVarSlope) > slopeThresh)
				step_fitted_1 = dq;
			else
				step_fitted_2 = dq;
		}else{
			printErr("Minimum step size reached, could not converge... exiting\n");
			return false;
		}
	}

	return true;
}//====================================================

/**
 *  @brief Increase the step size based on how quickly the previous solution converged
 *  @details the step size is increased by the <code>stepScaleFactor</code> factor.
 * 
 *  @param its the number of iterations the shooting algorithm used to converge
 *  the previous solution.
 */
void NatParamEngine::increaseStepSize(unsigned int its){
	// Check to see if we should update the step size
	if(its < stepCount_increase && orbitCount > numSimple){
		double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

		if(dq < maxStepSize){
			dq = stepScaleFactor*dq < maxStepSize ? stepScaleFactor*dq : maxStepSize;
			printColor(MAGENTA, "  Increased step size to %0.4e (max = %.4e)\n", dq, maxStepSize);
			if(std::abs(indVarSlope) > slopeThresh)
				step_fitted_1 = dq;
			else
				step_fitted_2 = dq;
		}
	}
}//====================================================

/**
 *  @brief Update the initial condition vector and time-of-flight using vanilla natural parameter continuation
 *  or fitted natural parameter continuation
 *  @details [long description]
 * 
 *  @param convSoln the previous converged solution
 *  @param ic reference to the initial condition vector
 *  @param tof reference to the time-of-flight value
 *  @param indVarIx vector of the independent variable indices
 *  @param depVarIx vector of the dependent variable indices
 *  
 *  @return whether or not the update was successful.
 */
bool NatParamEngine::updateIC(const Arcset &convSoln, std::vector<double>& ic, 
	double& tof, const std::vector<unsigned int>& indVarIx, const std::vector<unsigned int>& depVarIx,
	const std::vector<Arcset>& members){

	// Compute the slope
	if(orbitCount >= 2){
		deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVarIx[0]] - members[orbitCount-2].getStateByIx(0)[indVarIx[0]];
		deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVarIx[1]] - members[orbitCount-2].getStateByIx(0)[indVarIx[1]];
		indVarSlope = deltaVar1/deltaVar2;
	}

	if(orbitCount < numSimple){
		// Use simple continuation; copy the converged IC, step forward in the independent variable
		ic = convSoln.getStateByIx(0);
		ic.at(indVarIx[0]) += step_simple;
	}else{

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
			// Use continuation in indVarIx[0]
			ic.at(indVarIx[0]) = convSoln.getStateByIx(0).at(indVarIx[0]) + astrohelion::sign(deltaVar1)*step_fitted_1;
			fixStates.clear();
			fixStates.push_back(indVarIx[0]);

			// Use Least Squares to predict dependent vars and unusued ind. vars
			allDepVars = depVarIx;
			if(std::find(depVarIx.begin(), depVarIx.end(), indVarIx[1]) == depVarIx.end()){
				allDepVars.push_back(indVarIx[1]);	// only add if it isn't already part of depVarIx
			}
			predictedIC = familyCont_LS(indVarIx[0], ic.at(indVarIx[0]), allDepVars, prevStates);
		}else{
			// Use continuation in indVarIx[1]
			ic.at(indVarIx[1]) = convSoln.getStateByIx(0).at(indVarIx[1]) + astrohelion::sign(deltaVar2)*step_fitted_2;
			fixStates.clear();
			fixStates.push_back(indVarIx[1]);

			// Use Least Squares to predict dependent vars and unusued ind. vars
			allDepVars = depVarIx;
			if(std::find(depVarIx.begin(), depVarIx.end(), indVarIx[0]) == depVarIx.end()){
				allDepVars.push_back(indVarIx[0]);	// only add if it isn't already part of depVarIx
			}
			predictedIC = familyCont_LS(indVarIx[1], ic.at(indVarIx[1]), allDepVars, prevStates);
		}

		// Update ic with predicted variables
		for(unsigned int n = 0; n < allDepVars.size(); n++){
			unsigned int ix = allDepVars[n];
			if(ix < 6)
				ic[ix] = predictedIC[ix];
			else if(ix == 6)
				tof = predictedIC[ix];
			else{
				printErr("NatParamEngine::updateIC: Cannot update independent variable; index out of range");
				return false;
			}
		}
	}

	return true;
}//====================================================

/**
 *  @brief use least squares to predict new values of variables in a continuation process
 *
 *  This function uses a 2nd-order polynomial fit to predict a set of
 *  dependent variables using past relationships between an independent
 *  variable and the dependent variables. The number of points considered
 *  in those past relationships can be adjusted to vary the "stiffness"
 *  of the fit.
 *
 *  Ocasionally the 2nd-order fit does not work well and we encounter a 
 *  singular (or very near singular) matrix. In this case, the algorithm
 *  will apply linear regression, which generally solves the problem and
 *  results in a safe inversion.
 *
 *  For all these inputs, the "state vector" must take the following form:
 *      [x, y, z, xdot, ydot, zdot, Period, Jacobi Constant]
 *
 *  @param indVarIx the index of the state variable to use as the indpendent variable
 *  @param nextInd The next value of the independent variable
 *  @param depVars a vector specifying the indices of the states that will be
 *  dependent variables. The algorithm will predict fugure values for these
 *  variables based on how they have changed with the independent variable.
 *  @param varHistory a vector representing an n x 8 matrix which contains
 *  information about previous states, period, and JC. n should be at least 
 *  3. If it is larger than MemSize, only the first set of MemSize rows will
 *  be used.
 *
 *  @return an 8-element vector with predictions for the dependent variables.
 *  If a particular variable has not been predicted, its place will be kept with
 *  a NAN value.
 *
 *  An example may make things more clear:
 *
 *  Say I am continuing a family and am using x as the natural parameter in the
 *  continuation. indVarIx would be 0 to represent x. We input the value of x
 *  for the next orbit in the family (nextInd) and specify which variables (from
 *  the 8-element "state") we would like to have predicted by least-squares.
 *  @throws Exception if the `varHistory` vector contains fewer than 
 *  8 elements or if the `depVars` vector has no data
 */
std::vector<double> NatParamEngine::familyCont_LS(unsigned int indVarIx, double nextInd, std::vector<unsigned int> depVars, std::vector<double> varHistory){
    const unsigned int STATE_SIZE = 8;
    const double EPS = 1e-14;

    if(varHistory.size() < STATE_SIZE)
        throw Exception("Calculations::familyCont_LS: Not enough data to create A matrix\n");

    if(depVars.size() == 0)
        throw Exception("Calculations::familyCont_LS: Not enough data to create B matrix\n");

    // Form A and B matrices
    std::vector<double> A_data;
    std::vector<double> B_data;

    for(unsigned int n = 0; n < varHistory.size()/STATE_SIZE; n++){
        double d = varHistory[n*STATE_SIZE + indVarIx];
        A_data.push_back(d*d);  // ind. var^2
        A_data.push_back(d);    // ind. var^1
        A_data.push_back(1);    // ind. var^0

        for(unsigned int p = 0; p < depVars.size(); p++){
            // vector of dependent variables
            B_data.push_back(varHistory[n*STATE_SIZE + depVars[p]]);
        }
    }

    MatrixXRd A = Eigen::Map<MatrixXRd>(&(A_data[0]), varHistory.size()/STATE_SIZE, 3);
    MatrixXRd B = Eigen::Map<MatrixXRd>(&(B_data[0]), varHistory.size()/STATE_SIZE, depVars.size());

    // Generate coefficient matrix; these are coefficients for second-order
    // polynomials in the new independent variable
    MatrixXRd G(A.cols(), A.cols());
    G.noalias() = A.transpose()*A;
    
    Eigen::JacobiSVD<MatrixXRd> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd.setThreshold(1e-14);
    Eigen::VectorXd S = svd.singularValues();
    double smallestVal = S.minCoeff();
    
    Eigen::RowVectorXd P(depVars.size());

    if(smallestVal > EPS){
        // Use 2nd-order polynomial fit; solve GC = A.transpose()*B for C
        MatrixXRd C = G.fullPivLu().solve(A.transpose()*B);
        
        double indMatData[] = {nextInd*nextInd, nextInd, 1};
        Eigen::RowVector3d indMat = Eigen::Map<Eigen::RowVector3d>(indMatData, 1, 3);

        P.noalias() = indMat*C;
    }else{
        // User 1st-order polynomial fit
        std::vector<double> A_lin_data;
        for(unsigned int n = 0; n < varHistory.size()/STATE_SIZE; n++){
            A_lin_data.push_back(varHistory[n*STATE_SIZE + indVarIx]);
            A_lin_data.push_back(1);
        }

        MatrixXRd A_lin = Eigen::Map<MatrixXRd>(&(A_lin_data[0]), varHistory.size()/STATE_SIZE, 2);
        Eigen::Matrix2d G;
        G.noalias() = A_lin.transpose()*A_lin;

        MatrixXRd C = G.fullPivLu().solve(A_lin.transpose()*B);

        Eigen::RowVector2d indMat(nextInd, 1);

        P.noalias() = indMat*C;
    }

    // Insert NAN for states that have not been predicted
    std::vector<double> predicted;
    predicted.assign(STATE_SIZE, NAN);
    for(unsigned int i = 0; i < depVars.size(); i++)
        predicted[depVars[i]] = P(i);

    return predicted;
}//====================================================
//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Make a copy of the natural parameter engine
 * 
 *  @param engine the source of the copy
 */
void NatParamEngine::copyMe(const NatParamEngine &engine){
	ContinuationEngine::copyMe(engine);
	curveFitMem = engine.curveFitMem;
	numSimple = engine.curveFitMem;
	slopeThresh = engine.slopeThresh;
	step_simple = engine.step_simple;
	step_fitted_1 = engine.step_fitted_1;
	step_fitted_2 = engine.step_fitted_2;

	orbitCount = engine.orbitCount;
	indVarSlope = engine.indVarSlope;
	deltaVar1 = engine.deltaVar1;
	deltaVar2 = engine.deltaVar2;

	fixStates = engine.fixStates;
}//====================================================

/**
 *  @brief Reset all continuation parameters
 */
void NatParamEngine::reset(){
	ContinuationEngine::reset();
	curveFitMem = 5;
	numSimple = 3;
	slopeThresh = 1;
	step_simple = 5e-4;
	step_fitted_1 = 5e-3;
	step_fitted_2 = 5e-3;
	lineSearchMaxIts = 100;
	noSearchMaxIts = 20;
}//====================================================

/**
 *  @brief Reset any variables specific to an individual continuation process
 */
void NatParamEngine::cleanEngine(){
	ContinuationEngine::cleanEngine();

	orbitCount = 0;
	indVarSlope = NAN;
	deltaVar1 = 1;
	deltaVar2 = 1;

	fixStates.clear();
}//====================================================

}// End of astrohelion namespace


