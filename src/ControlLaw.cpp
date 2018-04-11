/**
 *  @file ControlLaw.cpp
 *  @author Andrew Cox
 *  @version March 3, 2017
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

#include <string>

#include "ControlLaw.hpp"
#include "Exceptions.hpp"
#include "SysData.hpp"

// temp include to do the law conversion
#include "ControlLaw_cr3bp_lt.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors
//-----------------------------------------------------------------------------

/**
 *  @brief Default constructor
 *  @details [long description]
 * 
 *  @param tp control law type; by default, set to NO_CTRL
 *  @param params vector of parameters used by the control law
 */
ControlLaw::ControlLaw(unsigned int tp, const std::vector<double> &params){
	this->lawType = tp;
	this->params = params;

	init();
}//====================================================

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

/**
 *  @brief Copy constructor
 * 
 *  @param law A constant reference to a ControlLaw object
 *  @return Assign this ControlLaw to be equal to the input ControlLaw
 */
ControlLaw& ControlLaw::operator =(const ControlLaw &law){
	copyMe(law);
	return *this;
}//====================================================

/**
 *  @brief Comparator
 * 
 *  @param lhs reference to left-hand-side ControlLaw
 *  @param rhs reference to right-hand-side ControlLaw
 * 
 *  @return Whether or not the two laws are equal
 */
bool operator ==(const ControlLaw &lhs, const ControlLaw &rhs){
	return lhs.lawType == rhs.lawType &&
		lhs.params == rhs.params;
}//====================================================

/**
 *  @brief Negative Comparator
 * 
 *  @param lhs reference to left-hand-side ControlLaw
 *  @param rhs reference to right-hand-side ControlLaw
 * 
 *  @return Whether or not the two laws are not equal
 */
bool operator !=(const ControlLaw &lhs, const ControlLaw &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve the control law ID
 *  @return the control law ID
 */
unsigned int ControlLaw::getType() const{ return lawType; }

/**
 *  @brief Retrieve the name of the law type as a string
 *  @return the name of the law type as a string
 */
std::string ControlLaw::getTypeString() const{
	return ControlLaw::typeToString(lawType);
}

/**
 *  @brief Retrieve the number of control variables that need to be included
 *  in an arcset propagation
 *  @return the number of control variables that need to be included
 *  in an arcset propagation
 */
unsigned int ControlLaw::getNumStates() const {return numStates; }

/**
 *  @brief Retrieve the number of outputs from the control Law
 *  @return the number of outputs from the control Law
 */
unsigned int ControlLaw::getNumOutputs() const {return numOutputs; }

/**
 *  @brief Retrieve the value of a parameter at the specified index.
 * 
 *  @param ix Index of the parameter. If negative, the index counts
 *  backwards from the end of the parameter vector, e.g., ix = -2
 *  returns the second-to-last element in params.
 *  @return the value stored in the parameter vector at the specified index.
 *  @throws Exception if the index is out of bounds
 */
double ControlLaw::getParam(int ix) const{
	int n = static_cast<int>(params.size());
	while(ix < 0){ ix += n; }

	if(ix >= n){
		char msg[128];
		sprintf(msg, "ControlLaw::getParam: Index ix = %d is larger"
			" than params vector (size %d)", ix, n);
		throw Exception(msg);
	}

	return params[ix];
}//====================================================

/**
 *  @brief Retrieve the parameter vector
 *  @return the parameter vector
 */
std::vector<double> ControlLaw::getParams() const { return params; }

/**
 *  @brief Retrieve a constant reference to the parameter vector
 *  @details Leverage this function when accessing parameters in 
 *  performance-sensitive operations
 *  @return a constant reference to the parameter vector
 */
const std::vector<double>& ControlLaw::getParamsRef_const() const {
	return params;
}

/**
 *  @brief Set the control law type
 * 
 *  @param id an ID number identifying the control law type
 */
void ControlLaw::setType(unsigned int id){
	lawType = id;
	init();
}//====================================================

/**
 *  @brief Set the value of a parameter at the specified index.
 * 
 *  @param ix Index of the parameter. If negative, the index counts
 *  backwards from the end of the parameter vector, e.g., ix = -2
 *  returns the second-to-last element in params.
 *  @param val the value to store in the parameter vector at the specified 
 *  index
 *  @throws Exception if the index is out of bounds
 */
void ControlLaw::setParam(int ix, double val){
	int n = static_cast<int>(params.size());
	while(ix < 0){ ix += n; }

	if(ix >= n){
		char msg[128];
		sprintf(msg, "ControlLaw::getParam: Index ix = %d is larger"
			" than params vector (size %d)", ix, n);
		throw Exception(msg);
	}

	params[ix] = val;
}//====================================================

/**
 *  @brief Set the control law parameters
 * 
 *  @param pParams array of parameters
 *  @param len number of elements in the parameter array
 */
void ControlLaw::setParams(const double *pParams, unsigned int len){
	std::copy(pParams, pParams+len, this->params.begin());
}//====================================================

/**
 *  @brief Set the control law parameters
 * 
 *  @param params reference to vector of parameters
 */
void ControlLaw::setParams(const std::vector<double> &params){
	this->params = params;
}//====================================================

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve the output of a control law.
 * 	@details For example, a set of accelerations incorporated into the velocity 
 * 	time derivative EOMs
 * 	
 *  @param t time parameter
 *  @param s state vector
 *  @param pSys system data object
 *  @param output initialized array of zeros in which to store the control 
 *  law outputs
 *  @param len number of elements in the `law` array
 *  
 *  @throws Exception if the control law ID, `lawType`, is not recognized
 */
void ControlLaw::getOutput(double t, const double *s, const SysData *pSys,
	double *output, unsigned int len) const{

	switch(lawType){
		case NO_CTRL:
			// Leave as a bunch of zeros
			break;
		default:
			throw Exception("ControlLaw::getLaw: Unrecognized lawType");
	}
	(void) t;
	(void) s;
	(void) pSys;
	(void) output;
	(void) len;
}//====================================================

/**
 *  @brief Retrieve the time derivatives of control states (if any exist)
 *  @details The behavior of this function defaults to returning
 *  all zeros for the derivatives, which is accurate for control laws that hold
 *  parameters constant along a control segment. By overriding this function 
 *  in derived classes, more specific behavior may be defined.
 * 
 *  @param t nondimensional integration time
 *  @param s integration state vector (from SimEngine / DynamicsModel)
 *  @param pSys pointer to the system data object
 *  @param deriv pointer to an array of zeros in which to store the 
 *  time-derivative of the control states
 *  @param len number of elements in the derivative array
 */
void ControlLaw::getTimeDeriv(double t, const double *s, const SysData *pSys,
	double *deriv, unsigned int len) const{

	switch(lawType){
		case NO_CTRL:
		default:
			// Leave as a bunch of zeros
			break;
	}
	(void) t;
	(void) s;
	(void) pSys;
	(void) deriv;
	(void) len;
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control state derivatives 
 *  with respect to all states (core states and control states).
 *  @details If a nontrivial set of control state derivatives exists, these 
 *  partial derivatives form the bottom set of rows of the A matrix
 * 
 *  @param t time parameter
 *  @param s state vector
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial 
 *  derivatives of the control state time derivatives
 *  @param len number of elements in the `partials` array
 */
void ControlLaw::getPartials_TimeDerivWRTAllState(double t, const double *s,
	const SysData *pSys, double *partials, unsigned int len) const{

	switch(lawType){
		case NO_CTRL:
		default:
			// Leave as zeros
			break;
	}
	(void) t;
	(void) s;
	(void) pSys;
	(void) partials;
	(void) len;
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the core state EOMs 
 *  with respect to the control states
 *  @details If a nontrivial set of control states exists, these partial
 *  derivatives form the right-hand block-column of the A matrix for the rows 
 *  associated with the core spacecraft state EOMs. I.e., these partial 
 *  derivatives do not include the partials of the control state time 
 *  derivatives w.r.t. the control states; those partial derivatives are 
 *  obtained from getPartials_TimeDerivWRTAllState()
 * 
 *  @param t time parameter
 *  @param s state vector
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial
 *  derivatives of the law output
 *  @param len number of elements in the `partials` array
 */
void ControlLaw::getPartials_EOMsWRTCtrlState(double t, const double *s,
	const SysData *pSys, double *partials, unsigned int len) const{

	switch(lawType){
		case NO_CTRL:
		default:
			// Leave as a bunch of zeros
			break;
	}
	(void) t;
	(void) s;
	(void) pSys;
	(void) partials;
	(void) len;
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control outputs with 
 *  respect to the core state variables
 *  @details Control outputs are leveraged in the EOMs, and the partial 
 *  derivatives of the control outputs w.r.t. core state variables are required 
 *  to form the linear relationship "A matrix."
 * 
 *  @param t time parameter
 *  @param s state vector
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial 
 *  derivatives of the control outputs
 *  @param len number of elements in the `partials` array
 */
void ControlLaw::getPartials_OutputWRTCoreState(double t, const double *s,
	const SysData *pSys, double *partials, unsigned int len) const{

	switch(lawType){
		case NO_CTRL:
		default:
			// Leave as a bunch of zeros
			break;
	}
	(void) t;
	(void) s;
	(void) pSys;
	(void) partials;
	(void) len;
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

void ControlLaw::print() const{
	printf("Control Law\n  Type = %s\n", typeToString(lawType).c_str());
	printf("  NumStates = %u\n  NumOutputs = %u\n", numStates, numOutputs);
	printf("  params = {");
	for(unsigned int i = 0; i < params.size(); i++){
		printf("%f", params[i]);
		if(i < params.size() - 1)
			printf(",  ");
	}
	printf("}\n");
}//====================================================

/**
 *  @brief Copy the ControlLaw object
 *  @param law reference to the source controlLaw
 */
void ControlLaw::copyMe(const ControlLaw &law){
	lawType = law.lawType;
	numStates = law.numStates;
	numOutputs = law.numOutputs;
	params = law.params;
}//====================================================

/**
 *  @brief Initialize all parameters that depend on the control law type
 *  @details [long description]
 */
void ControlLaw::init(){
	switch(lawType){
		case NO_CTRL:
			numStates = 0;
			numOutputs = 0;
			break;
		default: break;
	}
}//====================================================

/**
 *  @brief Retrieve a string that represents the law ID
 * 
 *  @param id control law ID
 *  @return a string that represents the law ID
 */
std::string ControlLaw::typeToString(unsigned int id){
	switch(id){
		case NO_CTRL: return "NONE";
		default: return "UNDEFINED";
	}
}//====================================================

/**
 * @brief Temporary funtion to convert low-thrust law IDs
 * @details [long description]
 * 
 * @param int [description]
 */
void ControlLaw::convertID(unsigned int &id){
	using L = ControlLaw_cr3bp_lt;
	switch(id){
		case 1:	// Const_FC_2D_LEFT
			id = L::CONST_F | L::CONST_C_2D | L::CSI_VAR_M | L::VEL_LEFT;
			break;
		case 2: // CONST_FC_2D_RIGHT
			id = L::CONST_F | L::CONST_C_2D | L::CSI_VAR_M | L::VEL_RIGHT;
			break;
		case 3: // CONST_F_PRO_VEL
			id = L::CONST_F | L::VEL_PT | L::CSI_VAR_M | L::VEL_PRO;
			break;
		case 4: // CONST_F_ANTI_VEL
			id = L::CONST_F | L::VEL_PT | L::CSI_VAR_M | L::VEL_ANTI;
			break;
		case 5: // CONST_F_GENERAL
			id = L::CONST_F | L::GENERAL | L::CSI_VAR_M;
			break;
		case 105: // CONST_MF_GENERAL
			id = L::CONST_F | L::GENERAL | L::CONST_M;
			break;
		case 1001: // VAR_F_CONST_C_2D_LEFT
			id = L::VAR_F_UBND | L::CONST_C_2D | L::CSI_VAR_M | L::VEL_LEFT;
			break;
		case 1002: // VAR_F_CONST_C_2D_RIGHT
			id = L::VAR_F_UBND | L::CONST_C_2D | L::CSI_VAR_M | L::VEL_RIGHT;
			break;
		case 1003: // VAR_F_PRO_VEL
			id = L::VAR_F_UBND | L::VEL_PT | L::CSI_VAR_M | L::VEL_PRO;
			break;
		case 1004: // VAR_F_ANTI_VEL
			id = L::VAR_F_UBND | L::VEL_PT | L::CSI_VAR_M | L::VEL_ANTI;
			break;
		case 1005: // VAR_F_GENERAL
			id = L::VAR_F_UBND | L::GENERAL | L::CSI_VAR_M;
			break;
		default: break;
	}
}

}// End of astrohelion namespace
