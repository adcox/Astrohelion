/**
 *  \file ControlLaw.cpp
 *  \author Andrew Cox
 *  \version March 3, 2017
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

#include "ControlLaw.hpp"
#include "Exceptions.hpp"
#include "SysData.hpp"

namespace astrohelion{

//------------------------------------------------------------------------------------------------------
//      Constructors
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Default constructor
 *  \details [long description]
 * 
 *  \param id control law ID; by default, set to NO_CTRL
 *  \param params vector of parameters used by the control law
 */
ControlLaw::ControlLaw(unsigned int id, std::vector<double> params){
	this->lawID = id;
	this->params = params;

	init();
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Operators
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Copy constructor
 * 
 *  \param law A constant reference to a ControlLaw object
 *  \return Assign this ControlLaw to be equal to the input ControlLaw
 */
ControlLaw& ControlLaw::operator =(const ControlLaw &law){
	copyMe(law);
}//====================================================

/**
 *  \brief Comparator
 * 
 *  \param lhs reference to left-hand-side ControlLaw
 *  \param rhs reference to right-hand-side ControlLaw
 * 
 *  \return Whether or not the two laws are equal
 */
bool operator ==(const ControlLaw &lhs, const ControlLaw &rhs){
	return lhs.lawID == rhs.lawID &&
		lhs.params == rhs.params;
}//====================================================

/**
 *  \brief Negative Comparator
 * 
 *  \param lhs reference to left-hand-side ControlLaw
 *  \param rhs reference to right-hand-side ControlLaw
 * 
 *  \return Whether or not the two laws are not equal
 */
bool operator !=(const ControlLaw &lhs, const ControlLaw &rhs){
	return !(lhs == rhs);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Set and Get Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the control law ID
 *  \return the control law ID
 */
unsigned int ControlLaw::getLawID() const{ return lawID; }

/**
 *  \brief Retrieve the number of control variables that need to be included
 *  in an arcset propagation
 *  \return the number of control variables that need to be included
 *  in an arcset propagation
 */
unsigned int ControlLaw::getNumStates() const {return numStates; }

/**
 *  \brief Retrieve the parameter vector
 *  \return the parameter vector
 */
std::vector<double> ControlLaw::getParams() const { return params; }

/**
 *  \brief Retrieve a constant reference to the parameter vector
 *  \details Leverage this function when accessing parameters in performance-sensitive operations
 *  \return a constant reference to the parameter vector
 */
const std::vector<double>& ControlLaw::getParamsRef_const() const { return params; }

void ControlLaw::setLawID(unsigned int id){
	lawID = id;
	init();
}//====================================================

void ControlLaw::setParams(double *params, unsigned int len){
	std::copy(params, params+len, this->params.begin());
}//====================================================

void ControlLaw::setParams(std::vector<double> params){
	this->params = params;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Switchboard Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the output of a control law
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param law empty, initialized array to store the control law output in
 *  \param len number of elements in the <tt>law</tt> array
 *  
 *  \throws Exception if the control law ID, <tt>lawID</tt>, is not recognized
 */
void ControlLaw::getLaw(double t, const double *s, const SysData *pSys, double *law, unsigned int len) const{
	switch(lawID){
		case NO_CTRL:
			// Handle default case with no control
			for(unsigned int i = 0; i < len; i++){
				law[i] = 0;
			}
			break;
		default:
			throw Exception("ControlLaw::getLaw: Unrecognized lawID");
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

/**
 *  \brief Retrieve the derivative of the control law
 *  \details The behavior of this function defaults to returning
 *  all zeros for the derivatives, which is accurate for control laws that hold
 *  parameters constant along a thrust arc. By overriding this function in derived
 *  classes, more specific behavior may be defined.
 * 
 *  \param t nondimensional integration time
 *  \param s integration state vector (from SimEngine / DynamicsModel)
 *  \param pSys pointer to the system data object
 *  \param deriv pointer to an array in which to store the derivative of the law
 *  \param int number of elements in the derivative array
 */
void ControlLaw::getLaw_deriv(double t, const double *s, const SysData *pSys, double *deriv, unsigned int len) const{
	switch(lawID){
		case NO_CTRL:
		default:
			// Handle default case with no control
			for(unsigned int i = 0; i < len; i++){
				deriv[i] = 0;
			}
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

/**
 *  \brief Retrieve the partial derivatives of the control law with respect to state variables
 *  \details A set of partial derivatives of the control law outputs are computed with respect to the 
 *  states at the given time, state, in the specified system
 * 
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param partials empty, initialized array to store the control law derivatives in
 *  \param len number of elements in the <tt>law</tt> array
 *  
 *  \throws Exception if the control law ID, <tt>lawID</tt>, is not recognized
 */
void ControlLaw::getPartials_State(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const{
	switch(lawID){
		case NO_CTRL:
			// Handle default case with no control
			for(unsigned int i = 0; i < len; i++){
				partials[i] = 0;
			}
			break;
		default:
			throw Exception("ControlLaw::GetPartials_State: Unrecognized lawID");
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Copy the ControlLaw object
 *  \param law reference to the source controlLaw
 */
void ControlLaw::copyMe(const ControlLaw &law){
	lawID = law.lawID;
	numStates = law.numStates;
	params = law.params;
}//====================================================

/**
 *  \brief Initialize all parameters that depend on the control law type
 *  \details [long description]
 */
void ControlLaw::init(){
	switch(lawID){
		case NO_CTRL:
			numStates = 0;
			break;
		default:
			throw Exception("ControlLaw::init: Unrecognized lawID");
	}
}//====================================================

/**
 *  \brief Retrieve a string that represents the law ID
 * 
 *  \param id control law ID
 *  \return a string that represents the law ID
 */
std::string ControlLaw::lawIDToString(unsigned int id) const{
	switch(id){
		case NO_CTRL: return "NONE";
		default: return "UNDEFINED";
	}
}//====================================================

}