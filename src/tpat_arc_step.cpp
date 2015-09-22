/**
 *  @file tpat_arc_step.cpp
 *	@brief Data object that stores information about a single integration state
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tpat.hpp"

#include "tpat_arc_step.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Default constructor
 */
tpat_arc_step::tpat_arc_step(){
	initArrays();
}//====================================================

/**
 *	@brief Copy constructor
 *	@param s 
 */
tpat_arc_step::tpat_arc_step(const tpat_arc_step &s){
	copyMe(s);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_arc_step::~tpat_arc_step(){
	extraParam.clear();
	flags.clear();
	constraints.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param s an arc step object
 *	@return set this object equal to s and return *this
 */
tpat_arc_step& tpat_arc_step::operator =(const tpat_arc_step &s){
	copyMe(s);
	return *this;
}//====================================================

/**
 *	@brief Determine if two steps are identical
 *
 *	Conditions for identicalness:
 *	* Exact same state vector
 *	* (Not Active) Exact same extra parameter vector (e.g. epoch, tof, time, mass)
 * 	* (Not Active) Exact same flag vector (e.g. velocity continuity)
 *	If these conditions are met, acceleration and the STM should also
 *	be identical.
 *
 *	@return whether or not two steps are identical
 */
bool operator ==(const tpat_arc_step &lhs, const tpat_arc_step &rhs){
	// Check state (implies accel is the same)
	for(int i = 0; i < 6; i++){
		if(lhs.posVelState[i] != rhs.posVelState[i])
			return false;
	}

	// // Check extra parameters
	// if(lhs.extraParam.size() != rhs.extraParam.size())
	// 	return false;

	// for(size_t i = 0; i < lhs.extraParam.size(); i++){
	// 	if(lhs.extraParam[i] != rhs.extraParam[i])
	// 		return false;
	// }

	// // Check flags
	// if(lhs.flags.size() != rhs.flags.size())
	// 	return false;

	// for(size_t i = 0; i < lhs.flags.size(); i++){
	// 	if(lhs.flags[i] != rhs.flags[i])
	// 		return false;
	// }

	return true;
}//====================================================

/**
 *	@brief Determine if two steps are different
 *	@return whether two steps are different
 *	@see operator==
 */
bool operator != (const tpat_arc_step &lhs, const tpat_arc_step &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a constraint to the current set for this step
 *	@param c a new constraint
 */
void tpat_arc_step::addConstraint(tpat_constraint c){
	constraints.push_back(c);
}//====================================================

/**
 *	@brief Clear all constraints associated with this step
 */
void tpat_arc_step::clearConstraints(){ constraints.clear(); }

/**
 *	@brief Remove the specified constraint
 *	@param ix the index of the constraint. If the ix < 0, it will
 *	count backwards from the end of the set
 */
void tpat_arc_step::removeConstraint(int ix){
	if(ix < 0)
		ix += constraints.size();
	
	constraints.erase(constraints.begin() + ix);
}//====================================================

/**
 *	@brief Set the list of constraints for this arc step
 *	@param cons a vector of constraints
 */
void tpat_arc_step::setConstraints(std::vector<tpat_constraint> cons){
	constraints = cons;
}//====================================================

/**
 *	@brief Get a three-element vector containing the accelerations
 *	at this step
 *	@return a vector of accelerations (non-dimensional)
 */
std::vector<double> tpat_arc_step::getAccel() const{
	return std::vector<double>(accel, accel+3);
}//====================================================

/**
 *	@brief Get all constraints for this step
 *	@return a vector containing all constraints applied to this step
 */
std::vector<tpat_constraint> tpat_arc_step::getConstraints() const{
	return constraints;
}//====================================================

/**
 *	@brief Access the value of the specified extra parameter
 *	@param ix the index of the parameter. If ix < 0, it will
 *	count backwards from the end of the array
 *	@return the value of the paramter associated with the 
 *	input index
 */
double tpat_arc_step::getExtraParam(int ix) const {
	if(ix < 0)
		ix += extraParam.size();

	if(ix < 0 || ix >= (int)(extraParam.size())){
		printErr("tpat_arc_step::getExtraParam: Attempting to access index %d\n", ix);
		throw tpat_exception("tpat_arc_step::getExtraParam: Cannot access extra param; index too high");
	}
	return extraParam[ix];
}//====================================================

/**
 *	@brief Get a vector containing all extra parameters for this step
 *	@return a vector containing all extra parameters for this step
 */
std::vector<double> tpat_arc_step::getExtraParams() const {
	return extraParam;
}//====================================================

/**
 *	@brief Get the 6-element non-dimensional position and velocity state vector
 *	@return the 6-element non-dimensional position and velocity state vector
 */
std::vector<double> tpat_arc_step::getPosVelState() const {
	return std::vector<double>(posVelState, posVelState+6);
}//====================================================

/**
 *	@brief Retrieve the STM for this step
 *	@return the STM for this step
 */
tpat_matrix tpat_arc_step::getSTM() const {
	double el[36];
	std::copy(stm, stm+36, el);
	return tpat_matrix(6, 6, el);
}//====================================================

/**
 *	@brief Retrieve the STM elements
 *	@return a vector in row-major order containing the 36 STM elements
 */
std::vector<double> tpat_arc_step::getSTMElements() const {
	return std::vector<double>(stm, stm+36);
}//====================================================

/**
 *	@brief Set the acceleration vector for this step
 *	@param a a 3-element array of non-dimensional accelerations. Note
 *	that if the input array has fewer than three elements, un-initialized
 *	memory will be accessed
 */
void tpat_arc_step::setAccel(double *a){
	std::copy(a, a+3, accel);
}//====================================================

/**
 *	@brief Set the acceleration vector for this step
 *	@param a a 3-element vector of non-dimensional accelerations
 */
void tpat_arc_step::setAccel(std::vector<double> a){
	if(a.size() != 3)
		throw tpat_exception("tpat_arc_step::setAccel: input acceleration must have three elements");
	std::copy(a.begin(), a.begin()+3, accel);
}//====================================================

/**
 *	@brief Set a specific extra parameter to a specific value
 *	
 *	If the extra parameter vector isn't large enough to store the
 * 	new value, it will be enlarged with all emtpy spaces filled
 *	with NAN values.
 *
 *	@param ix the index of the extra parameter
 *	@param val the value of the extra paramter
 */
void tpat_arc_step::setExtraParam(int ix, double val){
	// Make the vector bigger if need be
	if((int)(extraParam.size()) <= ix){
		std::vector<double> temp = extraParam;
		extraParam.clear();
		extraParam.assign(ix+1, NAN);
		for(size_t n = 0; n < temp.size(); n++)
			extraParam[n] = temp[n];
	}

	// Put the desired value in the desired spot
	extraParam[ix] = val;
}//====================================================

/**
 *	@brief Replace the extra parameter vector for this step
 *	@param p a new extra paremeter vector
 */
void tpat_arc_step::setExtraParams(std::vector<double> p){
	extraParam = p;
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s a 6-element array of non-dimensional position
 *	and velocity states. Note that if the input array has fewer
 *	than 6 states, un-initialized memory may be read.
 */
void tpat_arc_step::setPosVelState(double *s){
	std::copy(s, s+6, posVelState);
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s a 6-element vector of non-dimensional position
 *	and velocity states
 */
void tpat_arc_step::setPosVelState(std::vector<double> s){
	if(s.size() != 6)
		throw tpat_exception("tpat_arc_step::setPosVelState: input vector must have six elements");
	std::copy(s.begin(), s.begin()+6, posVelState);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param m a 6x6 state transition matrix (non-dim)
 */
void tpat_arc_step::setSTM(tpat_matrix m){
	if(m.getRows() != 6 || m.getCols() != 6)
		throw tpat_exception("tpat_arc_step::setSTM: STM must be 6x6");
	std::copy(m.getDataPtr(), m.getDataPtr()+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements an array of stm elements in 
 *	row-major order. Note that the array must have at least 36
 *	elements, or un-initialized memory may be accessed
 */
void tpat_arc_step::setSTM(double *elements){
	std::copy(elements, elements+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements a 36-element vector of STM elements
 *	(non-dimensional) in row-major order
 */
void tpat_arc_step::setSTM(std::vector<double> elements){
	if(elements.size() != 36)
		throw tpat_exception("tpat_arc_step::setSTM: input vector must have 36 elements");
	std::copy(elements.begin(), elements.begin()+36, stm);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy a step into this one
 *	@param s an arc step reference
 */
void tpat_arc_step::copyMe(const tpat_arc_step &s){
	initArrays();
	std::copy(s.posVelState, s.posVelState+6, posVelState);
	std::copy(s.accel, s.accel+3, accel);
	std::copy(s.stm, s.stm+36, stm);
	extraParam = s.extraParam;
	flags = s.flags;
	constraints = s.constraints;
}//====================================================

/**
 *	@brief Initialize all storage arrays
 */
void tpat_arc_step::initArrays(){
	for(int i = 0; i < 6; i++)
		posVelState[i] = NAN;

	for(int i = 0; i < 3; i++)
		accel[i] = NAN;

	for(int i = 0; i < 36; i++)
		stm[i] = i % 7 == 0 ? 1 : 0;
}//====================================================


