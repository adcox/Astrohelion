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

#include "tpat_exception.hpp"

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
 *	* Exact same extra parameter vector (e.g. epoch, tof, time, mass)
 * 	* Exact same flag vector (e.g. velocity continuity)
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

	// Check extra parameters
	if(lhs.extraParam.size() != rhs.extraParam.size())
		return false;

	for(size_t i = 0; i < lhs.extraParam.size(); i++){
		if(lhs.extraParam[i] != rhs.extraParam[i])
			return false;
	}

	// Check flags
	if(lhs.flags.size() != rhs.flags.size())
		return false;

	for(size_t i = 0; i < lhs.flags.size(); i++){
		if(lhs.flags[i] != rhs.flags[i])
			return false;
	}

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

std::vector<double> tpat_arc_step::getAccel() const{
	return std::vector<double>(accel, accel+3);
}//====================================================

double tpat_arc_step::getExtraParam(int ix) const {
	return extraParam.at(ix);
}//====================================================

std::vector<double> tpat_arc_step::getExtraParams() const {
	return extraParam;
}//====================================================

std::vector<double> tpat_arc_step::getPosVelState() const {
	return std::vector<double>(posVelState, posVelState+6);
}//====================================================

tpat_matrix tpat_arc_step::getSTM() const {
	return tpat_matrix(6,6,stm);
}//====================================================

std::vector<double> tpat_arc_step::getSTMElements() const {
	return std::vector<double>(stm, stm+36);
}//====================================================

void tpat_arc_step::setAccel(double *a){
	std::copy(a, a+3, accel);
}//====================================================

void tpat_arc_step::setAccel(std::vector<double> a){
	if(a.size() != 3)
		throw tpat_exception("tpat_arc_step::setAccel: input acceleration must have three elements");
	std::copy(a.begin(), a.begin()+3, accel);
}

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

void tpat_arc_step::setExtraParams(std::vector<double> p){
	extraParam = p;
}//====================================================

void tpat_arc_step::setPosVelState(double *s){
	std::copy(s, s+6, posVelState);
}//====================================================

void tpat_arc_step::setPosVelState(std::vector<double> s){
	if(s.size() != 6)
		throw tpat_exception("tpat_arc_step::setPosVelState: input vector must have six elements");
	std::copy(s.begin(), s.begin()+6, posVelState);
}//====================================================

void tpat_arc_step::setSTM(tpat_matrix m){
	if(m.getRows() != 6 || m.getCols() != 6)
		throw tpat_exception("tpat_arc_step::setSTM: STM must be 6x6");
	std::copy(m.getDataPtr(), m.getDataPtr()+36, stm);
}//====================================================

void tpat_arc_step::setSTM(double *elements){
	std::copy(elements, elements+36, stm);
}//====================================================

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


