/**
 *  @file tpat_node.cpp
 *
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

#include "tpat_node.hpp"

#include "tpat_exceptions.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Default constructor
 */
tpat_node::tpat_node() {
	initArrays();
}

/**
 *	@brief Create a nodeset with a specified state and TOF
 *	@param state an array of 6 position and velocity states (non-dim)
 *	@param tof the non-dimensional time-of-flight
 */
tpat_node::tpat_node(double *state, double tof){
	initArrays();
	std::copy(state, state+6, posVelState);
	this->tof = tof;
}//====================================================

/**
 *	@brief Construct a node from a state and time-of-flight
 *	@param state a 6-element vector of non-dimensional state values. 
 *	The vector can contain more than 6 elements, but only the first
 *	6 will be saved in the node
 *	@param tof the non-dimensional time of flight from this
 *	node to the next node. If this node is the last one in 
 *	a set, let the TOF be 0 or NAN.
 */
tpat_node::tpat_node(std::vector<double> state, double tof){
	if(state.size() < 6)
		throw tpat_exception("tpat_node: Cannot construct node from state with fewer than 6 elements");

	initArrays();
	std::copy(state.begin(), state.begin()+6, posVelState);
	this->tof = tof;
}//====================================================

/**
 *	@brief Copy constructor
 *	@param n a node reference
 */
tpat_node::tpat_node(const tpat_node &n){
	initArrays();
	copyMe(n);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_node::~tpat_node(){
	extraParam.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param n a node reference
 */
tpat_node& tpat_node::operator =(const tpat_node &n){
	copyMe(n);
	return *this;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs a node reference
 *	@param rhs a node reference
 *	@return true if the nodes are *exactly* the same. Numerical
 *	errors may cause this to evaluate to false even if two
 *	nodes are practically identical.
 */
bool operator ==(const tpat_node &lhs, const tpat_node &rhs){
	bool sameTOF = lhs.tof == rhs.tof;
	bool sameState = lhs.posVelState == rhs.posVelState;
	bool sameExtraParam = lhs.extraParam == rhs.extraParam;
	bool sameVelCon = lhs.velCon == rhs.velCon;

	return sameTOF && sameState && sameExtraParam && sameVelCon;
}//====================================================


//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the position and velocity states
 *	@return a vector containing three position states and
 *	three velocity states.
 */
std::vector<double> tpat_node::getPosVelState() const{
	std::vector<double> temp;
	temp.insert(temp.end(), posVelState, posVelState+6);
	return temp;
}//====================================================

/**
 *	@brief Retrieve one of the extra parameters for this node
 *	@return one of the extra parameters for this node
 */
double tpat_node::getExtraParam(int i) const { return extraParam.at(i); }

/**
 *	@brief Retrieve all extra parameters for this node
 *	@return all extra parameters for this node
 */
std::vector<double> tpat_node::getExtraParams() const { return extraParam; }

/**
 *	@brief Retrieve the time-of-flight for this node
 *	@return the time-of-flight for this node
 */
double tpat_node::getTOF() const { return tof; }

/**
 *	@brief Retrieve a vector describing which of the velocity states
 *	should be made continuous with the node before this one in a nodeset.
 *
 *	@return a 3-element boolean vector. The first element describes
 *	whether or not the x-velocity is continuous, the second describes
 *	the y-velocity continuity, etc.
 */
std::vector<bool> tpat_node::getVelCon() const {
	std::vector<bool> temp;
	temp.insert(temp.end(), velCon, velCon+3);
	return temp;
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
void tpat_node::setExtraParam(int ix, double val){
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
 *	@brief Set all extra parameters
 *	@param allExtra the new extra parameter vector
 */
void tpat_node::setExtraParams(std::vector<double> allExtra) {
	extraParam = allExtra;
}//====================================================

/**
 *	@brief Set the position and velocity states
 *	@param array a pointer to a 6-element array of non-dimensional states
 *
 *	If the input array has fewer than 6 elements, a memory problem
 *	will occur as the function will be reading from unknown memory blocks.
 */
void tpat_node::setPosVelState(double *array){
	std::copy(array, array+6, posVelState);
}//====================================================

/**
 *	@brief Set the position and velocity states
 *	@param vec a vector (size >= 6) containing the non-dimensional states
 */
void tpat_node::setPosVelState(std::vector<double> vec){
	if(vec.size() < 6)
		throw tpat_exception("tpat_node::setPosVelState: Cannot set state; input vector has too few elements");
	std::copy(vec.begin(), vec.begin()+6, posVelState);
}//====================================================

/**
 *	@brief set the time-of-flight
 *	@param t a non-dimensional time-of-flight
 */
void tpat_node::setTOF(double t){ tof = t; }

/**
 *	@brief Set all velocity states to be continuous
 */
void tpat_node::setVel_AllCon(){
	velCon[0] = true;
	velCon[1] = true;
	velCon[2] = true;
}//====================================================

/**
 *	@brief Set all velocity states to be discontinuous
 */
void tpat_node::setVel_AllDiscon(){
	velCon[0] = false;
	velCon[1] = false;
	velCon[2] = false;
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param data a three-element boolean array. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 */
void tpat_node::setVelCon(bool data[3]){
	std::copy(data, data+3, velCon);
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param data a three-element boolean vector. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 */
void tpat_node::setVelCon(std::vector<bool> data){
	if(data.size() < 3)
		throw tpat_exception("tpat_node::setVelCon: Need at least three velocity continuity booleans");

	std::copy(data.begin(), data.begin()+3, velCon);
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param xCon whether or not the x-velocity component should 
 *	be continuous with the node before this one
 *	@param yCon whether or not the y-velocity component should 
 *	be continuous with the node before this one
 *	@param zCon whether or not the z-velocity component should 
 *	be continuous with the node before this one
 */
void tpat_node::setVelCon(bool xCon, bool yCon, bool zCon){
	velCon[0] = xCon;
	velCon[1] = yCon;
	velCon[2] = zCon;
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy a node into this one
 *	@param n a node reference
 */
void tpat_node::copyMe(const tpat_node &n){
	std::copy(n.posVelState, n.posVelState+6, posVelState);
	std::copy(n.velCon, n.velCon+6, velCon);
	tof = n.tof;
	extraParam = n.extraParam;
}//====================================================

/**
 *	@brief Initialize storange ararys
 *	
 *	Velocity continuities are all set to TRUE and the position and
 *	velocity states are set to NAN
 */
void tpat_node::initArrays(){
	for(int i = 0; i < 6; i++)
		posVelState[i] = NAN;

	for(int i = 0; i < 3; i++)
		velCon[i] = true;
}//====================================================
