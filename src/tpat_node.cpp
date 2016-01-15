/**
 *  @file tpat_node.cpp
 *	@brief Data object that stores information about a node
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

#include "tpat_arc_step.hpp"
#include "tpat_exceptions.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a nodeset with a specified state and TOF
 *	@param state an array of 6 position and velocity states (non-dim)
 *	@param tof the non-dimensional time-of-flight
 */
tpat_node::tpat_node(const double *state, double tof){
	initArrays();
	std::copy(state, state+6, posVelState);
	extraParam[0] = tof;
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
	extraParam[0] = tof;
}//====================================================

/**
 *	@brief Copy constructor
 *	@param n a tpat_node reference
 */
tpat_node::tpat_node(const tpat_node &n) : tpat_arc_step(n) {}

/**
 *	@brief Create a node object from a generic arc
 *	step object.
 *
 *	This is permissible because tpat_node only defines
 *	new access methods, not new data objects
 *	@param s a tpat_arc_step reference
 */
tpat_node::tpat_node(const tpat_arc_step &s) : tpat_arc_step(s) {}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------


//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the time-of-flight for this node
 *	@return the time-of-flight for this node
 */
double tpat_node::getTOF() const { return extraParam[0]; }

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
	temp.insert(temp.end(), flags.begin(), flags.begin()+3);
	return temp;
}//====================================================

/**
 *	@brief set the time-of-flight
 *	@param t a non-dimensional time-of-flight
 */
void tpat_node::setTOF(double t){ extraParam[0] = t; }

/**
 *	@brief Set all velocity states to be continuous
 */
void tpat_node::setVel_AllCon(){
	flags[0] = true;
	flags[1] = true;
	flags[2] = true;
}//====================================================

/**
 *	@brief Set all velocity states to be discontinuous
 */
void tpat_node::setVel_AllDiscon(){
	flags[0] = false;
	flags[1] = false;
	flags[2] = false;
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param data a three-element boolean array. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 */
void tpat_node::setVelCon(bool data[3]){
	flags[0] = data[0];
	flags[1] = data[1];
	flags[2] = data[2];
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

	std::copy(data.begin(), data.begin()+3, flags.begin());
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
	flags[0] = xCon;
	flags[1] = yCon;
	flags[2] = zCon;
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize storange ararys
 *	
 *	Velocity continuities are all set to TRUE and the position and
 *	velocity states are set to NAN
 */
void tpat_node::initArrays(){
	tpat_arc_step::initArrays();

	for(int i = 0; i < 3; i++)
		flags.push_back(true);

	extraParam.assign(1,0);	// Definitely one spot for TOF, filled by constructor
}//====================================================
