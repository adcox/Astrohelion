/**
 *  @file Node.cpp
 *	@brief Stores information about a single node or integration state
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "Exceptions.hpp"
#include "Node.hpp"
#include "Utilities.hpp"

#include <cmath>


namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Default constructor
 */
Node::Node(){}

/**
 *  @brief Construct a node object
 * 
 *  @param state 6-element array of state variables
 *  @param epoch epoch associated with this node
 */
Node::Node(const double state[6], double epoch){
	std::copy(state, state+6, this->state);
	this->epoch = epoch;
}//====================================================

/**
 *  @brief Construct a node object
 * 
 *  @param state 6-element vector of state variables
 *  @param epoch epoch associated with this node
 *  @throw Exception if <tt>state</tt> does not have six elements
 */
Node::Node(std::vector<double> state, double epoch){
	if(state.size() != 6)
		throw Exception("Node::constructor: state must have six elements");

	std::copy(state.begin(), state.end(), this->state);
	this->epoch = epoch; 
}//====================================================

/**
 *  @brief Construct a node object
 * 
 *  @param state 6-element array of state variables
 *  @param accel 3-element array of acceleration values
 *  @param epoch epoch associated with this node
 */
Node::Node(const double state[6], const double accel[3], double epoch){
	std::copy(state, state+6, this->state);
	std::copy(accel, accel+3, this->accel);
	this->epoch = epoch;
}//====================================================

/**
 *  @brief Construct a node object
 * 
 *  @param state 6-element vector of state variables
 *  @param accel 3-element vector of acceleration values
 *  @param epoch epoch associated with this node
 *  @throw Exception if <tt>state</tt> does not have six elements
 *  @throw Exception if <tt>accel</tt> does not have three elements
 */
Node::Node(std::vector<double> state, std::vector<double> accel, double epoch){
	if(state.size() != 6)
		throw Exception("Node::constructor: state vector must have six elements");

	if(accel.size() != 3)
		throw Exception("Node::constructor: accel vector must have three elements");

	std::copy(state.begin(), state.end(), this->state);
	std::copy(accel.begin(), accel.end(), this->accel);
	this->epoch = epoch; 
}//====================================================

/**
 *  @brief Copy constructor
 * 
 *  @param n node object reference
 */
Node::Node(const Node &n) : Linkable(n){
	copyMe(n);
}//====================================================

/**
 *  @brief Destructor
 */
// Node::~Node(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param n a node object reference
 *	@return set this node equal to s and return *this
 */
Node& Node::operator =(const Node &n){
	Linkable::operator =(n);
	copyMe(n);
	return *this;
}//====================================================

/**
 *	@brief Determine if two nodes are identical
 *
 *	Conditions for identicalness:
 *	* Exact same state vector
 *	* (Not Active) Exact same extra parameter vector (e.g. epoch, tof, time, mass)
 * 	* (Not Active) Exact same flag vector (e.g. velocity continuity)
 *	If these conditions are met, acceleration and the STM should also
 *	be identical.
 *
 *	@return whether or not two nodes are identical
 */
bool operator ==(const Node &lhs, const Node &rhs){
	// Check state (implies accel is the same)
	for(int i = 0; i < 6; i++){
		if(lhs.state[i] != rhs.state[i])
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

	const Linkable link_lhs(lhs);
	const Linkable link_rhs(rhs);
	return link_lhs == link_rhs;
}//====================================================

/**
 *	@brief Determine if two nodes are different
 *	@return whether two nodes are different
 *	@see operator==
 */
bool operator != (const Node &lhs, const Node &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a constraint to the current set for this node
 *	@param c a new constraint
 */
void Node::addConstraint(Constraint c){
	cons.push_back(c);
}//====================================================

/**
 *	@brief Clear all constraints associated with this node
 */
void Node::clearConstraints(){ cons.clear(); }

/**
 *	@brief Remove the specified constraint
 *	@param ix the index of the constraint.
 *	@throw Exception if <tt>ix</tt> is out of bounds
 */
void Node::removeConstraint(int ix){
	if(ix < 0 || ix >= (int)(cons.size()))
		throw Exception("Node:removeConstraint: index out of bounds");
	cons.erase(cons.begin() + ix);
}//====================================================

/**
 *	@brief Set the list of constraints for this node
 *	@param constraints a vector of constraints
 */
void Node::setConstraints(std::vector<Constraint> constraints){
	cons = constraints;
}//====================================================

/**
 *	@brief Get a three-element vector containing the accelerations
 *	at this node
 *	@return a vector of accelerations (non-dimensional)
 */
std::vector<double> Node::getAccel() const{
	return std::vector<double>(accel, accel+3);
}//====================================================

/**
 *	@brief Get all constraints for this node
 *	@return a vector containing all constraints applied to this node
 */
std::vector<Constraint> Node::getConstraints() const{
	return cons;
}//====================================================

/**
 *  @brief Retrieve the epoch assocated with this node
 *  @return the epoch associated with this node, units consistent with the parent system
 */
double Node::getEpoch() const{ return epoch; }

/**
 *	@brief Access the value of the specified extra parameter
 *	@param ix the index of the parameter. If ix < 0, it will
 *	count backwards from the end of the array
 *	@return the value of the paramter associated with the 
 *	input index
 *	@throw Exception if <tt>ix</tt> is out of bounds
 */
double Node::getExtraParam(int ix) const {
	if(ix < 0)
		ix += extraParam.size();

	if(ix < 0 || ix >= (int)(extraParam.size())){
		astrohelion::printErr("Node::getExtraParam: Attempting to access index %d\n", ix);
		throw Exception("Node::getExtraParam: Cannot access extra param; index too high");
	}
	return extraParam[ix];
}//====================================================

/**
 *	@brief Get a vector containing all extra parameters for this node
 *	@return a vector containing all extra parameters for this node
 */
std::vector<double> Node::getExtraParams() const {
	return extraParam;
}//====================================================

/**
 *  @brief Retrieve the number of constraints stored by this object
 *  @return the number of constraints stored by this object
 */
int Node::getNumCons() const { return (int)(cons.size()); }

/**
 *	@brief Get the 6-element non-dimensional position and velocity state vector
 *	@return the 6-element non-dimensional position and velocity state vector
 */
std::vector<double> Node::getState() const {
	return std::vector<double>(state, state+6);
}//====================================================

/**
 *	@brief Set the acceleration vector for this node
 *	@param a a 3-element array of non-dimensional accelerations. Note
 *	that if the input array has fewer than three elements, un-initialized
 *	memory will be accessed
 */
void Node::setAccel(const double *a){
	std::copy(a, a+3, accel);
}//====================================================

/**
 *	@brief Set the acceleration vector for this node
 *	@param a a 3-element vector of non-dimensional accelerations
 *	@throw Exception if <tt>a</tt> does not have three elements
 */
void Node::setAccel(std::vector<double> a){
	if(a.size() != 3)
		throw Exception("Node::setAccel: input acceleration must have three elements");
	std::copy(a.begin(), a.begin()+3, accel);
}//====================================================

/**
 *  @brief Set the epoch associated with this node
 *  @param e the epoch, units consistent with parent system
 */
void Node::setEpoch(double e){ epoch = e; }

/**
 *  @brief Set the node number for all constraints stored in this node
 *  @param n The node number for the constraints; typically the index of this node
 */
void Node::setConstraintID(int n){
	for(size_t i = 0; i < cons.size(); i++){
		cons[i].setID(n);
	}
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
void Node::setExtraParam(int ix, double val){
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
 *	@brief Replace the extra parameter vector for this node
 *	@param p a new extra paremeter vector
 */
void Node::setExtraParams(std::vector<double> p){
	extraParam = p;
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s a 6-element array of non-dimensional position
 *	and velocity states. Note that if the input array has fewer
 *	than 6 states, un-initialized memory may be read.
 */
void Node::setState(const double *s){
	std::copy(s, s+6, state);
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s a 6-element vector of non-dimensional position
 *	and velocity states
 *	@throw Exception if <tt>s</tt> does not have six elements
 */
void Node::setState(std::vector<double> s){
	if(s.size() != 6)
		throw Exception("Node::setState: input vector must have six elements");
	std::copy(s.begin(), s.begin()+6, state);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy a node into this one
 *	@param n a node reference
 */
void Node::copyMe(const Node &n){
	std::copy(n.state, n.state+6, state);
	std::copy(n.accel, n.accel+3, accel);
	epoch = n.epoch;
	extraParam = n.extraParam;
	cons = n.cons;
	Linkable::copyMe(n);
}//====================================================

}// END of Astrohelion namespace