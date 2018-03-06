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
 *  Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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
 *  @param state array of state variables
 *  @param len length of the array
 *  @param epoch epoch associated with this node
 */
Node::Node(const double *state, unsigned int len, double epoch){
	this->state.assign(state, state + len);
	this->epoch = epoch;
}//====================================================

/**
 *  @brief Construct a node object
 * 
 *  @param state vector of state variables
 *  @param epoch epoch associated with this node
 */
Node::Node(std::vector<double> state, double epoch){
	this->state = state;
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
 *
 *	@return whether or not two nodes are identical
 */
bool operator ==(const Node &lhs, const Node &rhs){
	if(lhs.state.size() != rhs.state.size())
		return false;

	// Check state
	for(unsigned int i = 0; i < lhs.state.size(); i++){
		if(lhs.state[i] != rhs.state[i])
			return false;
	}

	// // Check extra parameters
	// if(lhs.extraParam.size() != rhs.extraParam.size())
	// 	return false;

	// for(unsigned int i = 0; i < lhs.extraParam.size(); i++){
	// 	if(lhs.extraParam[i] != rhs.extraParam[i])
	// 		return false;
	// }

	// // Check flags
	// if(lhs.flags.size() != rhs.flags.size())
	// 	return false;

	// for(unsigned int i = 0; i < lhs.flags.size(); i++){
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

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

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
 *	@throw Exception if `ix` is out of bounds
 */
void Node::removeConstraint(int ix){
	if(ix < 0 || ix >= static_cast<int>(cons.size()))
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
 *	@brief Get all constraints for this node
 *	@return a vector containing all constraints applied to this node
 */
std::vector<Constraint> Node::getConstraints() const{ return cons; }

/**
 * @brief Retrieve a reference to the constraint vector for in-place modifications
 * @return a reference to the vector containing all consraints applied to this node
 */
std::vector<Constraint>& Node::getConsRef() { return cons; }

/**
 *  @brief Retrieve the epoch assocated with this node
 *  @return the epoch associated with this node, units consistent with the parent system
 */
double Node::getEpoch() const{ return epoch; }

/**
 *  @brief Retrieve an extra parameter associated with the 
 *  specified key
 * 
 *  @param key A string that identifies an extra parameter
 *  @return the extra parameter value
 *  
 *  @throws Exception if the key is invalid
 */
double Node::getExtraParam(std::string key) const {
	if(extraParam.count(key) > 0){
		return extraParam.at(key);
	}else{
		throw Exception("Node::getExtraParam: Cannot access extra param; invalid key");
	}
}//====================================================

/**
 *	@brief Get a vector containing all extra parameters for this node
 *	@return a vector containing all extra parameters for this node
 */
std::map<std::string, double> Node::getExtraParams() const {
	return extraParam;
}//====================================================

/**
 *  @brief Retrieve an extra parameter vector associated with the
 *  specified key
 * 
 *  @param key String that identifies the extra parameter vector
 *  @return the extra parameter vector
 *  
 *  @throws Exception if the key is invalid
 */
std::vector<double> Node::getExtraParamVec(std::string key) const{
	if(extraParamVecs.count(key) > 0){
		return extraParamVecs.at(key);
	}else{
		throw Exception("Node::getExtraParamVec: Cannot access extra param vector; invalid key");
	}
}//====================================================

/**
 *  @brief Retrieve the extra parameter vector map
 *  @return the extra parameter vector map
 */
std::map<std::string, std::vector<double> > Node::getExtraParamVec() const{
	return extraParamVecs;
}//====================================================

/**
 *  @brief Retrieve the number of constraints stored by this object
 *  @return the number of constraints stored by this object
 */
unsigned int Node::getNumCons() const { return cons.size(); }

/**
 *  @brief Retrieve the event type that triggered the creation of this Node
 *  @return the event type that triggered the creation of this Node
 */
Event_tp Node::getTriggerEvent() const { return triggerEventTp; }

/**
 *	@brief Get the non-dimensional position and velocity state vector
 *	@return the non-dimensional position and velocity state vector
 */
std::vector<double> Node::getState() const { return state; }

/**
 *  @brief Retrieve a constant reference to the state
 *  @details This accessor is useful when computation time is important
 *  @return a constant reference to the state vector
 */
const std::vector<double>& Node::getStateRef_const() const{ return state; }

/**
 *  @brief Set the epoch associated with this node
 *  @param e the epoch, units consistent with parent system
 */
void Node::setEpoch(double e){ epoch = e; }

/**
 *  @brief Set an extra parameter value
 * 
 *  @param key A descriptive key identifying the extra parameter
 *  @param val value of the extra parameter
 */
void Node::setExtraParam(std::string key, double val){
	extraParam[key] = val;
}//====================================================

/**
 *	@brief Replace the extra parameter vector for this node
 *	@param p a new extra paremeter vector
 */
void Node::setExtraParams(std::map<std::string, double> p){
	extraParam = p;
}//====================================================

/**
 *  @brief Assign a vector to the designated key
 * 
 *  @param key string that identifies the extra parameter vector
 *  @param vec vector to assign to the key
 */
void Node::setExtraParamVec(std::string key, std::vector<double> vec){
	extraParamVecs[key] = vec;
}//====================================================

/**
 *  @brief Assign the entire extraParamVec map
 *  @param p A map that contains keys and extra parameter vectors
 */
void Node::setExtraParamVec(std::map<std::string, std::vector<double> > p){
	extraParamVecs = p;
}//====================================================

/**
 *  @brief Set the ID and also update the ID of any 
 *  associated constraints
 * 
 *  @param id ID that uniquely identifies the node
 */
void Node::setID(int id){
	Linkable::setID(id);
	for(unsigned int i = 0; i < cons.size(); i++){
		cons[i].setID(id);
	}
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s an array of non-dimensional position
 *	and velocity states
 *	@param len length of the array
 */
void Node::setState(const double *s, unsigned int len){
	state.assign(s, s+len);
}//====================================================

/**
 *	@brief Set the position-velocity state vector
 *	@param s a vector of non-dimensional position
 *	and velocity states
 */
void Node::setState(std::vector<double> s){ state = s; }

/**
 *  @brief Set the event that triggered the creation of this node
 *  @param tp Event type
 */
void Node::setTriggerEvent(Event_tp tp){ triggerEventTp = tp; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy a node into this one
 *	@param n a node reference
 */
void Node::copyMe(const Node &n){
	state = n.state;
	epoch = n.epoch;
	extraParam = n.extraParam;
	extraParamVecs = n.extraParamVecs;
	cons = n.cons;
	triggerEventTp = n.triggerEventTp;
	Linkable::copyMe(n);
}//====================================================

/**
 *  @brief Print a description of the node and its properties
 */
void Node::print() const{
	printf("Node | id = %d\n", ID);
	printf("\tEpoch = %.4f\n", epoch);
	printf("\tState =[");
	for(unsigned int s = 0; s < state.size(); s++){ printf("%.4f    ", state[s]); }
	printf("]\n");

	printf("\tExtra Parameters:\n");
	for(auto const& param : extraParam){
		printf("\t\t%s = %.4e\n", param.first.c_str(), param.second);
	}

	printf("\tExtra Parameter Vectors:\n");
	for(auto const& vec : extraParamVecs){
		printf("\t\t%s = [", vec.first.c_str());
		for(unsigned int i = 0; i < vec.second.size(); i++){
			printf("%.4f", vec.second[i]);
			if(i < vec.second.size() - 1)
				printf(",");
		}
		printf("]\n");
	}
}//====================================================

}// END of Astrohelion namespace