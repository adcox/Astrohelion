/**
 *	@file tpat_constraint.cpp
 *	@brief Data object that stores information about a node constraint
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_constraint.hpp"

#include <cmath>
#include <cstdio>

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Construct a constraint
 */
tpat_constraint::tpat_constraint(){
	data.clear();
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type
 *	@param type constraint type
 */
tpat_constraint::tpat_constraint(constraint_t type){
	this->type = type;
	data.clear();
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type and data values
 *	@param type constraint type
 *	@param node the node this constriant applies to
 *	@param data data vector (length n)
 */
tpat_constraint::tpat_constraint(constraint_t type, int node, std::vector<double> data){
	this->type = type;
	this->node = node;
	this->data = data;
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type, and data values
 *	@param type constraint type
 *	@param node the node this constriant applies to
 *	@param data data vector
 *	@param data_len the number of elements in d_len
 */
tpat_constraint::tpat_constraint(constraint_t type, int node, double* data, int data_len){
	this->type = type;
	this->node = node;
	this->data.reserve(data_len);

	for(int j = 0; j < data_len; j++){
		this->data.push_back(data[j]);
	}
}//============================================

/**
 *	@brief Create a copy of the specified constraint
 *	@param c a constraint
 */
tpat_constraint::tpat_constraint(const tpat_constraint& c){
	type = c.type;
	node = c.node;
	data = c.data;
}//============================================

/**
 *	@brief Destructor
 */
tpat_constraint::~tpat_constraint(){
	data.clear();
}

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param c a constraint
 *	@return this constraint, now equal to c
 */
tpat_constraint& tpat_constraint::operator =(const tpat_constraint& c){
	type = c.type;
	node = c.node;
	data = c.data;
	return *this;
}//============================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return what type of constraint this is
 */
tpat_constraint::constraint_t tpat_constraint::getType() const { return type; }

/**
 *	@return the node number associated with this constraint
 */
int tpat_constraint::getNode() const { return node; }

/**
 *	@return the data vector for this constraint
 */
std::vector<double> tpat_constraint::getData() const { return data; }

/**
 *	@return a count of the constrained states; certain constraint types, 
 *	like MATCH_CUST, give the option of constraining a subset of the entire
 *	node.
 */
int tpat_constraint::countConstrainedStates() const{
	int count = 0;
	for(int n = 0; n < ((int)data.size()); n++){
		count += !std::isnan(data.at(n));
	}
	return count;
}//============================================

/**
 *	@brief Set the constraint type
 *	@param t the type
 */
void tpat_constraint::setType(tpat_constraint::constraint_t t){ type = t; }

/**
 *	@brief Set the node index this constraint applies to
 *	@param n node index
 */
void tpat_constraint::setNode(int n){ node = n; }

/**
 *	Set the data for this node (should have nodeSize # elements)
 *	@param d the data, dimensions that match node dimensions
 */
void tpat_constraint::setData(std::vector<double> d){ data = d; }

/**
 *  @brief Set the data for this node
 * 
 *  @param dat an array of data values
 *  @param len number of elements in <tt>dat</tt>
 */
void tpat_constraint::setData(double *dat, int len){
	data.insert(data.begin(), dat, dat+len);
}

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@param t a constraint type
 *	@return a human-readable string representing a constraint type
 */
const char* tpat_constraint::getTypeStr(constraint_t t) const{
	switch(t){
		case NONE: { return "NONE"; break; }
		case STATE: { return "STATE"; break; }
		case MATCH_ALL: { return "MATCH_ALL"; break; }
		case MATCH_CUST: { return "MATCH_CUST"; break; }
		case DIST: { return "DIST"; break; }
		case MIN_DIST: { return "MIN_DIST"; break; }
		case MAX_DIST: { return "MAX_DIST"; break; }
		case MAX_DELTA_V: { return "MAX_DELTA_V"; break; }
		case DELTA_V: { return "DELTA_V"; break; }
		case JC: { return "JC"; break; }
		case SP: { return "SP"; break; }
		case SP_RANGE: { return "SP_RANGE"; break; }
		case TOF: { return "TOF"; break; }
		case APSE: {return "APS"; break; }
		case CONT_PV: {return "CONTINUOUS_POSITION_VELOCITY"; break; }
		case CONT_EX: { return "CONTINUOUS_EXTRA"; break; }
		case PSEUDOARC: { return "PSEUDO-ARCLENGTH"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//========================================

/**
 *	@brief Print this constraint and its data to the standard output.
 */
void tpat_constraint::print() const {
	printf("Constraint:\n  Type: %s\n  Node: %d\n  Data: ", getTypeStr(type), node);
	for(int n = 0; n < ((int)data.size()); n++){
		printf("%12.5f ", data[n]);
	}
	printf("\n");
}//========================================