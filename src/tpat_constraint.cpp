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
}//====================================================

/**
 *	@brief Construct a constraint with specified constraint type
 *	@param type constraint type
 */
tpat_constraint::tpat_constraint(tpat_constraint_tp type){
	this->type = type;
	data.clear();
	setAppType();
}//====================================================

/**
 *	@brief Construct a constraint with specified constraint type and data values
 *	@param type constraint type
 *	@param id the ID of the object this constriant applies to
 *	@param data data vector (length n)
 */
tpat_constraint::tpat_constraint(tpat_constraint_tp type, int id, std::vector<double> data){
	this->type = type;
	this->id = id;
	this->data = data;
	setAppType();
}//====================================================

/**
 *	@brief Construct a constraint with specified constraint type, and data values
 *	@param type constraint type
 *	@param id the ID of the object this constriant applies to
 *	@param data data vector
 *	@param data_len the number of elements in d_len
 */
tpat_constraint::tpat_constraint(tpat_constraint_tp type, int id, const double* data, int data_len){
	this->type = type;
	this->id = id;
	this->data.insert(this->data.begin(), data, data + data_len);
	setAppType();
}//====================================================

/**
 *	@brief Create a copy of the specified constraint
 *	@param c a constraint
 */
tpat_constraint::tpat_constraint(const tpat_constraint& c){
	copyMe(c);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_constraint::~tpat_constraint(){}

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param c a constraint
 *	@return this constraint, now equal to c
 */
tpat_constraint& tpat_constraint::operator =(const tpat_constraint& c){
	copyMe(c);
	return *this;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Retrieve the application type for this constraint, i.e., what type of
 *  objects it controls and can be applied to
 *  @return the application type
 */
tpat_constraint::tpat_conApp_tp tpat_constraint::getAppType() const{ return appType; }

/**
 *	@return what type of constraint this is
 */
tpat_constraint::tpat_constraint_tp tpat_constraint::getType() const { return type; }

/**
 *	@return the ID that identifies the object constrained by this constraint
 */
int tpat_constraint::getID() const { return id; }

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
}//====================================================

/**
 *	@brief Set the constraint type
 *	@param t the type
 */
void tpat_constraint::setType(tpat_constraint::tpat_constraint_tp t){
	type = t;
	setAppType(); 	// Update, if necessary
}//====================================================

/**
 *	@brief Set the object ID this constraint applies to
 *	@param n the ID
 */
void tpat_constraint::setID(int n){ id = n; }

/**
 *	Set the data for this id (should have nodeSize # elements)
 *	@param d the data, dimensions that match node dimensions
 */
void tpat_constraint::setData(std::vector<double> d){ data = d; }

/**
 *  @brief Set the data for this constraint
 * 
 *  @param dat an array of data values
 *  @param len number of elements in <tt>dat</tt>
 */
void tpat_constraint::setData(const double *dat, int len){
	data.clear();
	data.insert(data.begin(), dat, dat+len);
}//====================================================

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *  @brief Copy the constraint
 * 
 *  @param c reference to a constraint object
 */
void tpat_constraint::copyMe(const tpat_constraint &c){
	appType = c.appType;
	type = c.type;
	id = c.id;
	data = c.data;
}//====================================================

/**
 *  @brief Get a human-redable string representing the constraint type
 *  @return a human-redable string representing the constraint type
 */
const char* tpat_constraint::getTypeStr() const{ return getConTypeStr(type); }

/**
 *	@param t a constraint type
 *	@return a human-readable string representing a constraint type
 */
const char* tpat_constraint::getConTypeStr(tpat_constraint::tpat_constraint_tp t){
	switch(t){
		case tpat_constraint::NONE: { return "NONE"; break; }
		case tpat_constraint::STATE: { return "STATE"; break; }
		case tpat_constraint::MATCH_ALL: { return "MATCH_ALL"; break; }
		case tpat_constraint::MATCH_CUST: { return "MATCH_CUST"; break; }
		case tpat_constraint::DIST: { return "DIST"; break; }
		case tpat_constraint::MIN_DIST: { return "MIN_DIST"; break; }
		case tpat_constraint::MAX_DIST: { return "MAX_DIST"; break; }
		case tpat_constraint::MAX_DELTA_V: { return "MAX_DELTA_V"; break; }
		case tpat_constraint::DELTA_V: { return "DELTA_V"; break; }
		case tpat_constraint::JC: { return "JC"; break; }
		case tpat_constraint::SP: { return "SP"; break; }
		case tpat_constraint::SP_RANGE: { return "SP_RANGE"; break; }
		case tpat_constraint::SP_DIST: { return "SP_DIST"; break; }
		case tpat_constraint::SP_MAX_DIST: { return "SP_MAX_DIST"; break; }
		case tpat_constraint::TOF: { return "TOF"; break; }
		case tpat_constraint::APSE: {return "APSE"; break; }
		case tpat_constraint::CONT_PV: {return "CONTINUOUS SEG2NODE POSITION_VELOCITY"; break; }
		case tpat_constraint::CONT_EX: { return "CONTINUOUS SEG2NODE EXTRA"; break; }
		case tpat_constraint::SEG_CONT_PV: { return "CONTINUOUS SEG2SEG POSITION_VELOCITY"; break; }
		case tpat_constraint::SEG_CONT_EX: { return "CONTINUOUS SEG2SEG EXTRA";  break; }
		case tpat_constraint::PSEUDOARC: { return "PSEUDO-ARCLENGTH"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//====================================================

/**
 *  @param t a constraint application type
 *  @return a human-readable string representing an application type
 */
const char* tpat_constraint::getAppTypeStr(tpat_constraint::tpat_conApp_tp t){
	switch(t){
		case tpat_constraint::APP_TO_NODE: return "Nodes"; break;
		case tpat_constraint::APP_TO_ARC: return "Whole Arcset"; break;
		case tpat_constraint::APP_TO_SEG: return "Segments"; break;
		default: return "Undefined!"; break;
	}
}//====================================================

/**
 *	@brief Print this constraint and its data to the standard output.
 */
void tpat_constraint::print() const {
	printf("Constraint:\n  Type: %s\n  Applies to: %s (ID %d)\n  Data: ",
		getConTypeStr(type), getAppTypeStr(appType), id);
	for(int n = 0; n < ((int)data.size()); n++){
		printf("%12.5f ", data[n]);
	}
	printf("\n");
}//====================================================

/**
 *  @brief Set the application type appropriately based on the constraint type
 */
void tpat_constraint::setAppType(){
	switch(type){
		case tpat_constraint::PSEUDOARC:
		case tpat_constraint::TOF:
		case tpat_constraint::MAX_DELTA_V:
		case tpat_constraint::DELTA_V:
		case tpat_constraint::SEG_CONT_PV:
		case tpat_constraint::SEG_CONT_EX:
			appType = tpat_constraint::APP_TO_ARC;
			break;
		case tpat_constraint::CONT_PV:
		case tpat_constraint::CONT_EX:
			appType = tpat_constraint::APP_TO_SEG;
			break;
		default:
			appType = tpat_constraint::APP_TO_NODE;
			break;
	}
}//====================================================