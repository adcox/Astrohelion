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
TPAT_Constraint::TPAT_Constraint(){
	data.clear();
}//====================================================

/**
 *	@brief Construct a constraint with specified constraint type
 *	@param type constraint type
 */
TPAT_Constraint::TPAT_Constraint(TPAT_Constraint_Tp type){
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
TPAT_Constraint::TPAT_Constraint(TPAT_Constraint_Tp type, int id, std::vector<double> data){
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
TPAT_Constraint::TPAT_Constraint(TPAT_Constraint_Tp type, int id, const double* data, int data_len){
	this->type = type;
	this->id = id;
	this->data.insert(this->data.begin(), data, data + data_len);
	setAppType();
}//====================================================

/**
 *	@brief Create a copy of the specified constraint
 *	@param c a constraint
 */
TPAT_Constraint::TPAT_Constraint(const TPAT_Constraint& c){
	copyMe(c);
}//====================================================

/**
 *	@brief Destructor
 */
TPAT_Constraint::~TPAT_Constraint(){}

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param c a constraint
 *	@return this constraint, now equal to c
 */
TPAT_Constraint& TPAT_Constraint::operator =(const TPAT_Constraint& c){
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
TPAT_ConApp_Tp TPAT_Constraint::getAppType() const{ return appType; }

/**
 *	@return what type of constraint this is
 */
TPAT_Constraint_Tp TPAT_Constraint::getType() const { return type; }

/**
 *	@return the ID that identifies the object constrained by this constraint
 */
int TPAT_Constraint::getID() const { return id; }

/**
 *	@return the data vector for this constraint
 */
std::vector<double> TPAT_Constraint::getData() const { return data; }

/**
 *  @brief Retrieve the first data value that is not an NAN
 *  @return The value of the first data value. If no data
 *  value is located, NAN is returned.
 *  @see getFirstDataValue(int*)
 */
double TPAT_Constraint::getFirstDataValue() const{
	int ix = 0;
	return getFirstDataValue(&ix);
}//====================================================

/**
 *  @brief Retrieve the first data value that is not an NAN
 * 
 *  @param ix The index of the first data value will be
 *  stored in this integer. If no data value is located,
 *  ix will be set to -1.
 *  @return The value of the first data value. If no data
 *  value is located, NAN is returned.
 */
double TPAT_Constraint::getFirstDataValue(int *ix) const{
	for(size_t i = 0; i < data.size(); i++){
		if(!std::isnan(data[i])){
			*ix = (int)i;
			return data[i];
		}
	}

	*ix = -1;
	return NAN;
}//====================================================

/**
 *	@return a count of the constrained states; certain constraint types, 
 *	like TPAT_Constraint_Tp::MATCH_CUST, give the option of constraining a subset of the entire
 *	node.
 */
int TPAT_Constraint::countConstrainedStates() const{
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
void TPAT_Constraint::setType(TPAT_Constraint_Tp t){
	type = t;
	setAppType(); 	// Update, if necessary
}//====================================================

/**
 *	@brief Set the object ID this constraint applies to
 *	@param n the ID
 */
void TPAT_Constraint::setID(int n){ id = n; }

/**
 *	Set the data for this id (should have nodeSize # elements)
 *	@param d the data, dimensions that match node dimensions
 */
void TPAT_Constraint::setData(std::vector<double> d){ data = d; }

/**
 *  @brief Set the data for this constraint
 * 
 *  @param dat an array of data values
 *  @param len number of elements in <tt>dat</tt>
 */
void TPAT_Constraint::setData(const double *dat, int len){
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
void TPAT_Constraint::copyMe(const TPAT_Constraint &c){
	appType = c.appType;
	type = c.type;
	id = c.id;
	data = c.data;
}//====================================================

/**
 *  @brief Get a human-redable string representing the constraint type
 *  @return a human-redable string representing the constraint type
 */
const char* TPAT_Constraint::getTypeStr() const{ return getConTypeStr(type); }

/**
 *	@param t a constraint type
 *	@return a human-readable string representing a constraint type
 */
const char* TPAT_Constraint::getConTypeStr(TPAT_Constraint_Tp t){
	switch(t){
		case TPAT_Constraint_Tp::NONE: { return "NONE"; break; }
		case TPAT_Constraint_Tp::STATE: { return "TPAT_Constraint_Tp::STATE"; break; }
		case TPAT_Constraint_Tp::MATCH_ALL: { return "TPAT_Constraint_Tp::MATCH_ALL"; break; }
		case TPAT_Constraint_Tp::MATCH_CUST: { return "TPAT_Constraint_Tp::MATCH_CUST"; break; }
		case TPAT_Constraint_Tp::DIST: { return "TPAT_Constraint_Tp::DIST"; break; }
		case TPAT_Constraint_Tp::MIN_DIST: { return "TPAT_Constraint_Tp::MIN_DIST"; break; }
		case TPAT_Constraint_Tp::MAX_DIST: { return "TPAT_Constraint_Tp::MAX_DIST"; break; }
		case TPAT_Constraint_Tp::MAX_DELTA_V: { return "TPAT_Constraint_Tp::MAX_DELTA_V"; break; }
		case TPAT_Constraint_Tp::DELTA_V: { return "TPAT_Constraint_Tp::DELTA_V"; break; }
		case TPAT_Constraint_Tp::JC: { return "JC"; break; }
		case TPAT_Constraint_Tp::SP: { return "SP"; break; }
		case TPAT_Constraint_Tp::SP_RANGE: { return "SP_RANGE"; break; }
		case TPAT_Constraint_Tp::SP_DIST: { return "SP_DIST"; break; }
		case TPAT_Constraint_Tp::SP_MAX_DIST: { return "SP_MAX_DIST"; break; }
		case TPAT_Constraint_Tp::TOF: { return "TOF"; break; }
		case TPAT_Constraint_Tp::APSE: {return "APSE"; break; }
		case TPAT_Constraint_Tp::CONT_PV: {return "CONTINUOUS TPAT_Constraint_Tp::SEG_2NODE POSITION_VELOCITY"; break; }
		case TPAT_Constraint_Tp::CONT_EX: { return "CONTINUOUS TPAT_Constraint_Tp::SEG_2NODE EXTRA"; break; }
		case TPAT_Constraint_Tp::SEG_CONT_PV: { return "CONTINUOUS TPAT_Constraint_Tp::SEG_2TPAT_Constraint_Tp::SEG_ POSITION_VELOCITY"; break; }
		case TPAT_Constraint_Tp::SEG_CONT_EX: { return "CONTINUOUS TPAT_Constraint_Tp::SEG_2TPAT_Constraint_Tp::SEG_ EXTRA";  break; }
		case TPAT_Constraint_Tp::PSEUDOARC: { return "PSEUDO-ARCLENGTH"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//====================================================

/**
 *  @param t a constraint application type
 *  @return a human-readable string representing an application type
 */
const char* TPAT_Constraint::getAppTypeStr(TPAT_ConApp_Tp t){
	switch(t){
		case TPAT_ConApp_Tp::APP_TO_NODE: return "Nodes"; break;
		case TPAT_ConApp_Tp::APP_TO_ARC: return "Whole Arcset"; break;
		case TPAT_ConApp_Tp::APP_TO_SEG: return "Segments"; break;
		default: return "Undefined!"; break;
	}
}//====================================================

/**
 *	@brief Print this constraint and its data to the standard output.
 */
void TPAT_Constraint::print() const {
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
void TPAT_Constraint::setAppType(){
	switch(type){
		case TPAT_Constraint_Tp::PSEUDOARC:
		case TPAT_Constraint_Tp::TOF:
		case TPAT_Constraint_Tp::MAX_DELTA_V:
		case TPAT_Constraint_Tp::DELTA_V:
		case TPAT_Constraint_Tp::SEG_CONT_PV:
		case TPAT_Constraint_Tp::SEG_CONT_EX:
			appType = TPAT_ConApp_Tp::APP_TO_ARC;
			break;
		case TPAT_Constraint_Tp::CONT_PV:
		case TPAT_Constraint_Tp::CONT_EX:
			appType = TPAT_ConApp_Tp::APP_TO_SEG;
			break;
		default:
			appType = TPAT_ConApp_Tp::APP_TO_NODE;
			break;
	}
}//====================================================