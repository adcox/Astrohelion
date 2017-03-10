/**
 *	\file Constraint.cpp
 *	\brief Data object that stores information about a node constraint
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
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

#include "Constraint.hpp"

#include <cmath>
#include <cstdio>

namespace astrohelion{
//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	\brief Construct a constraint
 */
Constraint::Constraint(){
	data.clear();
}//====================================================

/**
 *	\brief Construct a constraint with specified constraint type
 *	\param type constraint type
 */
Constraint::Constraint(Constraint_tp type){
	this->type = type;
	data.clear();
	setAppType();
}//====================================================

/**
 *	\brief Construct a constraint with specified constraint type and data values
 *	\param type constraint type
 *	\param id the ID of the object this constriant applies to
 *	\param data data vector (length n)
 */
Constraint::Constraint(Constraint_tp type, int id, std::vector<double> data){
	this->type = type;
	this->id = id;
	this->data = data;
	setAppType();
}//====================================================

/**
 *	\brief Construct a constraint with specified constraint type, and data values
 *	\param type constraint type
 *	\param id the ID of the object this constriant applies to
 *	\param data data vector
 *	\param data_len the number of elements in d_len
 */
Constraint::Constraint(Constraint_tp type, int id, const double* data, int data_len){
	this->type = type;
	this->id = id;
	this->data.insert(this->data.begin(), data, data + data_len);
	setAppType();
}//====================================================

/**
 *	\brief Create a copy of the specified constraint
 *	\param c a constraint
 */
Constraint::Constraint(const Constraint& c){
	copyMe(c);
}//====================================================

/**
 *	\brief Destructor
 */
Constraint::~Constraint(){}

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	\brief Assignment operator
 *	\param c a constraint
 *	@return this constraint, now equal to c
 */
Constraint& Constraint::operator =(const Constraint& c){
	copyMe(c);
	return *this;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Retrieve the application type for this constraint, i.e., what type of
 *  objects it controls and can be applied to
 *  @return the application type
 */
ConstraintApp_tp Constraint::getAppType() const{ return appType; }

/**
 *	@return what type of constraint this is
 */
Constraint_tp Constraint::getType() const { return type; }

/**
 *	@return the ID that identifies the object constrained by this constraint
 */
int Constraint::getID() const { return id; }

/**
 *	@return the data vector for this constraint
 */
std::vector<double> Constraint::getData() const { return data; }

/**
 *  \brief Retrieve the first data value that is not an NAN
 *  @return The value of the first data value. If no data
 *  value is located, NAN is returned.
 *  @see getFirstDataValue(int*)
 */
double Constraint::getFirstDataValue() const{
	int ix = 0;
	return getFirstDataValue(&ix);
}//====================================================

/**
 *  \brief Retrieve the first data value that is not an NAN
 * 
 *  \param ix The index of the first data value will be
 *  stored in this integer. If no data value is located,
 *  ix will be set to -1.
 *  @return The value of the first data value. If no data
 *  value is located, NAN is returned.
 */
double Constraint::getFirstDataValue(int *ix) const{
	for(unsigned int i = 0; i < data.size(); i++){
		if(!std::isnan(data[i])){
			*ix = static_cast<int>(i);
			return data[i];
		}
	}

	*ix = -1;
	return NAN;
}//====================================================

/**
 *	@return a count of the constrained states; certain constraint types, 
 *	like Constraint_tp::MATCH_CUST, give the option of constraining a subset of the entire
 *	node.
 */
int Constraint::countConstrainedStates() const{
	int count = 0;
	for(unsigned int n = 0; n < data.size(); n++){
		count += !std::isnan(data.at(n));
	}
	return count;
}//====================================================

/**
 *	\brief Set the constraint type
 *	\param t the type
 */
void Constraint::setType(Constraint_tp t){
	type = t;
	setAppType(); 	// Update, if necessary
}//====================================================

/**
 *	\brief Set the object ID this constraint applies to
 *	\param n the ID
 */
void Constraint::setID(int n){ id = n; }

/**
 *	Set the data for this id (should have nodeSize # elements)
 *	\param d the data, dimensions that match node dimensions
 */
void Constraint::setData(std::vector<double> d){ data = d; }

/**
 *  \brief Set the data for this constraint
 * 
 *  \param dat an array of data values
 *  \param len number of elements in <tt>dat</tt>
 */
void Constraint::setData(const double *dat, int len){
	data.clear();
	data.insert(data.begin(), dat, dat+len);
}//====================================================

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *  \brief Copy the constraint
 * 
 *  \param c reference to a constraint object
 */
void Constraint::copyMe(const Constraint &c){
	appType = c.appType;
	type = c.type;
	id = c.id;
	data = c.data;
}//====================================================

/**
 *  \brief Get a human-redable string representing the constraint type
 *  @return a human-redable string representing the constraint type
 */
const char* Constraint::getTypeStr() const{ return getConTypeStr(type); }

/**
 *	\param t a constraint type
 *	@return a human-readable string representing a constraint type
 */
const char* Constraint::getConTypeStr(Constraint_tp t){
	switch(t){
		case Constraint_tp::NONE: { return "NONE"; break; }
		case Constraint_tp::STATE: { return "Constraint_tp::STATE"; break; }
		case Constraint_tp::MATCH_ALL: { return "Constraint_tp::MATCH_ALL"; break; }
		case Constraint_tp::MATCH_CUST: { return "Constraint_tp::MATCH_CUST"; break; }
		case Constraint_tp::DIST: { return "Constraint_tp::DIST"; break; }
		case Constraint_tp::MIN_DIST: { return "Constraint_tp::MIN_DIST"; break; }
		case Constraint_tp::MAX_DIST: { return "Constraint_tp::MAX_DIST"; break; }
		case Constraint_tp::MAX_DELTA_V: { return "Constraint_tp::MAX_DELTA_V"; break; }
		case Constraint_tp::DELTA_V: { return "Constraint_tp::DELTA_V"; break; }
		case Constraint_tp::JC: { return "JC"; break; }
		case Constraint_tp::SP: { return "SP"; break; }
		case Constraint_tp::SP_RANGE: { return "SP_RANGE"; break; }
		case Constraint_tp::SP_DIST: { return "SP_DIST"; break; }
		case Constraint_tp::SP_MAX_DIST: { return "SP_MAX_DIST"; break; }
		case Constraint_tp::TOF: { return "TOF"; break; }
		case Constraint_tp::APSE: {return "APSE"; break; }
		case Constraint_tp::CONT_PV: {return "CONTINUOUS Constraint_tp::SEG_2NODE POSITION_VELOCITY"; break; }
		case Constraint_tp::CONT_EX: { return "CONTINUOUS Constraint_tp::SEG_2NODE EXTRA"; break; }
		case Constraint_tp::SEG_CONT_PV: { return "CONTINUOUS Constraint_tp::SEG_2Constraint_tp::SEG_ POSITION_VELOCITY"; break; }
		case Constraint_tp::SEG_CONT_EX: { return "CONTINUOUS Constraint_tp::SEG_2Constraint_tp::SEG_ EXTRA";  break; }
		case Constraint_tp::PSEUDOARC: { return "PSEUDO-ARCLENGTH"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//====================================================

/**
 *  \param t a constraint application type
 *  @return a human-readable string representing an application type
 */
const char* Constraint::getAppTypeStr(ConstraintApp_tp t){
	switch(t){
		case ConstraintApp_tp::APP_TO_NODE: return "Node"; break;
		case ConstraintApp_tp::APP_TO_ARC: return "Whole Arcset"; break;
		case ConstraintApp_tp::APP_TO_SEG: return "Segment"; break;
		default: return "Undefined!"; break;
	}
}//====================================================

/**
 *	\brief Print this constraint and its data to the standard output.
 */
void Constraint::print() const {
	printf("Constraint:\n  Type: %s\n  Applies to: %s (ID %d)\n  Data: ",
		getConTypeStr(type), getAppTypeStr(appType), id);
	for(unsigned int n = 0; n < data.size(); n++){
		printf("%12.5f ", data[n]);
	}
	printf("\n");
}//====================================================

/**
 *  \brief Set the application type appropriately based on the constraint type
 */
void Constraint::setAppType(){
	switch(type){
		case Constraint_tp::PSEUDOARC:
		case Constraint_tp::TOF:
		case Constraint_tp::MAX_DELTA_V:
		case Constraint_tp::DELTA_V:
		case Constraint_tp::SEG_CONT_PV:
		case Constraint_tp::SEG_CONT_EX:
			appType = ConstraintApp_tp::APP_TO_ARC;
			break;
		case Constraint_tp::CONT_PV:
		case Constraint_tp::CONT_EX:
			appType = ConstraintApp_tp::APP_TO_SEG;
			break;
		default:
			appType = ConstraintApp_tp::APP_TO_NODE;
			break;
	}
}//====================================================
}// END of Astrohelion namespace