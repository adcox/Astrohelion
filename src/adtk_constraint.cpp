/**
 *	@file adtk_constraint.cpp
 */

#include "adtk_constraint.hpp"

#include <cmath>
#include <cstdio>
using namespace std;

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Construct a constraint
 */
adtk_constraint::adtk_constraint(){
	data.clear();
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type
 *	@param t constraint type
 */
adtk_constraint::adtk_constraint(constraint_t t){
	type = t;
	data.clear();
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type and data values
 *	@param t constraint type
 *	@param i the node this constriant applies to
 *	@param d data (length n)
 */
adtk_constraint::adtk_constraint(constraint_t t, int i, std::vector<double> d){
	type = t;
	node = i;
	data = d;
}//============================================

/**
 *	@brief Construct a constraint with specified constraint type, and data values
 *	@param t constraint type
 *	@param i the node this constriant applies to
 *	@param d data
 *	@param d_len the number of elements in d_len
 */
adtk_constraint::adtk_constraint(constraint_t t, int i, double* d, int d_len){
	type = t;
	node = i;
	data.reserve(d_len);

	for(int j = 0; j < d_len; j++){
		data.push_back(d[j]);
	}
}//============================================

/**
 *	@brief Create a copy of the specified constraint
 *	@param c a constraint
 */
adtk_constraint::adtk_constraint(const adtk_constraint& c){
	type = c.type;
	node = c.node;
	data = c.data;
}//============================================

/**
 *	@brief Destructor
 */
adtk_constraint::~adtk_constraint(){
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
adtk_constraint& adtk_constraint::operator =(const adtk_constraint& c){
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
adtk_constraint::constraint_t adtk_constraint::getType() const { return type; }

/**
 *	@return the node number associated with this constraint
 */
int adtk_constraint::getNode() const { return node; }

/**
 *	@return the data vector for this constraint
 */
std::vector<double> adtk_constraint::getData() const { return data; }

/**
 *	@return a count of the constrained states; certain constraint types, 
 *	like MATCH_CUST, give the option of constraining a subset of the entire
 *	node.
 */
int adtk_constraint::countConstrainedStates() const{
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
void adtk_constraint::setType(adtk_constraint::constraint_t t){ type = t; }

/**
 *	@brief Set the node index this constraint applies to
 *	@param n node index
 */
void adtk_constraint::setNode(int n){ node = n; }

/**
 *	Set the data for this node (should have nodeSize # elements)
 *	@param d the data, dimensions that match node dimensions
 */
void adtk_constraint::setData(std::vector<double> d){ data = d; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@param t a constraint type
 *	@return a human-readable string representing a constraint type
 */
const char* adtk_constraint::getTypeStr(constraint_t t) const{
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
		case SP: { return "SP"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//========================================

/**
 *	@brief Print this constraint and its data to the standard output.
 */
void adtk_constraint::print() const {
	printf("Constraint:\n  Type: %s\n  Node: %d\n  Data: ", getTypeStr(type), node);
	for(int n = 0; n < ((int)data.size()); n++){
		printf("%12.5f ", data[n]);
	}
	printf("\n");
}