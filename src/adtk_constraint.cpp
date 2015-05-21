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

adtk_constraint::adtk_constraint(int n) : nodeSize(n){
	data.assign(n, NAN);
}

adtk_constraint::adtk_constraint(int n, constraint_t t) : nodeSize(n){
	type = t;
	data.assign(n, NAN);
}

adtk_constraint::adtk_constraint(int n, constraint_t t, int i, std::vector<double> d) : nodeSize(n){
	type = t;
	node = i;
	data = d;
}

adtk_constraint::adtk_constraint(int n, constraint_t t, int i, double* d) : nodeSize(n){
	type = t;
	node = i;

	for(int j = 0; j < n; j++){
		data.push_back(d[j]);
	}
}

adtk_constraint::adtk_constraint(const adtk_constraint& c) : nodeSize(c.nodeSize){
	type = c.type;
	node = c.node;
	data = c.data;
}//============================================

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

adtk_constraint& adtk_constraint::operator =(const adtk_constraint& c){
	type = c.type;
	node = c.node;
	data = c.data;
	return *this;
}//============================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

adtk_constraint::constraint_t adtk_constraint::getType() const { return type; }
int adtk_constraint::getNode() const { return node; }
std::vector<double> adtk_constraint::getData() const { return data; }

void adtk_constraint::setType(adtk_constraint::constraint_t t){ type = t; }
void adtk_constraint::setNode(int n){ node = n; }
void adtk_constraint::setData(std::vector<double> d){ data = d; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

const char* adtk_constraint::getTypeStr(constraint_t t){
	switch(t){
		case NONE: { return "NONE"; break; }
		case STATE: { return "STATE"; break; }
		case MATCH_ALL: { return "MATCH_ALL"; break; }
		case MATCH_CUST: { return "MATCH_CUST"; break; }
		case MAX_DELTA_V: { return "MAX_DELTA_V"; break; }
		case DELTA_V: { return "DELTA_V"; break; }
		case SP: { return "SP"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//========================================

void adtk_constraint::print(){
	printf("Constraint:\n  Type: %s\n  Node: %d\n  Data: ", getTypeStr(type), node);
	for(int n = 0; n < nodeSize; n++){
		printf("%12.5f ", data[n]);
	}
	printf("\n");
}