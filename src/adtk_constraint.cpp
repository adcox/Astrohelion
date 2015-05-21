/**
 *	@file adtk_constraint.cpp
 */

#include "adtk_constraint.hpp"

#include <cmath>

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

adtk_constraint::adtk_constraint(const adtk_constraint& c) : nodeSize(c.nodeSize){
	type = c.type;
	node = c.node;
	data = c.data;
}//============================================

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

adtk_constraint& operator =(const adtk_constraint&){
	type = c.type;
	node = c.node;
	data = c.data;
	return *this;
}//============================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

constraint_t adtk_constraint::getType() const { return type; }
int adtk_constraint::getNode() const { return node; }
std::vector<double> adtk_constraint::getData() const { return data; }

void adtk_constraint::setType(constraint_t t){ type = t; }
void adtk_constraint::setNode(int n){ node = n; }
void adtk_constraint::setData(std::vector<double> d){ data = d; }