/**
 *  @file tpat_node.cpp
 *
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

#include "tpat_exceptions.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Default constructor
 */
tpat_node::tpat_node() {
	initArrays();
}

/**
 *	@brief Create a nodeset with a specified state and TOF
 *	@param state an array of 6 position and velocity states (non-dim)
 *	@param tof the non-dimensional time-of-flight
 */
tpat_node::tpat_node(double *state, double tof){
	initArrays();
	std::copy(state, state+6, posVelState);
	this->tof = tof;
}//====================================================

tpat_node::tpat_node(std::vector<double> state, double tof){
	if(state.size() < 6)
		throw tpat_exception("tpat_node: Cannot construct node from state with fewer than 6 elements");

	initArrays();
	std::copy(state.begin(), state.begin()+6, posVelState);
	this->tof = tof;
}//====================================================

tpat_node::tpat_node(const tpat_node &n){
	initArrays();
	copyMe(n);
}//====================================================

tpat_node::~tpat_node(){
	extraParam.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

tpat_node& tpat_node::operator =(const tpat_node &n){
	copyMe(n);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

std::vector<double> tpat_node::getPosVelState() const{
	std::vector<double> temp;
	temp.insert(temp.end(), posVelState, posVelState+6);
	return temp;
}//====================================================

double tpat_node::getExtraParam(int i) const { return extraParam.at(i); }

double tpat_node::getTOF() const { return tof; }

std::vector<bool> tpat_node::getVelCon() const {
	std::vector<bool> temp;
	temp.insert(temp.end(), velCon, velCon+3);
	return temp;
}//====================================================

void tpat_node::setPosVelState(double *array){
	std::copy(array, array+6, posVelState);
}//====================================================

void tpat_node::setPosVelState(std::vector<double> vec){
	if(vec.size() < 6)
		throw tpat_exception("tpat_node::setPosVelState: Cannot set state; input vector has too few elements");
	std::copy(vec.begin(), vec.begin()+6, posVelState);
}//====================================================

void tpat_node::setTOF(double t){ tof = t; }

void tpat_node::setVel_AllCon(){
	velCon[0] = true;
	velCon[1] = true;
	velCon[2] = true;
}//====================================================

void tpat_node::setVel_AllDiscon(){
	velCon[0] = false;
	velCon[1] = false;
	velCon[2] = false;
}//====================================================

void tpat_node::setVelCon(bool data[3]){
	std::copy(data, data+3, velCon);
}//====================================================

void tpat_node::setVelCon(bool xCon, bool yCon, bool zCon){
	velCon[0] = xCon;
	velCon[1] = yCon;
	velCon[2] = zCon;
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_node::copyMe(const tpat_node &n){
	std::copy(n.posVelState, n.posVelState+6, posVelState);
	std::copy(n.velCon, n.velCon+6, velCon);
	tof = n.tof;
	extraParam = n.extraParam;
}//====================================================

void tpat_node::initArrays(){
	for(int i = 0; i < 6; i++)
		posVelState[i] = NAN;

	for(int i = 0; i < 3; i++)
		velCon[i] = true;
}//====================================================
