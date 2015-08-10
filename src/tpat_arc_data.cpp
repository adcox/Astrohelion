/**
 *  @file tpat_arc_data.cpp
 *	@brief Data object that stores information about an integrated arc
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

#include "tpat_arc_data.hpp"

#include "tpat_sys_data.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

tpat_arc_data::tpat_arc_data(tpat_sys_data *sys){
	sysData = sys;
}//====================================================

tpat_arc_data::tpat_arc_data(const tpat_arc_data &d){
	copyMe(d);
}//====================================================

tpat_arc_data::~tpat_arc_data(){
	steps.clear();
	extraParamRowSize.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

tpat_arc_data& tpat_arc_data:operator =(const tpat_arc_data &d){
	copyMe(d);
	return *this;
}//====================================================

tpat_arc_data& operator +(const tpat_arc_data &lhs, const tpat_arc_data &rhs){
	if(*(lhs.sysData) != *(rhs.sysData))
		throw tpat_exception("tpat_arc_data::+: Cannot concatenate data sets from different systems");

	if(lhs.steps.size() == 0)
		return rhs;

	if(rhs.steps.size() == 0)
		return lhs;

	int skipShift = 0;
	double tol = lhs.tol > rhs.tol ? lhs.tol : rhs.tol;
	if(tol == 0)
		tol = 1e-9;

	if(lhs.steps[0] == rhs.steps[0])
		skipShift = 1;

	tpat_arc_data temp = lhs; // Copy all data from lhs

	// Copy data from rhs, optionally skipping the first step if it matches the final step from lhs
	temp.steps.insert(temp.steps.end(), rhs.steps.begin()+skipShift, rhs.steps.end());

	return temp;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------
/**
 *	@brief Retrieve an acceleration on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last acceleration)
 *	@return the acceleration associated with the specified index
 */
std::vector<double> tpat_arc_data::getAccel(int ix) const{
	if(ix < 0)
		ix += steps.size();
	return steps[ix].getAccel();
}//====================================================

/**
 *	@brief Get a vector of one coordinate for all steps
 *	@param ix the index of the coordinate: 0 = x, 1 = y, etc.
 *	@return a vector containing the specified coordinate for all
 *	integration steps
 */
std::vector<double> tpat_arc_data::getCoord(int ix) const{
	if(i >= 6)
		throw tpat_exception("tpat_arc_data::getCoord: Index Out of Range");

	std::vector<double> coord;
	for(size_t n = 0; n < steps.size(); n++){
		coord.push_back(steps[n].getPosVelState()[ix]);
	}

	return coord;
}//====================================================

/**
 *	@brief Retrieve the length of this arc, in number of steps
 *	@return the length of this arc, in number of steps
 */
int tpat_arc_data::getLength() const {return (int)(steps.size()); }

/**
 *	@brief Retrieve a set of extra parameters for the specified step
 *	@param step the step index. If it is negative, this index will count 
 *	backwards from the end of the arc (e.g. step = -1 will return something from
 *	the final integration step)
 *	@param ix the index of the extra parameter
 *	@return a vector containing the extra parameter at the specified step and index
 */
std::vector<double> tpat_arc_data::getExtraParam(int step, int ix) const{
	if(step < 0)
		step += steps.size();

	if(ix < 0 || ix >= extraParamRowSize.size())
		throw tpat_exception("tpat_arc_data::getExtraParam: parameter index out of bounds");

	int startIx = 0;
	for(int i = 0; i < ix; i++)
		startIx += exrraParamRowSize[i];

	int size = extraParamRowSize[ix];
	std::vector<double> extraParam = steps[step].getExtraParams();
	return std::vector<double>(extraParam.begin()+startIx, extraParam.begin()+startIx + size);
}//====================================================

/**
 *	@brief Retrieve a position-velocity state on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last state)
 *	@return the state associated with the specified index
 */
std::vector<double> tpat_arc_data::getState(int ix) const{
	if(ix < 0)
		ix += steps.size();
	return steps[ix].getPosVelState();
}//====================================================

/**
 *	@brief Retrieve data for one integration step
 *	@param ix the step index; if the index is negative, it will count backwards
 *	from the end of the steps (e.g. ix = -1 will return the last step)
 */
tpat_arc_step tpat_arc_data::getStep(int ix) const{
	if(ix < 0)
		ix += steps.size();
	return steps[ix];
}//====================================================

/**
 *	@brief Retrieve an STM on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last STM)
 *	@return the STM associated with the specified index
 */
tpat_matrix tpat_arc_data::getSTM(int ix) const{
	if(ix < 0)
		ix += steps.size();
	return steps[ix].getSTM();
}//====================================================

/**
 *	@brief Retrieve the tolerance with which data in this object was computed
 *	@return the tolerance with which data in this object was computed
 */
double tpat_arc_data::getTol() const { return tol; }

/**
 *	@brief Set the computational tolerance for this data object
 *	@param d the tolerance
 */
void tpat_arc_data::setTol(double d){ tol = d; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy all data from the input arc data to this one
 *	@param d an arc data object reference
 */
void tpat_arc_data::copyMe(const tpat_arc_data &d){
	steps = d.steps;
	sysData = d.sysData; // Copying ADDRESS of sys_data object
	numExtraParam = d.numExtraParam;
	extraParamRowSize = d.extraParamRowSize;
	tol = d.tol;
}//====================================================





