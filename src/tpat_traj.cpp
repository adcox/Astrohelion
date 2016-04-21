/**
 *  @file tpat_traj.cpp
 *	@brief Stores information about a trajectory
 *
 *	@author Andrew Cox
 *	@version August 30, 2015
 *	@copyright GNU GPL v3.0
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

#include "tpat_traj.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset.hpp"
#include "tpat_traj_step.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param data a pointer to a system data object
 */
tpat_traj::tpat_traj(const tpat_sys_data *data) : tpat_arc_data(data) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
tpat_traj::tpat_traj(const tpat_traj &t) : tpat_arc_data(t) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
tpat_traj::tpat_traj(const tpat_arc_data &a) : tpat_arc_data(a) {
	initExtraParam();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------


//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the time along the trajectory at a specific step
 *	@param ix step index; if < 0, it will count backwards from end of trajectory
 *	@return the non-dimensional time along the trajectory at the specified step
 */
double tpat_traj::getTime(int ix) const {
	if(ix < 0)
		ix += steps.size();
	
	if(ix < 0 || ix > ((int)steps.size()))
		throw tpat_exception("tpat_traj::getTime: invalid index");

	tpat_traj_step step(steps[ix]);
	return step.getTime();
}//====================================================

/**
 *  @brief Get the time of flight on this trajectory
 *  @return time of flight on this trajectory, non-dimensional time
 */
double tpat_traj::getTOF() const {
	try{
		return getTime(-1) - getTime(0);
	}catch(tpat_exception &e){
		printErr("Could not compute TOF; index error:\n");
		throw(e);
	}
}//====================================================

/**
 *	@brief Retrieve the specified step
 *	@param ix step index; if < 0, it will count backwards from end of trajectory
 *	@return the requested trajectory step object
 */
tpat_traj_step tpat_traj::getStep(int ix) const{
	if(ix < 0)
		ix += steps.size();

	if(ix < 0 || ix > ((int)steps.size()))
		throw tpat_exception("tpat_traj::getStep: invalid index");

	return tpat_traj_step(steps[ix]);
}//====================================================

/**
 *	@brief Append a step to the end of the trajectory
 *	@param s a new step
 */
void tpat_traj::appendStep(tpat_traj_step s){ steps.push_back(s); }

/**
 *	@brief Set the time for the specified node
 *	@param ix node index; if < 0, counts backwards from end
 *	@param val non-dimensional time value for the specified node
 */
void tpat_traj::setTime(int ix, double val){
	if(ix < 0)
		ix += steps.size();

	if(ix < 0 || ix > ((int)steps.size()))
		throw tpat_exception("tpat_traj::setTime: invalid index");

	tpat_traj_step *step = static_cast<tpat_traj_step*>(&(steps[ix]));
	step->setTime(val);
}//====================================================

/**
 *  @brief Shift all time values by a constant amount
 *  @details This can be useful for use with the EM2SE and SE2EM functions
 * 
 *  @param amount a constant, non-dimensional time shift to apply to 
 *  all time values for points on this trajectory
 */
void tpat_traj::shiftAllTimes(double amount){
	for(size_t i = 0; i < steps.size(); i++){
		tpat_traj_step *step = static_cast<tpat_traj_step*>(&(steps[i]));
		step->setTime(step->getTime() + amount);
	}
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Discretize a trajectory into a set of nodes without using integration
 *
 *	This method uses the existing steps in a trajectory to create nodes. Arc segments
 *	are created with equal number of steps regardless of the time or arc length separating
 *	the nodes.
 *
 *	@param numNodes number of nodes, including both the initial and final states on 
 *	the trajectory, which are always included
 *	@return a nodeset with the specified number of nodes
 */
tpat_nodeset tpat_traj::discretize(int numNodes) const{
	if(numNodes < 2)
		throw tpat_exception("tpat_traj::discretize: Cannot split a trajectory into fewer than 2 nodes");

	if(numNodes > (int)(steps.size())){
		printWarn("tpat_traj::discretize: User requested more nodes than there are states; returning one node per step, will not meet requested number of nodes\n");
		numNodes = steps.size();
	}

	double stepSize = (double)(steps.size()-1)/((double)numNodes - 1.0);

	tpat_nodeset nodes(sysData);
	int n = 0;
	while(n < numNodes){
		// Round the step number
		int ix = std::floor(n*stepSize);
		int nextIx = std::floor((n+1)*stepSize);

		if(n == numNodes-2)
			printf("");

		// Create a node from this step
		std::vector<double> state = getState(ix);
		double tof = n < numNodes-1 ? getTime(nextIx) - getTime(ix) : NAN;
		tpat_node node(state, tof);
		nodes.appendNode(node);

		n++;
	}

	return nodes;
}//=================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj::saveToMat(const char* filename) const{
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("tpat_traj::saveToMat: Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
		sysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Print a useful message describing this trajectory to the standard output
 */
void tpat_traj::print() const {
	printf("This is a trajectory\n\tTODO - MAKE THIS MESSAGE USEFUL\n");
}//====================================================

/**
 *	@brief Save the times at each integration step
 *	@param file a pointer to an open Mat file
 */
void tpat_traj::saveTime(mat_t *file) const{
	saveExtraParam(file, 0, "Time");
}//====================================================

/**
 *	@brief Initialize the extra param vector for trajectory-specific info
 */
void tpat_traj::initExtraParam(){
	// ExtraParam = [time]
	numExtraParam = 1;
	extraParamRowSize.push_back(1);
}//====================================================

/**
 *  @brief Populate data in this trajectory from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 */
void tpat_traj::readFromMat(const char *filepath){

	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw tpat_exception("tpat_traj: Could not load data from file");
	}
	initStepVectorFromMat(matfp, "State");
	readStateFromMat(matfp, "State");
	readAccelFromMat(matfp);
	readSTMFromMat(matfp);
	readExtraParamFromMat(matfp, 0, "Time");

	Mat_Close(matfp);
}//====================================================