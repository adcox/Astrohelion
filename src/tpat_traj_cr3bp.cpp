/**
 *  @file tpat_traj_cr3bp.cpp
 *	@brief Derivative of tpat_traj, specific to CR3BP
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
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

#include "tpat.hpp"

#include "tpat_traj_cr3bp.hpp"

#include "tpat_arc_step.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
 
/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(tpat_sys_data_cr3bp *sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(const tpat_traj_cr3bp &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(const tpat_arc_data &a) : tpat_traj(a){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from a nodeset
 *
 *	This algorithm will concatenate trajectories integrated from each node in 
 *	the nodeset. It does not check to make sure the arcs are continuous; that
 *	is up to you
 *
 *	@param nodes a nodeset
 *	@return a trajectory formed from the integrated nodeset
 */
tpat_traj_cr3bp tpat_traj_cr3bp::fromNodeset(tpat_nodeset_cr3bp nodes){
	tpat_sys_data_cr3bp *sys = static_cast<tpat_sys_data_cr3bp*>(nodes.getSysData());
	tpat_simulation_engine simEngine(sys);
	simEngine.clearEvents();	// don't trigger crashes; assume this has been taken care of already
	tpat_traj_cr3bp totalTraj(sys);

	for(int n = 0; n < nodes.getNumNodes()-1; n++){
		simEngine.runSim(nodes.getNode(n).getPosVelState(), nodes.getTOF(n));

		if(n == 0){
			totalTraj = simEngine.getCR3BP_Traj();
		}else{
			tpat_traj_cr3bp temp = totalTraj + simEngine.getCR3BP_Traj();
			totalTraj = temp;
		}
	}

	return totalTraj;
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the value of Jacobi's Constant at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@return Jacobi at the specified step
 */
double tpat_traj_cr3bp::getJacobi(int ix) const{
	if(ix < 0)
		ix += steps.size();
	tpat_arc_step step = steps[ix];
	return step.getExtraParam(1);
}//====================================================

/**
 *	@brief Set Jacobi at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@param val value of Jacobi
 */
void tpat_traj_cr3bp::setJacobi(int ix, double val){
	if(ix < 0)
		ix += steps.size();

	steps[ix].setExtraParam(1, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void tpat_traj_cr3bp::initExtraParam(){
	// This function in tpat_traj was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the time variable

	// Add another variable for Jacobi Constant
	numExtraParam = 2;
	extraParamRowSize.push_back(1);
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_cr3bp::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file @version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
		saveExtraParam(matfp, 1, "Jacobi");
		sysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================
