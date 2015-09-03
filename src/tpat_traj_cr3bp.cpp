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
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
 
tpat_traj_cr3bp::tpat_traj_cr3bp(tpat_sys_data_cr3bp *sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

tpat_traj_cr3bp::tpat_traj_cr3bp(const tpat_traj_cr3bp &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

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

double tpat_traj_cr3bp::getJacobi(int ix) const{
	if(ix < 0)
		ix += steps.size();
	tpat_arc_step step = steps[ix];
	return step.getExtraParam(2);
}//====================================================

void tpat_traj_cr3bp::setJacobi(int ix, double val){
	if(ix < 0)
		ix += steps.size();

	steps[ix].setExtraParam(ix, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_traj_cr3bp::initExtraParam(){
	// This function in tpat_traj was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the time variable

	// Add another variable for Jacobi Constant
	numExtraParam = 2;
	extraParamRowSize.push_back(1);
}//====================================================
