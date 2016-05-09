/**
 *  @file tpat_traj_cr3bp_ltvp.cpp
 *	@brief 
 *
 *	@author Andrew Cox
 *	@version 
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

#include "tpat_traj_cr3bp_ltvp.hpp"


#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_utilities.hpp"
 
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_sys_data_cr3bp_ltvp* sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_traj_cr3bp_ltvp &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_arcset &a) : tpat_traj(a){
	initExtraParam();
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  CR3BP LTVP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created trajectory
 */
tpat_traj_cr3bp_ltvp* tpat_traj_cr3bp_ltvp::create( const tpat_sys_data *sys) const{
	const tpat_sys_data_cr3bp_ltvp *crSys = static_cast<const tpat_sys_data_cr3bp_ltvp*>(sys);
	return new tpat_traj_cr3bp_ltvp(crSys);
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned trajectory
 */
tpat_traj_cr3bp_ltvp* tpat_traj_cr3bp_ltvp::clone() const{
	return new tpat_traj_cr3bp_ltvp(*this);
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
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
double tpat_traj_cr3bp_ltvp::getJacobiByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_traj_cr3bp_ltvp::getJacobiByIx: invalid index");

	return nodes[ix].getExtraParam(0);
}//====================================================

/**
 *	@brief Retrieve the mass at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@return mass at the specified step (non-dim)
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
double tpat_traj_cr3bp_ltvp::getMassByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_traj_cr3bp_ltvp::getMassByIx: invalid index");

	return nodes[ix].getExtraParam(1);
}//====================================================

/**
 *	@brief Set Jacobi at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@param val value of Jacobi
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
void tpat_traj_cr3bp_ltvp::setJacobiByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_traj_cr3bp_ltvp::setJacobiByIx: invalid index");

	nodes[ix].setExtraParam(0, val);
}//====================================================

/**
 *	@brief Set mass at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@param val mass value (non-dim)
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
void tpat_traj_cr3bp_ltvp::setMassByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_traj_cr3bp_ltvp::setMassByIx: invalid index");

	nodes[ix].setExtraParam(1, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void tpat_traj_cr3bp_ltvp::initExtraParam(){
	// Add another variable for Jacobi Constant, and one for mass
	numExtraParam = 2;
	extraParamRowSize.push_back(1);	// add var for Jacobi
	extraParamRowSize.push_back(1); // add var for Mass
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_cr3bp_ltvp::saveToMat(const char* filename) const{
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
		saveEpoch(matfp, "Time");
		saveSTMs(matfp);
		saveExtraParam(matfp, 0, "Jacobi");
		saveExtraParam(matfp, 1, "Mass");
		sysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================