/**
 *  @file Traj_cr3bp_ltvp.cpp
 *	@brief 
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "Traj_cr3bp_ltvp.hpp"


#include "SysData_cr3bp_ltvp.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"
 

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
Traj_cr3bp_ltvp::Traj_cr3bp_ltvp(const SysData_cr3bp_ltvp* sys) : Traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
Traj_cr3bp_ltvp::Traj_cr3bp_ltvp(const Traj_cr3bp_ltvp &t) : Traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
Traj_cr3bp_ltvp::Traj_cr3bp_ltvp(const BaseArcset &a) : Traj(a){
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
baseArcsetPtr Traj_cr3bp_ltvp::create( const SysData *sys) const{
	const SysData_cr3bp_ltvp *crSys = static_cast<const SysData_cr3bp_ltvp*>(sys);
	return baseArcsetPtr(new Traj_cr3bp_ltvp(crSys));
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
baseArcsetPtr Traj_cr3bp_ltvp::clone() const{
	return baseArcsetPtr(new Traj_cr3bp_ltvp(*this));
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
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
double Traj_cr3bp_ltvp::getJacobiByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_cr3bp_ltvp::getJacobiByIx: invalid index");

	return nodes[ix].getExtraParam("J");
}//====================================================

/**
 *	@brief Retrieve the mass at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@return mass at the specified step (non-dim)
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
double Traj_cr3bp_ltvp::getMassByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_cr3bp_ltvp::getMassByIx: invalid index");

	return nodes[ix].getExtraParam("m");
}//====================================================

/**
 *	@brief Set Jacobi at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@param val value of Jacobi
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
void Traj_cr3bp_ltvp::setJacobiByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_cr3bp_ltvp::setJacobiByIx: invalid index");

	nodes[ix].setExtraParam("J", val);
}//====================================================

/**
 *	@brief Set mass at the specified step
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@param val mass value (non-dim)
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
void Traj_cr3bp_ltvp::setMassByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_cr3bp_ltvp::setMassByIx: invalid index");

	nodes[ix].setExtraParam("m", val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void Traj_cr3bp_ltvp::initExtraParam(){
	// Add another variable for Jacobi Constant, and one for mass
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void Traj_cr3bp_ltvp::saveToMat(const char* filename) const{
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
		astrohelion::printErr("Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveEpoch(matfp, "Time");
		saveTOF(matfp, "TOFs");
		saveSTMs(matfp);
		saveExtraParam(matfp, "J", "Jacobi");
		saveExtraParam(matfp, "m", "Mass");
		pSysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================



}