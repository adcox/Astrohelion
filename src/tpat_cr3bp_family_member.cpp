/**
 *	@file tpat_cr3bp_family_member
 *	@brief Data object for CR3BP family members
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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
 
#include "tpat_cr3bp_family_member.hpp"

#include "tpat_cr3bp_traj.hpp"
#include "tpat_cr3bp_nodeset.hpp"

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

tpat_cr3bp_family_member::tpat_cr3bp_family_member(double *ic, double tof,
	double jc, double xWid, double yWid, double zWid){
	IC.clear();
	IC.insert(IC.begin(), ic, ic+6);
	TOF = tof;
	JC = jc;
	xWidth = xWid;
	yWidth = yWid;
	zWidth = zWid;
}//====================================================

/**
 *	@brief Create a family member from a trajectory object
 */
tpat_cr3bp_family_member::tpat_cr3bp_family_member(const tpat_cr3bp_traj traj){
	IC = traj.getState(0);
	TOF = traj.getTime(-1);
	JC = traj.getJC(0);

	std::vector<double> x = traj.getCoord(0);
	std::vector<double> y = traj.getCoord(0);
	std::vector<double> z = traj.getCoord(0);

	xWidth = *std::max_element(x.begin(), x.end()) - *std::min_element(x.begin(), x.end());
	yWidth = *std::max_element(y.begin(), y.end()) - *std::min_element(y.begin(), y.end());
	zWidth = *std::max_element(z.begin(), z.end()) - *std::min_element(z.begin(), z.end());
}//===================================================

/**
 *	@brief Copy constructor
 */
tpat_cr3bp_family_member::tpat_cr3bp_family_member(const tpat_cr3bp_family_member& mem){
	copyMe(mem);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_cr3bp_family_member::~tpat_cr3bp_family_member(){
	IC.clear();
}//===================================================


//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 */
tpat_cr3bp_family_member& tpat_cr3bp_family_member::operator= (const tpat_cr3bp_family_member& mem){
	copyMe(mem);
	return *this;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the initial state for this trajectory (non-dim)
 */
std::vector<double> tpat_cr3bp_family_member::getIC() const { return IC; }

/**
 *	@breif Retrieve the Time-Of-Flight along this trajectory (non-dim)
 */
double tpat_cr3bp_family_member::getTOF() const { return TOF; }

/**
 *	@brief Retrieve the Jacobi Constant for this trajectory
 */
double tpat_cr3bp_family_member::getJC() const { return JC; }

/**
 *	@brief Retrieve the maximum width in the x-direction
 */
double tpat_cr3bp_family_member::getXWidth() const { return xWidth; }

/**
 *	@brief Retrieve the maximum width in the y-direction
 */
double tpat_cr3bp_family_member::getYWidth() const { return yWidth; }

/**
 *	@brief Retrieve the maximum width in the z-direction
 */
double tpat_cr3bp_family_member::getZWidth() const { return zWidth; }

/**
 *	@brief Set the initial state
 *	@param ic The initial state (non-dim)
 */
void tpat_cr3bp_family_member::setIC( std::vector<double> ic ){ IC = ic; }

/**
 *	@brief Set the time-of-flight
 *	@param tof The time-of-flight (non-dim)
 */
void tpat_cr3bp_family_member::setTOF( double tof ){ TOF = tof; }

/**
 *	@brief Set the Jacobi Constant
 *	@param jc The Jacobi Constant (non-dim)
 */
void tpat_cr3bp_family_member::setJC( double jc ){ JC = jc; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void tpat_cr3bp_family_member::setXWidth(double w){ xWidth = w; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void tpat_cr3bp_family_member::setYWidth(double w){ yWidth = w; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void tpat_cr3bp_family_member::setZWidth(double w){ zWidth = w; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy an input family member into this one
 *	@param mem some other family member
 */
void tpat_cr3bp_family_member::copyMe(const tpat_cr3bp_family_member& mem){
	IC = mem.IC;
	JC = mem.JC;
	TOF = mem.TOF;
	xWidth = mem.xWidth;
	yWidth = mem.yWidth;
	zWidth = mem.zWidth;
}//===================================================