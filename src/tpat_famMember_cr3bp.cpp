/**
 *	@file tpat_fammember_cr3bp.cpp
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
 
#include "tpat_famMember_cr3bp.hpp"

#include "tpat_constants.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Construct a family member given all the required parameters
 *
 *	@param ic a 6-element vector describing the initial state on the trajectory (non-dim)
 *	@param tof the time-of-flight along the trajectory (non-dim)
 *	@param jc the Jacobi Constant on this trajectory
 *	@param xAmp the maximum amplitude in the x-direction (non-dim)
 *	@param yAmp the maximum amplitude in the y-direction (non-dim)
 *	@param zAmp the maximum amplitude in the z-direction (non-dim)
 */
TPAT_FamMember_CR3BP::TPAT_FamMember_CR3BP(double *ic, double tof,
	double jc, double xAmp, double yAmp, double zAmp){
	
	IC.clear();
	IC.insert(IC.begin(), ic, ic+6);
	TOF = tof;
	JC = jc;
	xAmplitude = xAmp;
	yAmplitude = yAmp;
	zAmplitude = zAmp;
}//====================================================

/**
 *	@brief Create a family member from a trajectory object
 *	@param traj a trajectory reference
 */
TPAT_FamMember_CR3BP::TPAT_FamMember_CR3BP(const TPAT_Traj_CR3BP traj){
	IC = traj.getStateByIx(0);
	TOF = traj.getTotalTOF();
	JC = traj.getJacobiByIx(0);

	std::vector<double> x = traj.getCoord(0);
	std::vector<double> y = traj.getCoord(1);
	std::vector<double> z = traj.getCoord(2);

	double xMax = *std::max_element(x.begin(), x.end());
	double xMin = *std::min_element(x.begin(), x.end());
	double yMax = *std::max_element(y.begin(), y.end());
	double yMin = *std::min_element(y.begin(), y.end());
	double zMax = *std::max_element(z.begin(), z.end());
	double zMin = *std::min_element(z.begin(), z.end());

	xAmplitude = std::max(xMax, -1*xMin);
	yAmplitude = std::max(yMax, -1*yMin);
	zAmplitude = std::max(zMax, -1*zMin);

	// printf("\nxAmp = %.4f\nyAmp = %.4f\nzAmp = %.4f\n", xAmplitude, yAmplitude, zAmplitude);
}//===================================================

/**
 *	@brief Copy constructor
 *	@param mem a family member reference
 */
TPAT_FamMember_CR3BP::TPAT_FamMember_CR3BP(const TPAT_FamMember_CR3BP& mem){
	copyMe(mem);
}//====================================================

/**
 *	@brief Destructor
 */
TPAT_FamMember_CR3BP::~TPAT_FamMember_CR3BP(){}


//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param mem a family member reference
 */
TPAT_FamMember_CR3BP& TPAT_FamMember_CR3BP::operator= (const TPAT_FamMember_CR3BP& mem){
	copyMe(mem);
	return *this;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve a vector of eigenvalues (of the final STM, likely the Monodromy matrix)
 */
std::vector<cdouble> TPAT_FamMember_CR3BP::getEigVals() const { return eigVals; }

/**
 *	@brief Retrieve the initial state for this trajectory (non-dim)
 */
std::vector<double> TPAT_FamMember_CR3BP::getIC() const { return IC; }

/**
 *	@brief Retrieve the Time-Of-Flight along this trajectory (non-dim)
 */
double TPAT_FamMember_CR3BP::getTOF() const { return TOF; }

/**
 *	@brief Retrieve the Jacobi Constant for this trajectory
 */
double TPAT_FamMember_CR3BP::getJacobi() const { return JC; }

/**
 *	@brief Retrieve the maximum amplitude in the x-direction
 */
double TPAT_FamMember_CR3BP::getXAmplitude() const { return xAmplitude; }

/**
 *	@brief Retrieve the maximum amplitude in the y-direction
 */
double TPAT_FamMember_CR3BP::getYAmplitude() const { return yAmplitude; }

/**
 *	@brief Retrieve the maximum amplitude in the z-direction
 */
double TPAT_FamMember_CR3BP::getZAmplitude() const { return zAmplitude; }

/**
 *	@brief Set the eigenvalues for this orbit
 *
 *	These should be the eigenvalues of the final STM and/or Monodromy matrix
 *	@param vals the eigenvalues
 *	@throws TPAT_Exception if <tt>vals</tt> does not have six elements
 */
void TPAT_FamMember_CR3BP::setEigVals(std::vector<cdouble> vals) {
	if(vals.size() != 6)
		throw TPAT_Exception("TPAT_FamMember_CR3BP::setEigVals: There must be 6 eigenvalues");
	eigVals = vals;
}//====================================================

/**
 *	@brief Set the initial state
 *	@param ic The initial state (non-dim)
 *	@throws TPAT_Exception if <tt>ic</tt> does not have six elements
 */
void TPAT_FamMember_CR3BP::setIC( std::vector<double> ic ){
	if(ic.size() != 6)
		throw TPAT_Exception("TPAT_FamMember_CR3BP::setIC: There must be 6 elements!");
	IC = ic;
}//====================================================

/**
 *	@brief Set the time-of-flight
 *	@param tof The time-of-flight (non-dim)
 */
void TPAT_FamMember_CR3BP::setTOF( double tof ){ TOF = tof; }

/**
 *	@brief Set the Jacobi Constant
 *	@param jc The Jacobi Constant (non-dim)
 */
void TPAT_FamMember_CR3BP::setJacobi( double jc ){ JC = jc; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void TPAT_FamMember_CR3BP::setXAmplitude(double w){ xAmplitude = w; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void TPAT_FamMember_CR3BP::setYAmplitude(double w){ yAmplitude = w; }

/**
 *	@brief Set the width of this trajectory in the x-direction (non-dim)
 */
void TPAT_FamMember_CR3BP::setZAmplitude(double w){ zAmplitude = w; }

/**
 *  @brief Convert the family member object to a trajectory object
 *  @details This function simply numerically integrates the IC for the
 *  specified TOF to produce the trajectory object
 * 
 *  @param sys The system the family member exists in
 *  @return A trajectory object
 */
TPAT_Traj_CR3BP TPAT_FamMember_CR3BP::toTraj(const TPAT_Sys_Data_CR3BP *sys){
	TPAT_Sim_Engine sim;
	TPAT_Traj_CR3BP traj(sys);
	sim.runSim(IC, TOF, &traj);
	return traj;
}//====================================================

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy an input family member into this one
 *	@param mem some other family member
 */
void TPAT_FamMember_CR3BP::copyMe(const TPAT_FamMember_CR3BP& mem){
	eigVals = mem.eigVals;
	IC = mem.IC;
	JC = mem.JC;
	TOF = mem.TOF;
	xAmplitude = mem.xAmplitude;
	yAmplitude = mem.yAmplitude;
	zAmplitude = mem.zAmplitude;
}//===================================================