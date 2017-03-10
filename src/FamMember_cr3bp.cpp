/**
 *	\file FamMember_cr3bp.cpp
 *	\brief Data object for CR3BP family members
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
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


#include <Eigen/Eigenvalues>
#include <iostream>

#include "FamMember_cr3bp.hpp"

#include "Common.hpp"
#include "Exceptions.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"

namespace astrohelion{
//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	\brief Construct a family member given all the required parameters
 *
 *	\param ic a 6-element vector describing the initial state on the trajectory (non-dim)
 *	\param tof the time-of-flight along the trajectory (non-dim)
 *	\param jc the Jacobi Constant on this trajectory
 *	\param xAmp the maximum amplitude in the x-direction (non-dim)
 *	\param yAmp the maximum amplitude in the y-direction (non-dim)
 *	\param zAmp the maximum amplitude in the z-direction (non-dim)
 */
FamMember_cr3bp::FamMember_cr3bp(double *ic, double tof,
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
 *	\brief Create a family member from a trajectory object
 *	\param traj a trajectory reference
 */
FamMember_cr3bp::FamMember_cr3bp(const Traj_cr3bp traj){
	IC = traj.getStateByIx(0);
	TOF = traj.getTotalTOF();
	JC = traj.getJacobiByIx(0);
	stm = traj.getSTMByIx(-1);

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

	// Compute Eigenvalues and Eigenvectors from the STM	
	Eigen::EigenSolver<MatrixXRd> eigensolver(stm);
    if(eigensolver.info() != Eigen::Success)
        throw Exception("FamMember_cr3bp::FamMember_cr3bp: Could not compute eigenvalues of STM");

    Eigen::VectorXcd vals = eigensolver.eigenvalues();
    // std::vector<cdouble> vals(vals.data(), vals.data()+6);
    eigVals = std::vector<cdouble>(vals.data(), vals.data()+6);
    eigVecs = eigensolver.eigenvectors();

    // Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    // std::cout << eigVecs.format(OctaveFmt) << std::endl;

	// printf("\nxAmp = %.4f\nyAmp = %.4f\nzAmp = %.4f\n", xAmplitude, yAmplitude, zAmplitude);
}//===================================================

/**
 *	\brief Copy constructor
 *	\param mem a family member reference
 */
FamMember_cr3bp::FamMember_cr3bp(const FamMember_cr3bp& mem){
	copyMe(mem);
}//====================================================

/**
 *	\brief Destructor
 */
FamMember_cr3bp::~FamMember_cr3bp(){}


//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	\brief Assignment operator
 *	\param mem a family member reference
 */
FamMember_cr3bp& FamMember_cr3bp::operator= (const FamMember_cr3bp& mem){
	copyMe(mem);
	return *this;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	\brief Retrieve a vector of eigenvalues (of the final STM, likely the Monodromy matrix)
 */
std::vector<cdouble> FamMember_cr3bp::getEigVals() const { return eigVals; }

MatrixXRcd FamMember_cr3bp::getEigVecs() const { return eigVecs; }

/**
 *	\brief Retrieve the initial state for this trajectory (non-dim)
 */
std::vector<double> FamMember_cr3bp::getIC() const { return IC; }

/**
 *	\brief Retrieve the Time-Of-Flight along this trajectory (non-dim)
 */
double FamMember_cr3bp::getTOF() const { return TOF; }

/**
 *  \brief Retrieve the STM for this family member
 */
MatrixXRd FamMember_cr3bp::getSTM() const { return stm; }

/**
 *	\brief Retrieve the Jacobi Constant for this trajectory
 */
double FamMember_cr3bp::getJacobi() const { return JC; }

/**
 *	\brief Retrieve the maximum amplitude in the x-direction
 */
double FamMember_cr3bp::getXAmplitude() const { return xAmplitude; }

/**
 *	\brief Retrieve the maximum amplitude in the y-direction
 */
double FamMember_cr3bp::getYAmplitude() const { return yAmplitude; }

/**
 *	\brief Retrieve the maximum amplitude in the z-direction
 */
double FamMember_cr3bp::getZAmplitude() const { return zAmplitude; }

/**
 *	\brief Set the eigenvalues for this orbit
 *
 *	These should be the eigenvalues of the final STM and/or Monodromy matrix
 *	\param vals the eigenvalues
 *	\throws Exception if <tt>vals</tt> does not have six elements
 */
void FamMember_cr3bp::setEigVals(std::vector<cdouble> vals) {
	if(vals.size() != 6)
		throw Exception("FamMember_cr3bp::setEigVals: There must be 6 eigenvalues");
	eigVals = vals;
}//====================================================

void FamMember_cr3bp::setEigVecs(MatrixXRcd vecs){ eigVecs = vecs; }

/**
 *	\brief Set the initial state
 *	\param ic The initial state (non-dim)
 *	\throws Exception if <tt>ic</tt> does not have six elements
 */
void FamMember_cr3bp::setIC( std::vector<double> ic ){
	if(ic.size() != 6)
		throw Exception("FamMember_cr3bp::setIC: There must be 6 elements!");
	IC = ic;
}//====================================================

/**
 *	\brief Set the time-of-flight
 *	\param tof The time-of-flight (non-dim)
 */
void FamMember_cr3bp::setTOF( double tof ){ TOF = tof; }

/**
 *  \brief Set the STM
 *  \param STM the STM
 */
void FamMember_cr3bp::setSTM( MatrixXRd STM){ stm = STM; }

/**
 *	\brief Set the Jacobi Constant
 *	\param jc The Jacobi Constant (non-dim)
 */
void FamMember_cr3bp::setJacobi( double jc ){ JC = jc; }

/**
 *	\brief Set the width of this trajectory in the x-direction (non-dim)
 */
void FamMember_cr3bp::setXAmplitude(double w){ xAmplitude = w; }

/**
 *	\brief Set the width of this trajectory in the x-direction (non-dim)
 */
void FamMember_cr3bp::setYAmplitude(double w){ yAmplitude = w; }

/**
 *	\brief Set the width of this trajectory in the x-direction (non-dim)
 */
void FamMember_cr3bp::setZAmplitude(double w){ zAmplitude = w; }

/**
 *  \brief Convert the family member object to a trajectory object
 *  @details This function simply numerically integrates the IC for the
 *  specified TOF to produce the trajectory object
 * 
 *  \param sys The system the family member exists in
 *  @return A trajectory object
 */
Traj_cr3bp FamMember_cr3bp::toTraj(const SysData_cr3bp *sys){
	SimEngine sim;
	Traj_cr3bp traj(sys);
	sim.runSim(IC, TOF, &traj);
	return traj;
}//====================================================

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	\brief Copy an input family member into this one
 *	\param mem some other family member
 */
void FamMember_cr3bp::copyMe(const FamMember_cr3bp& mem){
	eigVals = mem.eigVals;
	eigVecs = mem.eigVecs;
	IC = mem.IC;
	JC = mem.JC;
	TOF = mem.TOF;
	xAmplitude = mem.xAmplitude;
	yAmplitude = mem.yAmplitude;
	zAmplitude = mem.zAmplitude;
	stm = mem.stm;
}//===================================================


}// END of Astrohelion namespace