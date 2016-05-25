/**
 *	@file tpat_fam_generator.cpp
 *	@brief Generate families of orbits
 *
 *	So far, this only applies to the CR3BP
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
 
#include "tpat_fam_generator.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_body_data.hpp"
#include "tpat_constraint.hpp"
#include "tpat_event.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_fam_cr3bp.hpp"
#include "tpat_famMember_cr3bp.hpp"
#include "tpat_linMotion_engine.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <Eigen/Dense>

#include <cmath>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief simple, do-nothing constructor
 */
TPAT_Fam_Generator::TPAT_Fam_Generator(){}

/** 
 *	@brief Copy constructor
 *	@param f a family generator reference
 */
TPAT_Fam_Generator::TPAT_Fam_Generator(const TPAT_Fam_Generator &f){
	copyMe(f);
}//====================================================

TPAT_Fam_Generator::~TPAT_Fam_Generator(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Copy operator
 *	@param f a family generator reference
 */
TPAT_Fam_Generator& TPAT_Fam_Generator::operator =(const TPAT_Fam_Generator &f){
	copyMe(f);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Methods
//-----------------------------------------------------

/**
 *	@brief Set the continuation type/method
 *	@param t continuation type
 */
void TPAT_Fam_Generator::setContType(TPAT_Continuation_Tp t){ contType = t; }

/**
 *	@brief Set the slope threshold
 *
 *	This quantity is used to determine which independent variable to fix. If 
 *	one variable is changeing quickly while the other is not, the algorithm changes
 *	the continuation scheme to vary the more slowly-changing variable, increasing
 *	accuracy and allowing the code to move around corners in the independent
 *	variable space.
 *
 *	@param d the slope threshold
 */
void TPAT_Fam_Generator::setSlopeThresh(double d){ slopeThresh = std::abs(d); }

/** 
 *	@brief Set the step size we take when performing simple continuation
 *
 *	The default value is 0.0005
 *	@param d the step size, non-dimensional units
 */
void TPAT_Fam_Generator::setStep_simple(double d){ step_simple = d; }

/** 
 *	@brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	@param d the step size, non-dimensional units
 */
void TPAT_Fam_Generator::setStep_fitted_1(double d){ step_fitted_1 = d; }

/** 
 *	@brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	@param d the step size, non-dimensional units
 */
void TPAT_Fam_Generator::setStep_fitted_2(double d){ step_fitted_2 = d; }

/**
 * @brief  Set the corrector tolerance for the family generator
 * @details This is the tolerance that a periodic orbit will be
 * judged by; if constraints are not met to this tolerance, no
 * periodic orbit will be returned. Note that a looser tolerance may
 * allow continuation to make more prorgress
 * 
 * @param t corrector tolerance, non-dimensional
 */
void TPAT_Fam_Generator::setTol(double t){ tol = t; }

/**
 *	@brief Set the number of nodes used for corrections processes
 *
 *	The default value is 3
 *	@param n the number of nodes
 */
void TPAT_Fam_Generator::setNumNodes(int n){ numNodes = n; }

/**
 *	@brief Set the maximum number of orbits for this family
 *
 *	The default value is 500
 *	@param n the number of orbits
 */
void TPAT_Fam_Generator::setNumOrbits(int n){ numOrbits = n; }

/**
 *  @brief Set the minimum step size for any of the stepping values
 * 
 *  @param s Nondimensional minimum step size
 */
void TPAT_Fam_Generator::setMinStepSize(double s){ minStepSize = s;}

/**
 *  @brief Set the maximum step size for any of the stepping values
 * 
 *  @param s Nondimensional maximum step size
 */
void TPAT_Fam_Generator::setMaxStepSize(double s){ maxStepSize = s;}

//-----------------------------------------------------
//      Operations and Utility
//-----------------------------------------------------

/** 
 *	@brief Copy this object from the specified guy
 *	@param f the source family generator
 */
void TPAT_Fam_Generator::copyMe(const TPAT_Fam_Generator &f){
	numOrbits = f.numOrbits;
	numSimple = f.numSimple;
	step_simple = f.step_simple;
	step_fitted_1 = f.step_fitted_1;
	step_fitted_2 = f.step_fitted_2;
	curveFitMem = f.curveFitMem;
	tol = f.tol;
}//====================================================

/**
 *	@brief Generate an Axial family in the CR3BP
 *
 *	The axial family is generated by finding the second Lyapunov bifurcation
 *	and perturbing it in the z-dot direction.
 *
 *	Tuning:
 *	
 *	- To generate the Northern Axial families, choose initStepSize > 0. In the
 *	  Earth-Moon system, a magnitude of 1e-4 works well
 *
 *	@param lyapFamFile the filepath to a lyapunov family file
 *	@param initStepSize the initial step-off in the z-dot direction
 *	@return an axial family
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateAxial(const char* lyapFamFile, double initStepSize){
	TPAT_Fam_CR3BP lyapFam(lyapFamFile);

	TPAT_Fam_CR3BP axialFam(lyapFam.getSysData());

	// Try to find bifurcations
	lyapFam.sortMembers();
	lyapFam.sortEigs();
	std::vector<int> bifs = lyapFam.findBifurcations();

	if(bifs.size() == 0){
		printErr("Could not locate any bifurcations in the Lyapunov family; extiting...\n");
		return axialFam;
	}

	if(bifs.size() != 3)
		printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!\n");

	std::vector<double> IC = lyapFam.getMember(bifs[1]).getIC();
	double period = lyapFam.getMember(bifs[1]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {5};	// force z-dot to be non-zero
	IC[5] += initStepSize;

	TPAT_Traj_CR3BP firstAxial = cr3bp_getPeriodic(axialFam.getSysDataPtr(), IC, period,
		numNodes, 1, TPAT_Mirror_Tp::MIRROR_X_AX_H, fixStates, tol);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		std::vector<int> indVars {5,4};	// begin stepping in z-dot, optionally use y-dot
		std::vector<int> depVars {0,6};	// Predict x and period with least squares
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_X_AX_H, TPAT_Mirror_Tp::MIRROR_X_AX_V};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(&axialFam, firstAxial, mirrorTypes, indVars, depVars, 1);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(firstAxial, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;
		std::vector<int> initDir {0, 0, 0, 0, 0, sign};
		cr3bp_pseudoArcCont(&axialFam, initGuess, TPAT_Mirror_Tp::MIRROR_X_AX_H, initDir);
	}
	return axialFam;
}//====================================================

/**
 *  @brief Generate a family of vertical orbits
 * 
 *  @param axialFamFile a pointer to a file containing the axial family at the same 
 *  collinear point as the desired vertical family
 *  @param initStepSize initial step size from the bifurcating axial orbit
 * 
 *  @return a family of vertical orbits
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateVertical(const char* axialFamFile, double initStepSize){
	TPAT_Fam_CR3BP axialFam(axialFamFile);

	TPAT_Fam_CR3BP vertFam(axialFam.getSysData());

	// Try to find bifurcations
	axialFam.sortMembers();
	axialFam.sortEigs();
	std::vector<int> bifs = axialFam.findBifurcations();

	if(bifs.size() == 0){
		printErr("Could not locate any bifurcations in the Axial family; exiting...\n");
		return vertFam;
	}

	if(bifs.size() > 2){
		printWarn("Axial family has more than 2 bifurcations... may be incorrect and yield unexpected results\n");
	}

	for(size_t i = 0; i < bifs.size(); i++){
		printf("Axial Bifurcation at orbit %d\n", bifs[i]);
	}

	std::vector<double> IC = axialFam.getMember(bifs[0]).getIC();
	double period = axialFam.getMember(bifs[0]).getTOF();

	// The axial family has ICs at the x-axis; I want the vertical family to have ICs at the top of their figure 8's,
	// so the first step is to integrate to that point
	TPAT_Sim_Engine sim;
	sim.addEvent(TPAT_Event(vertFam.getSysDataPtr(), TPAT_Event_Tp::XZ_PLANE, 0, true));	// Stop integrated at XZ plane, going opposite direction as initial state
	TPAT_Traj_CR3BP quarterArc(vertFam.getSysDataPtr());
	sim.runSim(IC, period/3, &quarterArc);	// 1/3 period should be long enough to fly 1/4 of the trajectory

	IC = quarterArc.getStateByIx(-1);

	int numNodes = 3;
	std::vector<int> fixStates {2}; // force z-dot to be non-zero
	IC[2] += initStepSize;

	TPAT_Traj_CR3BP firstVertical = cr3bp_getPeriodic(vertFam.getSysDataPtr(), IC, period,
		numNodes, 2, TPAT_Mirror_Tp::MIRROR_XZ, fixStates, tol);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		std::vector<int> indVars {2, 0};	// Begin stepping in z, optionally use x
		std::vector<int> depVars {5, 6}; 	// Predict y-dot and period with least squares
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(&vertFam, firstVertical, mirrorTypes, indVars, depVars, 2);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(firstVertical, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;

		std::vector<int> initDir {0, 0, sign, 0, 0, 0};
		cr3bp_pseudoArcCont(&vertFam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}

	return vertFam;
}//====================================================

/**
 *	@brief Generate a Halo family in the CR3BP
 *
 *	The halo family is generated by finding the first Lyapunov bifurcation and 
 *	perturbing it into a 3D family.
 *
 *	Tuning:
 *	
 *	- To generate the Northern Halo families, choose initStepSize > 0 for L2 and < 0
 *	for L1 and L3. In the Earth-Moon system, 20 km (5.2e-5 non-dim) is a good magnitude.
 *
 *	@param lyapFamFile the location of a Lyapunov family file. 
 *	@param initStepSize the size of the initial step away from the bifurcating
 *	Lyapunov orbit (non-dimensional units)
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateHalo(const char* lyapFamFile, double initStepSize){
	TPAT_Fam_CR3BP lyapFam(lyapFamFile);
	
	TPAT_Fam_CR3BP haloFam(lyapFam.getSysData());	

	// Try to find bifurcations
	lyapFam.sortMembers();
	lyapFam.sortEigs();
	std::vector<int> bifs = lyapFam.findBifurcations();

	if(bifs.size() == 0){
		printErr("Could not locate any bifurcations in the Lyapunov family; extiting...\n");
		return haloFam;
	}

	if(bifs.size() != 3)
		printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!\n");

	std::vector<double> IC = lyapFam.getMember(bifs[0]).getIC();
	double period = lyapFam.getMember(bifs[0]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {2};	// force z to be out of plane
	IC[2] += initStepSize;

	TPAT_Traj_CR3BP firstHalo = cr3bp_getPeriodic(haloFam.getSysDataPtr(), IC, period,
		numNodes, 1, TPAT_Mirror_Tp::MIRROR_XZ, fixStates, tol);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		std::vector<int> indVars {2,0};	// begin stepping in z, optionally using x
		std::vector<int> depVars {4};	// Predict y-dot with least-squares
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(&haloFam, firstHalo, mirrorTypes, indVars, depVars, 1);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){

		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(firstHalo, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;
		std::vector<int> initDir {0, 0, sign, 0, 0, 0};
		cr3bp_pseudoArcCont(&haloFam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}
	return haloFam;
}//=======================================================

/**
 *	@brief Generate a Lyapunov family in the CR3BP
 *
 *	Tuning:
 *
 *	- The L1 region is typically a little more sensitive, so try smaller step sizes
 *	- A good initial guess for r0 is 0.001; smaller values may be possible, but the
 *	  corrector often has trouble.
 *	
 *	@param sysData represents the system the Lyapunov exists in
 *	@param LPt The Lagrange point number [1-5]
 *	@param x0 the initial displacement from the Lagrange point along the x-axis.
 *
 *	@return a family of orbits
 *	@throws TPAT_Exception if <tt>LPt</tt> is invalid
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateLyap(TPAT_Sys_Data_CR3BP sysData, int LPt, double x0){
	if(LPt < 1 || LPt > 3)
		throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateLyap: Invalid LPt number");

	// Initialize variables and containers for data
	TPAT_Fam_CR3BP fam(sysData);

	// Get initial guess from linearization
	double LPt_data[] = {0,0,0};
	TPAT_Model_CR3BP::getEquilibPt(fam.getSysDataPtr(), LPt, 1e-14, LPt_data);

	// Begin solving - get linear approximation at ICs
	double r0[] = {x0, 0, 0};
	TPAT_LinMotion_Engine linEngine;
	TPAT_Traj_CR3BP linTraj = linEngine.getCR3BPLinear(LPt, r0,
		TPAT_LinMotion_Tp::ELLIP, fam.getSysDataPtr());

	fam.setSortType(TPAT_SortFam_Tp::SORT_X);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
	
		std::vector<int> indVars;
		indVars.push_back(0);	// We're going to fix the x-coordinate in the corrector to keep it from slipping
		indVars.push_back(4);	// Optionally, allow y-dot to be an independent variable if x is changing too quickly
		std::vector<int> depVars {4}; // Predict y-dot with least squares in the algorithm
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};
		cr3bp_natParamCont(&fam, linTraj, mirrorTypes, indVars, depVars, 1);

	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){

		// Make a copy of sysData so we can pass in a pointer
		TPAT_Sys_Data_CR3BP sys(sysData);

		// Get the initial state and tof from the linearization
		std::vector<double> IC = linTraj.getStateByIx(0);
		double tof = linTraj.getTimeByIx(-1);

		// Correct the initial guess to a true periodic orbit; we need a full DF matrix
		// for a CONVERGED family member to start PAC
		std::vector<int> fixStates {0};
		int order = 1;
		TPAT_Traj_CR3BP perOrbit = cr3bp_getPeriodic(&sys, IC, tof, numNodes, order, TPAT_Mirror_Tp::MIRROR_XZ, fixStates, tol);

		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(perOrbit, 2*numNodes-1);

		// Apply Pseudo Arclength Continuation: Ignore y (ix = 0) for periodicity, force y to equal 0 at node 0
		int sign = IC[0] - LPt_data[0] < 0 ? -1 : 1;	// force the first step to be away from Lagrange point
		std::vector<int> initDir {sign, 0, 0, 0, 0, 0};
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}

	return fam;
}//====================================================

/**
 *	@brief Generate a Butterfly family in the CR3BP
 *	
 *	@param sysData represents the system the Lyapunov exists in
 *	@param LPt The Lagrange point number [1-5]
 *
 *	@return a family of orbits
 *	@throws TPAT_Exception if <tt>LPt</tt> is not equal to two (others not implemented)
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateButterfly(TPAT_Sys_Data_CR3BP *sysData, int LPt){
	if(LPt != 2)
		throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_butterfly: LPts != 2 are not implemented");

	double LPt_data[] = {0,0,0};
	TPAT_Model_CR3BP::getEquilibPt(sysData, LPt, 1e-14, LPt_data);

	// The butterfly orbits bifurcate from the Halo Family, but I don't have good enough data
	// and/or bifurcation detection algorithms to find the proper bifurcation. For now,
	// use this IC for the bifucating Halo from Dan Grebow's Thesis
	double ic[] = {1.0406, 0, 0.1735, 0, -0.0770, 0};
	std::vector<double> icVec (ic, ic+6);
	double tof = 2.8077;

	printf("Correcting Butterfly...\n");
	// Correct to a periodic orbit
	std::vector<int> fixed {4};
	TPAT_Traj_CR3BP perOrbit = cr3bp_getPeriodic(sysData, icVec, tof, 8, 2, TPAT_Mirror_Tp::MIRROR_XZ, fixed, tol);

	printf("Creating Family...\n");
	// Initialize variables and containers for data
	TPAT_Fam_CR3BP fam(*sysData);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		// Butterfly-specific settings
		fam.setSortType(TPAT_SortFam_Tp::SORT_X);
		// std::vector<int> indVars {0,2};
		std::vector<int> indVars {0, 2};
		// std::vector<int> depVars {4,6};
		std::vector<int> depVars {4, 6};
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, 2);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(perOrbit, 2*numNodes-1);

		std::vector<int> initDir {1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}
	return fam;
}//====================================================

/**
 *  @brief Generate the Distant Retrograde Orbit (DRO) family 
 *  @details This family is initialized from a conic orbit with an orbital 
 *  radius that corresponds to the minimum flyby altitude, or, if that value is 
 *  less than 1 km, an altitude of 100 km.
 * 
 *  @param sysData System data object describing the CR3BP
 *  @return A family of DROs
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateDRO(TPAT_Sys_Data_CR3BP *sysData){
	TPAT_Body_Data P2Data = TPAT_Body_Data(sysData->getPrimary(1));
	double orbR = P2Data.getRadius() + 
		(P2Data.getMinFlyBy() > P2Data.getRadius() ? P2Data.getMinFlyBy() : P2Data.getRadius());	// minimum acceptable orbital radius, km
	double orbV = sqrt(P2Data.getGravParam()/orbR);							// Circular velocity at orbR, km/s
	double orbT = 2*PI*sqrt(pow(orbR, 3)/P2Data.getGravParam());					// Orbital period, sec

	double IC[] {1 - sysData->getMu() - orbR/sysData->getCharL(), 0, 0,
				 0, orbV*sysData->getCharT()/sysData->getCharL(), 0};		// IC for a DRO from the conic
	std::vector<double> icVec (IC, IC+6);

	// printf("Conic Arc State = [%.6f, 0, 0, 0, %.6f, 0], Period = %.6f\n", IC[0], IC[4], orbT/sysData->getCharT());

	// waitForUser();
	// Correct to be periodic
	printf("Correcting initial DRO from conic...\n");
	TPAT_Traj_CR3BP perOrbit = cr3bp_getPeriodic(sysData, icVec, orbT/sysData->getCharT(), TPAT_Mirror_Tp::MIRROR_XZ, tol);

	printf("Creating Family...\n");
	TPAT_Fam_CR3BP fam(*sysData);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		fam.setSortType(TPAT_SortFam_Tp::SORT_X);
		std::vector<int> indVars {0, 4};			// Vary x and vy
		std::vector<int> depVars {0, 4, 6};			// Predict x, vy, and period
		std::vector<TPAT_Mirror_Tp> mirrorTypes{TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		// Set the simple step size to be negative so that it moves away from P2
		if(step_simple > 0)
			step_simple *= -1;

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, 1);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		TPAT_Nodeset_CR3BP initGuess(perOrbit, 2*numNodes-1);
		std::vector<int> initDir {-1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}

	return fam;
}//====================================================

/**
 *  @brief Generate the Low Prograde Orbit (LPO) family 
 *  @details This family is initialized from a conic orbit with an orbital 
 *  radius that corresponds to the minimum flyby altitude, or, if that value is 
 *  less than 1 km, an altitude of 100 km.
 * 
 *  @param sysData System data object describing the CR3BP
 *  @return A family of DROs
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateLPO(TPAT_Sys_Data_CR3BP *sysData){
	TPAT_Body_Data P2Data = TPAT_Body_Data(sysData->getPrimary(1));
	double orbR = P2Data.getRadius() + 
		(P2Data.getMinFlyBy() > P2Data.getRadius() ? P2Data.getMinFlyBy() : P2Data.getRadius());	// minimum acceptable orbital radius, km
	double orbV = sqrt(P2Data.getGravParam()/orbR);							// Circular velocity at orbR, km/s
	double orbT = 2*PI*sqrt(pow(orbR, 3)/P2Data.getGravParam());					// Orbital period, sec

	double IC[] {1 - sysData->getMu() - orbR/sysData->getCharL(), 0, 0,
				 0, -orbV*sysData->getCharT()/sysData->getCharL(), 0};		// IC for a DRO from the conic
	std::vector<double> icVec (IC, IC+6);

	// printf("Conic Arc State = [%.6f, 0, 0, 0, %.6f, 0], Period = %.6f\n", IC[0], IC[4], orbT/sysData->getCharT());

	// waitForUser();
	// Correct to be periodic
	printf("Correcting initial LPO from conic...\n");
	TPAT_Traj_CR3BP perOrbit = cr3bp_getPeriodic(sysData, icVec, orbT/sysData->getCharT(), TPAT_Mirror_Tp::MIRROR_XZ, tol);

	printf("Creating Family...\n");
	TPAT_Fam_CR3BP fam(*sysData);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		fam.setSortType(TPAT_SortFam_Tp::SORT_X);
		std::vector<int> indVars {0, 4};			// Vary x and vy
		std::vector<int> depVars {0, 4, 6};			// Predict x, vy, and period
		std::vector<TPAT_Mirror_Tp> mirrorTypes{TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		// Set the simple step size to be negative so that it moves away from P2
		if(step_simple > 0)
			step_simple *= -1;

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, 1);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		TPAT_Nodeset_CR3BP initGuess(perOrbit, 2*numNodes-1);
		std::vector<int> initDir {-1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}

	return fam;
}//====================================================

/**
 *  @brief Compute a family of p:q resonant orbits
 * 
 *  @param sysData Earth-Moon CR3BP system data object (other systems not implemented)
 *  @param p Resonance ratio numerator; the orbit completes p revolutions in an inertial frame
 *  in the same amount of time as the CR3BP system completes q revolutions in an inertial frame.
 *  @param q Resonance ratio denominator
 *  @return A family of resonant orbits
 *  @throws TPAT_Exception of the resonance ratio p:q is not implemented or recognized
 */
TPAT_Fam_CR3BP TPAT_Fam_Generator::cr3bp_generateRes(TPAT_Sys_Data_CR3BP *sysData, int p, int q){
	double x = 0, vy = 0, T = 0;
	int order = 0;
	switch(p){
		case 1:
		{
			switch(q){
				case 1: x = 0.6339739688; vy = 0.8390456686; T = 5.43075997; order = 1; break;
				case 2: x = 0.5140368985; vy = 1.2740998965; T = 12.32249803; order = 2; break;
				case 3: x = 0.6097558619; vy = 1.0695212172; T = 18.52105931; order = 3; break;
				case 4: x = 0.6647328846; vy = 0.9735457114; T = 24.75827998; order = 4; break;
				case 5: x = 0.7165422857; vy = 0.8894024945; T = 30.98350896; order = 5; break;
				default: throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 2:
		{
			switch(q){
				case 1: x = 1.1237405523; vy = -0.8620117202; T = 6.04619177; order = 3; break;
				case 3: x = 0.5259693391; vy = 1.2027315054; T = 18.53612002; order = 3; break;
				case 5: x = 0.6502418226; vy = 0.9609312003; T = 31.00065761; order = 5; break;
				default: throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 3:
		{
			switch(q){
				case 1: x = 0.8525476977; vy = -0.4790301584; T = 6.35120050; order = 3; break;
				case 2: x = 1.0929978357; vy = -0.6758673511; T = 11.74717118; order = 2; break;
				case 4: x = 0.5522611666; vy = 1.1178451576; T = 24.74772340; order = 4; break;
				case 5: x = 0.6047408300; vy = 1.0072641104; T = 30.97096176; order = 5; break;
				default: throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 4:
		{
			switch(q){
				case 1: x = 0.531016298978725; vy = 0.529364977382337; T = 6.20854994765688; order = 3; break;
				case 3: x = 1.13067423947448; vy = -0.646972098815793; T = 17.99449516; order = 3; break;
				case 5: x = 0.7502059802; vy = 0.6870566313; T = 30.21938914; order = 3; break;
				default: throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 5:
		{
			switch(q){
				case 1: x = 0.5464336946; vy = 0.2531922385; T = 6.30641755; order = 4; break;
				case 2: x = 0.8484521141; vy = -0.2213536734; T = 13.03725695; order = 4; break;
				case 3: x = 1.1076664385; vy = -0.7339807287; T = 18.38756543; order = 5; break;
				case 4: x = 1.1514231706; vy = -0.6411329235; T = 24.29007804; order = 4; break;
			}
			break;
		}
		default: throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_generateRes: Unsupported resonance ratio\n");
	}

	std::vector<double> ic {x, 0, 0, 0, vy, 0};		// Inititial state orthogonal to XZ plane
	std::vector<int> fixed {0};						// Fix x to begin with
	numNodes = T > 4 ? floor(T/2) : 2;

	if(p == 4 && q == 3)
		numNodes *= 2;
	
	printf("Correcting %d:%d Resonant Orbit...\n", p, q);
	// Correct to a periodic orbit
	TPAT_Traj_CR3BP perOrbit = cr3bp_getPeriodic(sysData, ic, T, numNodes, order, TPAT_Mirror_Tp::MIRROR_XZ, fixed, tol);
	// perOrbit.saveToMat("resOrbit_initSoln.mat");
	// waitForUser();

	printf("Creating Family...\n");
	// Initialize variables and containers for data
	TPAT_Fam_CR3BP fam(*sysData);
	fam.setSortType(TPAT_SortFam_Tp::SORT_X);

	if(contType == TPAT_Continuation_Tp::NAT_PARAM){
		// Butterfly-specific settings
		fam.setSortType(TPAT_SortFam_Tp::SORT_X);
		// std::vector<int> indVars {0,2};
		std::vector<int> indVars {0, 4};
		// std::vector<int> depVars {4,6};
		std::vector<int> depVars {4, 6};
		std::vector<TPAT_Mirror_Tp> mirrorTypes {TPAT_Mirror_Tp::MIRROR_XZ, TPAT_Mirror_Tp::MIRROR_XZ};

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, order);
		
		// Run the other direction too
		step_simple *= -1;
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, order);
	}else if(contType == TPAT_Continuation_Tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		TPAT_Nodeset_CR3BP initGuess(perOrbit, 2*numNodes-1);

		std::vector<int> initDir {1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);

		// Run the other direction too
		initDir[0] *= -1;
		cr3bp_pseudoArcCont(&fam, initGuess, TPAT_Mirror_Tp::MIRROR_XZ, initDir);
	}

	return fam;
}

/**
 *	@brief Continue a family of orbits in the CR3BP
 *	
 * 	@param fam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	@param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	@param mirrorTypes a vector of variables that describe how the family mirrors in the rotating 
 *	reference frame. Each entry corresponds to an independent variable in <tt>indVarIx</tt>
 *	@param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	@param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 *	@param order the multiplicity or order of the family; i.e. the number of revs around the primary
 *	or system before the orbit repeats itself. For example, a Period-3 DRO has order 3, and a butterfly
 *	has order 2
 *	@throws TPAT_Exception if <tt>indVarIx</tt> has fewer than two elements
 *	@throws TPAT_Exception if <tt>mirrorTypes</tt> does not have the same size as <tt>indVarIx</tt>
 *	@throws TPAT_Exception if the eigenvalues of the monodromy matrix cannot be computed
 *	@throws TPAT_Exception if one of the indices stored in <tt>indVarIx</tt> or <tt>depVarIx</tt> is
 *	out of range
 */
void TPAT_Fam_Generator::cr3bp_natParamCont(TPAT_Fam_CR3BP *fam, TPAT_Traj_CR3BP initialGuess,
	std::vector<TPAT_Mirror_Tp> mirrorTypes, std::vector<int> indVarIx, std::vector<int> depVarIx, int order){

	TPAT_Sys_Data_CR3BP sys = fam->getSysData();

	if(indVarIx.size() < 2)
		throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_natParamCont: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_natParamCont: there must be an equal number of ind. vars and mirror types");

	int indVar1 = indVarIx[0];
	int indVar2 = indVarIx[1];
	TPAT_Mirror_Tp mirrorType = mirrorTypes[0];

	// Initially assume that we're fixing indVar1
	std::vector<int> fixStates;
	fixStates.push_back(indVar1);

	// Get info from the initial guess trajectory
	std::vector<double> IC = initialGuess.getStateByIx(0);
	double tof = initialGuess.getTimeByIx(-1);

	// Initialize counters and storage containers
	int orbitCount = 0;
	double indVarSlope = NAN;
	double deltaVar1 = 1;
	double deltaVar2 = 1;

	std::vector<TPAT_Traj_CR3BP> members;
	bool diverged = false;
	TPAT_MultShoot_Data *itData = NULL;

	while(orbitCount < numOrbits){
		TPAT_Traj_CR3BP perOrbit(&sys);
		try{
			printf("Guess for IC: [%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f] %.4f\n", IC[0], IC[1], IC[2], IC[3],
				IC[4], IC[5], tof);
			printf("Fix States: ");
			for(size_t i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			// Simulate the orbit
			perOrbit = cr3bp_getPeriodic(&sys, IC, tof, numNodes, order, mirrorType, fixStates, tol, itData);

			diverged = false;
			printf("Orbit %03d converged!\n", ((int)members.size()));
		}catch(TPAT_Diverge &e){
			diverged = true;
		}catch(TPAT_LinAlg_Err &e){
			printErr("There was a linear algebra error during family continuation...\n");
			break;
		}

		// Check for large changes in period to detect leaving family
		if(!diverged && orbitCount > 2){
			// difference in TOF; use abs() because corrector may employ reverse time and switch to forward time
			double dTOF = std::abs(perOrbit.getTimeByIx(-1)) - std::abs(members[members.size()-1].getTimeByIx(-1));
			double percChange = std::abs(dTOF/perOrbit.getTimeByIx(-1));
			if(percChange > 0.25){
				printf("percChange = %.4f\n", percChange);
				printWarn("Period jumped (now = %.5f)! Left the family! Trying smaller step size...\n", perOrbit.getTimeByIx(-1));
				diverged = true;
			}
		}

		if(diverged && orbitCount == 0){
			printErr("Could not converge on the first family member; try a smaller step size\n");
			break;
		}else if(diverged && orbitCount > 0){
			perOrbit = members.back();	// Use previous solution as the converged solution to get correct next guess

			if(orbitCount <= numSimple){
				if(step_simple > minStepSize){
					step_simple = step_simple/2 > minStepSize ? step_simple/2 : minStepSize;
					printf("  Decreased step size to %0.4e (min = %.4e)\n", step_simple, minStepSize);
				}else{
					printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}else{
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq > minStepSize){
					dq = dq/2 > minStepSize ? dq/2 : minStepSize;
					printf("  Decreased step size to %0.4e (min = %.4e)\n", dq, minStepSize);

					if(std::abs(indVarSlope) > slopeThresh)
						step_fitted_1 = dq;
					else
						step_fitted_2 = dq;
				}else{
					printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}
		}else{	// Did not diverge

			// Save the computed orbit
			members.push_back(perOrbit);
			orbitCount++;

			// Check to see if we should update the step size
			if(itData->count < 4 && orbitCount > numSimple){
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq < maxStepSize){
					dq = 2*dq < maxStepSize ? 2*dq : maxStepSize;
					printf("  Increased step size to %0.4e (max = %.4e)\n", dq, maxStepSize);
					if(std::abs(indVarSlope) > slopeThresh)
						step_fitted_1 = dq;
					else
						step_fitted_2 = dq;
				}
			}

			// Compute eigenvalues
			MatrixXRd mono = perOrbit.getSTMByIx(-1);

			double monoErr = std::abs(1.0 - mono.determinant());
			if(monoErr > 1e-5)
				printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);
			
			Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
		    if(eigensolver.info() != Eigen::Success)
		        throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_natParamCont: Could not compute eigenvalues of monodromy matrix");

		    Eigen::VectorXcd vals = eigensolver.eigenvalues();
		    std::vector<cdouble> eigVals(vals.data(), vals.data()+6);

			// For debugging:
			// printf("Eigenvalues:\n");
			// for(int i = 0; i < 6; i++){ printf("  %s\n", complexToStr(eigVals[i]).c_str()); }
			// waitForUser();

			// Add orbit to family
			TPAT_FamMember_CR3BP child(perOrbit);
			child.setEigVals(eigVals);
			fam->addMember(child);
		}

		// Create next initial guess
		tof = perOrbit.getTimeByIx(-1);

		if(orbitCount < numSimple){
			// Use simple continuation; copy the converged IC, step forward in the independent variable
			IC = perOrbit.getStateByIx(0);
			IC.at(indVar1) += step_simple;
		}else{

			// Compute the slope for the first time
			if(orbitCount == numSimple){
				deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVar1] - members[orbitCount-2].getStateByIx(0)[indVar1];
				deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVar2] - members[orbitCount-2].getStateByIx(0)[indVar2];
				indVarSlope = deltaVar1/deltaVar2;
			}

			// Use least-squares fitting to predict the value for the independent and dependent variables
			std::vector<double> prevStates;
			int first = ((int)members.size()) - curveFitMem < 0 ? 0 : ((int)members.size()) - curveFitMem;

			for(size_t n = first; n < members.size(); n++){
				std::vector<double> ic = members[n].getStateByIx(0);
				prevStates.insert(prevStates.end(), ic.begin(), ic.begin()+6);
				prevStates.push_back(members[n].getTimeByIx(-1));
				prevStates.push_back(members[n].getJacobiByIx(0));
			}

			// This will hold the input depVars plus the unused independent variable
			std::vector<int> allDepVars;
			std::vector<double> predictedIC;
			if(std::abs(indVarSlope) > slopeThresh){
				mirrorType = mirrorTypes[0];
				// Use continuation in indVar1
				IC.at(indVar1) = perOrbit.getStateByIx(0).at(indVar1) + TPAT_Util::sign(deltaVar1)*step_fitted_1;
				fixStates.clear();
				fixStates.push_back(indVar1);

				// Use Least Squares to predict dependent vars and unusued ind. vars
				allDepVars = depVarIx;
				if(std::find(depVarIx.begin(), depVarIx.end(), indVar2) == depVarIx.end()){
					allDepVars.push_back(indVar2);	// only add if it isn't already part of depVarIx
				}
				predictedIC = familyCont_LS(indVar1, IC.at(indVar1),
					allDepVars, prevStates);
			}else{
				mirrorType = mirrorTypes[1];
				// Use continuation in indVar2
				IC.at(indVar2) = perOrbit.getStateByIx(0).at(indVar2) + TPAT_Util::sign(deltaVar2)*step_fitted_2;
				fixStates.clear();
				fixStates.push_back(indVar2);

				// Use Least Squares to predict dependent vars and unusued ind. vars
				allDepVars = depVarIx;
				if(std::find(depVarIx.begin(), depVarIx.end(), indVar1) == depVarIx.end()){
					allDepVars.push_back(indVar1);	// only add if it isn't already part of depVarIx
				}
				predictedIC = familyCont_LS(indVar2, IC.at(indVar2),
					allDepVars, prevStates);
			}

			// Update IC with predicted variables
			for(size_t n = 0; n < allDepVars.size(); n++){
				int ix = allDepVars[n];
				if(ix < 6)
					IC[ix] = predictedIC[ix];
				else if(ix == 6)
					tof = predictedIC[ix];
				else
					throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_natParamCont: Cannot update independent variable; index out of range");
			}

			// Update slope
			deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVar1] - members[orbitCount-2].getStateByIx(0)[indVar1];
			deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVar2] - members[orbitCount-2].getStateByIx(0)[indVar2];
			indVarSlope = deltaVar1/deltaVar2;
		}
	}// end of while loop
}//==================================================

/**
 *	@brief Continue a family of orbits in the CR3BP
 *	
 *	This algorithm is specifically tailored to periodic orbits; see the commented out
 *	code about constraint method 1 for a more general approach that can be applied to 
 *	families of non-periodic orbits
 *	
 * 	@param fam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	@param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	@param mirrorType describes how this orbit mirrors
 *	@param initDir a vector that contains one non-zero element that indicates the sign of 
 *	the desired step for the family. The index corresponds to the state index. For example,
 *	if I wish the family to continue with a step in the negative z direction for the first node,
 *	I would input a vector of the form {0, 0, -1, 0, 0, 0, ...}. Technically, you can constrain
 * 	a step on any of the free variables, but only one step direction will be considered (the first one
 * 	as the algorithm reads through the vector)
 * 	@throws TPAT_Exception if the mirrorType is not recognized
 * 	@throws TPAT_Exception if the free variable vector contains fewer than six states
 * 	@throws TPAT_Exception if the eigenvalues of the monodromy matrix cannot be computed
 */
void TPAT_Fam_Generator::cr3bp_pseudoArcCont(TPAT_Fam_CR3BP *fam, TPAT_Nodeset_CR3BP initialGuess,
	TPAT_Mirror_Tp mirrorType, std::vector<int> initDir){

	// Check inputs (Only applies to Constraint Method 1)
	// if(periodicityIgnoreIx < 0 || periodicityIgnoreIx > 5)
	// 	throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Periodicity Ignore Index out of range");

	// if(fixToVal_ix < 0 || fixToVal_ix > 7)
	// 	throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_pseudoArcCont: FixToVal Index out of range");

	// TODO - Make these editable?
	double stepSize = 0.001;
	double maxStepSize = 0.5;
	double minStepSize = 1e-7;

	TPAT_Sys_Data_CR3BP sys = fam->getSysData();
	TPAT_Nodeset_CR3BP familyMember(initialGuess);	// Copy the input initial guess

	printf("Correcting Initial Guess...\n");
	
	// Clear constraints and add new ones
	familyMember.clearAllConstraints();

	/*	Constraint Method 1:
	 *	* Apply a periodicity constraint that forces the first and final node to be collocated,
	 *	  	ignoring one to avoid numerical troubles. It is best to ignore one of the planar 
	 *		position components, especially for planar families - ignoring z or z-dot can shift the
	 *		null vector to have non-zero elements in the z and z-dot spots, effectively stepping 
	 *		out of the plane, which is not desireable for a planar family...
	 *	* Constrain one extra state to be zero; I used something that makes sense for a perpendicular
	 *		plane crossing here (avoid constraining z or z-dot to be zero for planar families)
	 *
	 *	These were the descriptions for three input ints for this method:
	 *
	 *	@param periodicityIgnoreIx (int) index of a state variable to ignore when enforcing periodicity. It is best
	 *	to ignore one of the planar components (i.e. x or y, corresponding to indices 0 and 1, respectively)
	 *	@param fixToVal_ix (int) index of a variable to fix to <tt>fixToVal_val</tt> at the first node. If the 
	 *	index is between 0 and 5, it will constraint one of the usual state variables. An index of 6 will 
	 *	constrain total TOF, and an index of 7 will constrain Jacobi at the first node
	 *	@param fixToVal_val (double) the value to constrain state <tt>fixToVal_ix</tt> to
	 */
	 /*
	// Create a periodicity constraint
	double periodicConData[] = {0,0,0,0,0,0};
	periodicConData[periodicityIgnoreIx] = NAN;
	TPAT_Constraint periodicCon(TPAT_Constraint_Tp::MATCH_CUST, familyMember.getNumNodes()-1, periodicConData, 6);

	TPAT_Constraint extraCon;
	if(fixToVal_ix < 6){
		double extraCon_data[] = {NAN, NAN, NAN, NAN, NAN, NAN};
		extraCon_data[fixToVal_ix] = fixToVal_val;
		extraCon = TPAT_Constraint(TPAT_Constraint_Tp::STATE, 0, extraCon_data, 6);
	}else if(fixToVal_ix == 6){
		double val = fixToVal_val;
		extraCon = TPAT_Constraint(TPAT_Constraint_Tp::TOF, 0, &val, 1);
	}else if(fixToVal_ix == 7){
		double val = fixToVal_val;
		extraCon = TPAT_Constraint(TPAT_Constraint_Tp::JC, 0, &val, 1);
	}

	familyMember.addConstraint(periodicCon);
	familyMember.addConstraint(extraCon);
	*/

	/* Constraint Method 2:
	 *
	 *	Apply two constraints that enforce perpendicular crossings at the initial and final nodes
	 *	This gives MUCH better performance for the Lyapunov families
	 */
	double perpCross_data[] = {NAN,NAN,NAN,NAN,NAN,NAN};
	switch(mirrorType){
		case TPAT_Mirror_Tp::MIRROR_XZ:
			perpCross_data[1] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case TPAT_Mirror_Tp::MIRROR_YZ:
            perpCross_data[0] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case TPAT_Mirror_Tp::MIRROR_XY:
            perpCross_data[2] = 0;
			perpCross_data[3] = 0;
			perpCross_data[4] = 0;
            break;
        case TPAT_Mirror_Tp::MIRROR_X_AX_H:
        case TPAT_Mirror_Tp::MIRROR_X_AX_V:
        	perpCross_data[1] = 0;
			perpCross_data[2] = 0;
			perpCross_data[3] = 0;
            break;
        default:
            throw TPAT_Exception("Mirror type either not defined or not implemented");
	}
	// constrain perpendicular crossings and periodicity
	TPAT_Constraint perpCross1_Con(TPAT_Constraint_Tp::STATE, 0, perpCross_data, 6);
	TPAT_Constraint perpCross2_Con(TPAT_Constraint_Tp::STATE, familyMember.getNumNodes()-1, perpCross_data, 6);

	familyMember.addConstraint(perpCross1_Con);
	familyMember.addConstraint(perpCross2_Con);

	std::vector<TPAT_Constraint> constraints {perpCross1_Con, perpCross2_Con};

	// Correct the nodeset to retrieve a free-variable vector for a family member
	TPAT_Correction_Engine corrector;
	corrector.setVarTime(true);			// Variable time MUST be enabled for PAC
	corrector.setEqualArcTime(true);	// MUST use equal arc time to get propper # of constraints
	corrector.setTol(tol);
	corrector.setIgnoreCrash(true);		// Ignore crashes into primary
	TPAT_MultShoot_Data familyItData(&familyMember);
	TPAT_Nodeset_CR3BP perNodes(static_cast<const TPAT_Sys_Data_CR3BP *>(initialGuess.getSysData()));

	try{
		familyItData = corrector.multShoot(&familyMember, &perNodes);
	}catch(TPAT_Diverge &e){
		printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Could not converge initial guess!\n");
	}catch(TPAT_LinAlg_Err &e){
		printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
	}

	printf("Applying continuation to compute family...\n");
	
	// Initialize counters and storage containers
	int orbitCount = 0;
	Eigen::VectorXd convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);
	Eigen::VectorXd prevN(familyItData.totalFree, 1);
	
	std::vector<TPAT_Traj_CR3BP> members;

	while(orbitCount < numOrbits){

		/* 
		 *	The first iteration should have a DF matrix that is (n-1) x n, but all further iterations will
		 * 	have an extra row for the pseudo-arc-length constraint; we want to remove that row and take the
		 * 	nullspace of the submatrix
		 */
		std::vector<double> DF_data;
		if(familyItData.totalCons == familyItData.totalFree){
			DF_data.insert(DF_data.begin(), familyItData.DF.begin(), familyItData.DF.begin() + familyItData.totalFree * (familyItData.totalCons - 1));
		}else{
			DF_data = familyItData.DF;
		}

		// Compute null space of previously computed member's Jacobian Matrix
		MatrixXRd DF = Eigen::Map<MatrixXRd>(&(DF_data[0]), familyItData.totalFree-1, familyItData.totalFree);
		Eigen::FullPivLU<MatrixXRd> lu(DF);
		lu.setThreshold(1e-14);
		MatrixXRd N = lu.kernel();

		printf("DF has dimensions %ld x %ld\n", DF.rows(), DF.cols());
		// Check to make sure the IS a nullspace
		if(N.rows() == 1){
			printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Nullspace is zero-dimensional; cannot proceed...\n");
			return;
		}		

		// // For debugging, save nullspace vectors to file
		// char filename[16];
		// sprintf(filename, "N%02d.csv", orbitCount);
		// N.toCSV(filename);

		/**
		 *	Choose the nullspace vector that is closest to the previous one (which converged)
		 */
		printf("Choosing Nullspace Vector (%ldD, %ld elements)\n", N.cols(), N.rows());
		if(orbitCount == 0){
			if(N.cols() > 1){
				printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Nullspace is multidimensional on first iteration; unsure how to proceed...\n");
				return;
			}

			bool sameDir = true;
			for(size_t i = 0; i < initDir.size(); i++){
				// Find a non-zero element
				if(initDir[i] != 0 && i < (size_t)(N.rows())){
					// If signs are different, assume direction is different
					if(N(i,0)*initDir[i] < 0){
						sameDir = false;
						break;
					}
				}
			}

			// Reverse direction of nullspace
			if(!sameDir)
				N *= -1;

		}else{
			/* Make sure nullspace direction stays consistent by choosing the 
			 * null vector that is closest to the same direction as the previous one
			 */
			int best_ix = 0;	// index of the column of the best vector option
			int best_sign = 1;	// sign associated with best vector
			double best_angle = 181;	// the best (smallest) angle found
			for(int i = 0; i < N.cols(); i++){
				// Compute angle from dot product
				Eigen::VectorXd col_i = N.cols() > 1 ? N.col(i) : N;
				Eigen::VectorXd dotProd = prevN.transpose()*col_i / (prevN.norm()*N.norm());
				double angle = std::acos(dotProd(0));
				int sign = 1;

				// Flip the sign if the angle is greater than 90 and change the value to 180 - angle
				printf("dot product = %.4f\n", dotProd(0));
				printf("Angle = %.4f deg\n", angle*180/PI);
				if(angle > PI/2.0){
					printColor(CYAN, "Angle is %.4f > pi/2; changing sign and angle\n", angle);
					angle = PI - angle;
					sign = -1;
				}

				// Keep track of the best option
				if(angle < best_angle){
					best_angle = angle;
					best_sign = sign;
					best_ix = i;
				}
			}

			printf("best ix = %d, sign = %d\n", best_ix, best_sign);
			Eigen::VectorXd temp = N.cols() > 1 ? N.col(best_ix) : N;
			N = best_sign*temp;	// Apply sign change, if needed
		}

		prevN = N;	// Update memory
		printf("Chose N with first elements = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...]\n",
			N(0), N(1), N(2), N(3), N(4), N(5));

		TPAT_Nodeset_CR3BP newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData, constraints);

		/*
		 *	Apply multiple shooting to converge the new guess to be a member of the family
		 */
		bool killLoop = false;
		try{
			while(stepSize >= minStepSize){
				try{
					familyItData = corrector.multShoot(&newMember, &perNodes);

					// If convergence was fast, take bigger steps
					if(familyItData.count < 5){
						stepSize = stepSize*2 < maxStepSize ? stepSize*2 : maxStepSize;
						printColor(MAGENTA, "Increased Step Size to %.4e (max %.4e)!\n", stepSize, maxStepSize);
					}else if(familyItData.count > 10){
						// If convergence was slow, take smaller steps
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);
					}

					// Exit this loop; we converged!
					break;
				}catch(TPAT_Diverge &e){
					if(stepSize > minStepSize){
						printWarn("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Corrector diverged... trying smaller step size\n");

						// Decrease step size and try again
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

						// Re-Create the initial guess using the new step size
						newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData, constraints);
					}else{
						printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Could not converge new family member!\n");
						killLoop = true;
						break;
					}
				}
			}
		}catch(TPAT_LinAlg_Err &e){
			printErr("TPAT_Fam_Generator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
			killLoop = true;
		}

		if(killLoop)
			break;

		printf("Orbit %03d converged!\n", ((int)members.size()));

		
		if(familyItData.totalFree < 6)
			throw TPAT_Exception("TPAT_Fam_Generator::PAC algorithm expects at least 6 states in the free variable vector");
		// Check to see if the converged family vector is significantly different from previously computed family member
		// Note that only the first 6 states (the IC for the trajectory) is checked; differences in other nodes 
		// are assumed to be insignificant (if IC is the same, only possible change is change in TOF)
		Eigen::VectorXd newX_init = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), 6, 1);
		Eigen::VectorXd oldX_init = Eigen::Map<Eigen::VectorXd>(convergedFreeVarVec.data(), 6, 1);

		Eigen::VectorXd diff = newX_init - oldX_init;
		double diffX = diff.norm();
		printErr("||diff in X(1:6)|| = %.4e\n", diffX);

		if(diffX < tol){
			printErr("Solutions from pseudo-arc-length have ceased changing; ending continuation\n");
			break;
		}

		if(diffX > 2*maxStepSize){
			printErr("Solution changed by amount greater than maxStepSize; solution jumped, ending continuation\n");
			break;
		}
		
		// Save new converged family vector
		convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);

		// Convert converged nodeset to an orbit to save; TODO - could be improved to be much faster!
		TPAT_Traj_CR3BP perOrbit = TPAT_Traj_CR3BP::fromNodeset(perNodes);

		members.push_back(perOrbit);
		orbitCount++;

		// Compute eigenvalues
		MatrixXRd mono = perOrbit.getSTMByIx(-1);

		double monoErr = std::abs(1.0 - mono.determinant());
		if(monoErr > 1e-5)
			printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);
		
		Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
	    if(eigensolver.info() != Eigen::Success)
	        throw TPAT_Exception("TPAT_Fam_Generator::cr3bp_pseudoArcCont: Could not compute eigenvalues of monodromy matrix");

	    Eigen::VectorXcd vals = eigensolver.eigenvalues();
	    std::vector<cdouble> eigVals(vals.data(), vals.data()+6);

		// Add orbit to family
		TPAT_FamMember_CR3BP child(perOrbit);
		child.setEigVals(eigVals);
		fam->addMember(child);
	}// end of while loop
}//==================================================

/**
 *	@brief Create a nodeset that contains an initial guess for a family member
 *	using pseudo arclength continuation
 *
 *	@param convergedFreeVarVec a matrix containing the free variable vector of the previous
 *	(nearest) converged family member
 *	@param N a 1D nullspace vector that lies tangent to the family
 *	@param stepSize scales the size of the step by scaling the nullspace vector
 *	@param familyItData an TPAT_MultShoot_Data object containing corrections information about the
 *	previous (nearest) converged family member
 *	@param cons a vector of constraints to place on the nodeset
 */
TPAT_Nodeset_CR3BP TPAT_Fam_Generator::cr3bp_getNextPACGuess(Eigen::VectorXd convergedFreeVarVec,
	Eigen::VectorXd N, double stepSize, TPAT_MultShoot_Data familyItData, std::vector<TPAT_Constraint> cons){

	/**
	 *	Step forwards away from previously converged solution
	 */
	Eigen::VectorXd newFreeVarVec = convergedFreeVarVec + stepSize*N;
	double *X = newFreeVarVec.data();

	// Convert into a new nodeset (TODO: Make this more flexible by putting conversion code in a model?)
	const TPAT_Sys_Data_CR3BP *sys = static_cast<const TPAT_Sys_Data_CR3BP *>(familyItData.sysData);
	TPAT_Nodeset_CR3BP newMember(sys);

	for(int n = 0; n < familyItData.numNodes; n++){
		newMember.addNode(TPAT_Node(X+6*n, 0));
		if(n > 0)
			newMember.addSeg(TPAT_Segment(n-1, n, X[6*familyItData.numNodes]/(familyItData.numNodes - 1)));
	}

	// Add same constraints
	for(size_t c = 0; c < cons.size(); c++){
		newMember.addConstraint(cons[c]);
	}

	/* 
	 *	Form the Pseudo-Arclength Continuation constraint
	 */
	std::vector<double> pacCon_data = familyItData.X;
	// Append the null vector (i.e. step direction)
	pacCon_data.insert(pacCon_data.end(), N.data(), N.data()+N.rows());
	// Append the step size
	pacCon_data.insert(pacCon_data.end(), stepSize);
	// Create the actual constraint
	TPAT_Constraint pacCon(TPAT_Constraint_Tp::PSEUDOARC, familyItData.numNodes-1, pacCon_data);
	newMember.addConstraint(pacCon);

	// Outputs for debugging and sanity checks
	printf("New IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...] tof = %.4f\n",
		X[0], X[1], X[2], X[3], X[4], X[5], X[newFreeVarVec.rows()-1]);

	return newMember;
}//====================================================

/**
 *	@brief Reset all parameters to their default values
 */
void TPAT_Fam_Generator::reset(){
	numOrbits = 500;
	numSimple = 3;
	step_simple = 0.0005;
	step_fitted_1 = 0.005;
	step_fitted_2 = 0.005;
	curveFitMem = 5;
	numNodes = 3;
	slopeThresh = 1;
	tol = 1e-12;
}//====================================================

//