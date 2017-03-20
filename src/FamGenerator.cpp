/**
 *	\file FamGenerator.cpp
 *	\brief Generate families of orbits
 *
 *	So far, this only applies to the CR3BP
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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
 
#include "FamGenerator.hpp"

#include "AsciiOutput.hpp"
#include "BodyData.hpp"
#include "Constraint.hpp"
#include "Event.hpp"
#include "EigenDefs.hpp"
#include "Fam_cr3bp.hpp"
#include "FamMember_cr3bp.hpp"
#include "LinMotionEngine.hpp"
#include "MultShootData.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Utilities.hpp"

#include <Eigen/Dense>

#include <cmath>


namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief simple, do-nothing constructor
 */
FamGenerator::FamGenerator(){}

/** 
 *	\brief Copy constructor
 *	\param f a family generator reference
 */
FamGenerator::FamGenerator(const FamGenerator &f){
	copyMe(f);
}//====================================================

/**
 *  \brief Destructor
 */
FamGenerator::~FamGenerator(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	\brief Copy operator
 *	\param f a family generator reference
 */
FamGenerator& FamGenerator::operator =(const FamGenerator &f){
	copyMe(f);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Methods
//-----------------------------------------------------

/**
 *	\brief Set the continuation type/method
 *	\param t continuation type
 */
void FamGenerator::setContType(Continuation_tp t){ contType = t; }

/**
 *	\brief Set the slope threshold
 *
 *	This quantity is used to determine which independent variable to fix. If 
 *	one variable is changeing quickly while the other is not, the algorithm changes
 *	the continuation scheme to vary the more slowly-changing variable, increasing
 *	accuracy and allowing the code to move around corners in the independent
 *	variable space.
 *
 *	\param d the slope threshold
 */
void FamGenerator::setSlopeThresh(double d){ slopeThresh = std::abs(d); }

/** 
 *	\brief Set the step size we take when performing simple continuation
 *
 *	The default value is 0.0005
 *	\param d the step size, non-dimensional units
 */
void FamGenerator::setStep_simple(double d){ step_simple = d; }

/** 
 *	\brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	\param d the step size, non-dimensional units
 */
void FamGenerator::setStep_fitted_1(double d){ step_fitted_1 = d; }

/** 
 *	\brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	\param d the step size, non-dimensional units
 */
void FamGenerator::setStep_fitted_2(double d){ step_fitted_2 = d; }

/**
 * \brief  Set the corrector tolerance for the family generator
 * \details This is the tolerance that a periodic orbit will be
 * judged by; if constraints are not met to this tolerance, no
 * periodic orbit will be returned. Note that a looser tolerance may
 * allow continuation to make more prorgress
 * 
 * \param t corrector tolerance, non-dimensional
 */
void FamGenerator::setTol(double t){ tol = t; }

/**
 *	\brief Set the number of nodes used for corrections processes
 *
 *	The default value is 3
 *	\param n the number of nodes
 */
void FamGenerator::setNumNodes(int n){ numNodes = n; }

/**
 *	\brief Set the maximum number of orbits for this family
 *
 *	The default value is 500
 *	\param n the number of orbits
 */
void FamGenerator::setNumOrbits(int n){ numOrbits = n; }

/**
 *  \brief Set the minimum step size for any of the stepping values
 * 
 *  \param s Nondimensional minimum step size
 */
void FamGenerator::setMinStepSize(double s){ minStepSize = s;}

/**
 *  \brief Set the maximum step size for any of the stepping values
 * 
 *  \param s Nondimensional maximum step size
 */
void FamGenerator::setMaxStepSize(double s){ maxStepSize = s;}

//-----------------------------------------------------
//      Operations and Utility
//-----------------------------------------------------

/** 
 *	\brief Copy this object from the specified guy
 *	\param f the source family generator
 */
void FamGenerator::copyMe(const FamGenerator &f){
	Engine::copyBaseEngine(f);
	numOrbits = f.numOrbits;
	numSimple = f.numSimple;
	step_simple = f.step_simple;
	step_fitted_1 = f.step_fitted_1;
	step_fitted_2 = f.step_fitted_2;
	curveFitMem = f.curveFitMem;
	tol = f.tol;
}//====================================================

/**
 *	\brief Generate a Lyapunov family in the CR3BP
 *
 *	Tuning:
 *
 *	- The L1 region is typically a little more sensitive, so try smaller step sizes
 *	- A good initial guess for r0 is 0.001; smaller values may be possible, but the
 *	  corrector often has trouble.
 *	
 *	\param LPt The Lagrange point number [1-5]
 *	\param x0 the initial displacement from the Lagrange point along the x-axis.
 *	\param pFam pointer to system data object the orbit exists in
 *	
 *	\return a family of orbits
 *	\throws Exception if <tt>LPt</tt> is invalid
 */
void FamGenerator::cr3bp_generateLyap(int LPt, double x0, Fam_cr3bp *pFam){
	if(LPt < 1 || LPt > 3)
		throw Exception("FamGenerator::cr3bp_generateLyap: Invalid LPt number");

	// Get initial guess from linearization
	double LPt_data[] = {0,0,0};
	DynamicsModel_cr3bp::getEquilibPt(pFam->getSysDataPtr(), LPt, 1e-14, LPt_data);

	// Begin solving - get linear approximation at ICs
	double r0[] = {x0, 0, 0};
	LinMotionEngine linEngine;
	Traj_cr3bp linTraj = linEngine.getCR3BPLinear(LPt, r0,
		LinMotion_tp::ELLIP, pFam->getSysDataPtr());

	pFam->setSortType(FamSort_tp::SORT_X);

	if(contType == Continuation_tp::NAT_PARAM){
	
		std::vector<int> indVars;
		indVars.push_back(0);	// We're going to fix the x-coordinate in the corrector to keep it from slipping
		indVars.push_back(4);	// Optionally, allow y-dot to be an independent variable if x is changing too quickly
		std::vector<int> depVars {4}; // Predict y-dot with least squares in the algorithm
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};
		cr3bp_natParamCont(pFam, linTraj, mirrorTypes, indVars, depVars, 1);

	}else if(contType == Continuation_tp::PSEUDO_ARC){
		// Get the initial state and tof from the linearization
		std::vector<double> IC = linTraj.getStateByIx(0);
		double tof = linTraj.getTimeByIx(-1);

		// Correct the initial guess to a true periodic orbit; we need a full DF matrix
		// for a CONVERGED family member to start PAC
		std::vector<int> fixStates {0};
		int order = 1;
		Traj_cr3bp perOrbit = cr3bp_getPeriodic(pFam->getSysDataPtr(), IC, tof, numNodes, order, Mirror_tp::MIRROR_XZ, fixStates, tol);

		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);

		// Apply Pseudo Arclength Continuation: Ignore y (ix = 0) for periodicity, force y to equal 0 at node 0
		int sign = IC[0] - LPt_data[0] < 0 ? -1 : 1;	// force the first step to be away from Lagrange point
		std::vector<int> initDir {sign, 0, 0, 0, 0, 0};
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *	\brief Generate a Halo family in the CR3BP
 *
 *	The halo family is generated by finding the first Lyapunov bifurcation and 
 *	perturbing it into a 3D family.
 *
 *	Tuning:
 *	
 *	- To generate the Northern Halo families, choose initStepSize > 0 for L2 and < 0
 *	for L1 and L3. In the Earth-Moon system, 20 km (5.2e-5 non-dim) is a good magnitude.
 *
 *	\param lyapFamFile the location of a Lyapunov family file. 
 *	\param initStepSize the size of the initial step away from the bifurcating
 *	Lyapunov orbit (non-dimensional units)
 *	\param pHaloFam pointer to the halo orbit family object
 */
void FamGenerator::cr3bp_generateHalo(const char* lyapFamFile, double initStepSize, Fam_cr3bp *pHaloFam){
	Fam_cr3bp lyapFam(lyapFamFile);
	
	if(lyapFam.getSysData() != pHaloFam->getSysData())
		throw Exception("FamGenerator::cr3bp_generateHalo: Halo family must have same system data as lyapunov it bifurcates from");

	Fam_cr3bp haloFam(lyapFam.getSysData());	

	// Try to find bifurcations
	lyapFam.sortMembers();
	lyapFam.sortEigs();
	std::vector<int> bifs = lyapFam.findBifurcations();

	if(bifs.size() == 0){
		astrohelion::printErr("Could not locate any bifurcations in the Lyapunov family; extiting...\n");
	}

	if(bifs.size() != 3)
		astrohelion::printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!\n");

	std::vector<double> IC = lyapFam.getMember(bifs[0]).getIC();
	double period = lyapFam.getMember(bifs[0]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {2};	// force z to be out of plane
	IC[2] += initStepSize;

	Traj_cr3bp firstHalo = cr3bp_getPeriodic(pHaloFam->getSysDataPtr(), IC, period,
		numNodes, 1, Mirror_tp::MIRROR_XZ, fixStates, tol);

	if(contType == Continuation_tp::NAT_PARAM){
		std::vector<int> indVars {2,0};	// begin stepping in z, optionally using x
		std::vector<int> depVars {4};	// Predict y-dot with least-squares
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(pHaloFam, firstHalo, mirrorTypes, indVars, depVars, 1);
	}else if(contType == Continuation_tp::PSEUDO_ARC){

		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(firstHalo, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;
		std::vector<int> initDir {0, 0, sign, 0, 0, 0};
		cr3bp_pseudoArcCont(pHaloFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//=======================================================

/**
 *	\brief Generate an Axial family in the CR3BP
 *
 *	The axial family is generated by finding the second Lyapunov bifurcation
 *	and perturbing it in the z-dot direction.
 *
 *	Tuning:
 *	
 *	- To generate the Northern Axial families, choose initStepSize > 0. In the
 *	  Earth-Moon system, a magnitude of 1e-4 works well
 *
 *	\param lyapFamFile the filepath to a lyapunov family file
 *	\param initStepSize the initial step-off in the z-dot direction
 *	\param pAxialFam pointer to the axial orbit family object
 */
void FamGenerator::cr3bp_generateAxial(const char* lyapFamFile, double initStepSize, Fam_cr3bp *pAxialFam){
	Fam_cr3bp lyapFam(lyapFamFile);

	if(lyapFam.getSysData() != pAxialFam->getSysData()){
		throw Exception("FamGenerator::cr3bp_generateAxial: Axial family must have same system data as the lyapunov it bifurcates from");
	}

	// Try to find bifurcations
	lyapFam.sortMembers();
	lyapFam.sortEigs();
	std::vector<int> bifs = lyapFam.findBifurcations();

	if(bifs.size() == 0){
		astrohelion::printErr("Could not locate any bifurcations in the Lyapunov family; extiting...\n");
	}

	if(bifs.size() != 3)
		astrohelion::printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!\n");

	std::vector<double> IC = lyapFam.getMember(bifs[1]).getIC();
	double period = lyapFam.getMember(bifs[1]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {5};	// force z-dot to be non-zero
	IC[5] += initStepSize;

	Traj_cr3bp firstAxial = cr3bp_getPeriodic(pAxialFam->getSysDataPtr(), IC, period,
		numNodes, 1, Mirror_tp::MIRROR_X_AX_H, fixStates, tol);

	if(contType == Continuation_tp::NAT_PARAM){
		std::vector<int> indVars {5,4};	// begin stepping in z-dot, optionally use y-dot
		std::vector<int> depVars {0,6};	// Predict x and period with least squares
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_X_AX_H, Mirror_tp::MIRROR_X_AX_V};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(pAxialFam, firstAxial, mirrorTypes, indVars, depVars, 1);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(firstAxial, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;
		std::vector<int> initDir {0, 0, 0, 0, 0, sign};
		cr3bp_pseudoArcCont(pAxialFam, initGuess, Mirror_tp::MIRROR_X_AX_H, initDir);
	}
}//====================================================

/**
 *  \brief Generate a family of vertical orbits
 * 
 *  \param axialFamFile a pointer to a file containing the axial family at the same 
 *  collinear point as the desired vertical family
 *  \param initStepSize initial step size from the bifurcating axial orbit
 *  \param pVertFam pointer to the vertical orbit family object
 *  
 *  \return a family of vertical orbits
 */
void FamGenerator::cr3bp_generateVertical(const char* axialFamFile, double initStepSize, Fam_cr3bp *pVertFam){
	Fam_cr3bp axialFam(axialFamFile);

	if(axialFam.getSysData() != pVertFam->getSysData()){
		throw Exception("FamGenerator::cr3bp_generateVertical: Vertical family must have the same system data object as the axial family it bifurcates from");
	}

	// Try to find bifurcations
	axialFam.sortMembers();
	axialFam.sortEigs();
	std::vector<int> bifs = axialFam.findBifurcations();

	if(bifs.size() == 0){
		astrohelion::printErr("Could not locate any bifurcations in the Axial family; exiting...\n");
	}

	if(bifs.size() > 2){
		astrohelion::printWarn("Axial family has more than 2 bifurcations... may be incorrect and yield unexpected results\n");
	}

	for(unsigned int i = 0; i < bifs.size(); i++){
		printf("Axial Bifurcation at orbit %d\n", bifs[i]);
	}

	std::vector<double> IC = axialFam.getMember(bifs[0]).getIC();
	double period = axialFam.getMember(bifs[0]).getTOF();

	// The axial family has ICs at the x-axis; I want the vertical family to have ICs at the top of their figure 8's,
	// so the first step is to integrate to that point
	SimEngine sim;
	sim.addEvent(Event(Event_tp::XZ_PLANE, 0, true));	// Stop integrated at XZ plane, going opposite direction as initial state
	Traj_cr3bp quarterArc(pVertFam->getSysDataPtr());
	sim.runSim(IC, period/3, &quarterArc);	// 1/3 period should be long enough to fly 1/4 of the trajectory

	IC = quarterArc.getStateByIx(-1);

	int numNodes = 3;
	std::vector<int> fixStates {2}; // force z-dot to be non-zero
	IC[2] += initStepSize;

	Traj_cr3bp firstVertical = cr3bp_getPeriodic(pVertFam->getSysDataPtr(), IC, period,
		numNodes, 2, Mirror_tp::MIRROR_XZ, fixStates, tol);

	if(contType == Continuation_tp::NAT_PARAM){
		std::vector<int> indVars {2, 0};	// Begin stepping in z, optionally use x
		std::vector<int> depVars {5, 6}; 	// Predict y-dot and period with least squares
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		// Set the simple step size to be negative if the user inputs a negative step-off distance
		if(step_simple > 0 && initStepSize < 0)
			step_simple *= -1;

		cr3bp_natParamCont(pVertFam, firstVertical, mirrorTypes, indVars, depVars, 2);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(firstVertical, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;

		std::vector<int> initDir {0, 0, sign, 0, 0, 0};
		cr3bp_pseudoArcCont(pVertFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *	\brief Generate a Butterfly family in the CR3BP
 *	
 *	\param LPt The Lagrange point number [1-5]
 *	\param pFam pointer to a family object
 *
 *	\throws Exception if <tt>LPt</tt> is not equal to two (others not implemented)
 */
void FamGenerator::cr3bp_generateButterfly(int LPt, Fam_cr3bp *pFam){
	if(LPt != 2)
		throw Exception("FamGenerator::cr3bp_butterfly: LPts != 2 are not implemented");

	SysData_cr3bp *pSys = pFam->getSysDataPtr();
	double LPt_data[] = {0,0,0};
	DynamicsModel_cr3bp::getEquilibPt(pSys, LPt, 1e-14, LPt_data);

	// The butterfly orbits bifurcate from the Halo Family, but I don't have good enough data
	// and/or bifurcation detection algorithms to find the proper bifurcation. For now,
	// use this IC for the bifurcating Halo from Dan Grebow's Thesis
	double ic[] = {1.0406, 0, 0.1735, 0, -0.0770, 0};
	std::vector<double> icVec (ic, ic+6);
	double tof = 2.8077;

	printf("Correcting Butterfly...\n");
	// Correct to a periodic orbit
	std::vector<int> fixed {4};
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(pSys, icVec, tof, 8, 2, Mirror_tp::MIRROR_XZ, fixed, tol);

	printf("Creating Family...\n");
	if(contType == Continuation_tp::NAT_PARAM){
		// Butterfly-specific settings
		pFam->setSortType(FamSort_tp::SORT_X);
		// std::vector<int> indVars {0,2};
		std::vector<int> indVars {0, 2};
		// std::vector<int> depVars {4,6};
		std::vector<int> depVars {4, 6};
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, 2);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);

		std::vector<int> initDir {1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *  \brief Generate the Distant Retrograde Orbit (DRO) family 
 *  \details This family is initialized from a conic orbit with an orbital 
 *  radius that corresponds to the minimum flyby altitude, or, if that value is 
 *  less than 1 km, an altitude of 100 km.
 * 
 *  \param pFam pointer to a family object
 */
void FamGenerator::cr3bp_generateDRO(Fam_cr3bp *pFam){
	SysData_cr3bp *pSys = pFam->getSysDataPtr();
	BodyData P2Data = BodyData(pSys->getPrimID(1));
	double orbR = P2Data.getBodyRad() + 
		(P2Data.getMinFlyBy() > P2Data.getBodyRad() ? P2Data.getMinFlyBy() : P2Data.getBodyRad());	// minimum acceptable orbital radius, km
	double orbV = sqrt(P2Data.getGravParam()/orbR);							// Circular velocity at orbR, km/s
	double orbT = 2*PI*sqrt(pow(orbR, 3)/P2Data.getGravParam());					// Orbital period, sec

	double IC[] {1 - pSys->getMu() - orbR/pSys->getCharL(), 0, 0,
				 0, orbV*pSys->getCharT()/pSys->getCharL(), 0};		// IC for a DRO from the conic
	std::vector<double> icVec (IC, IC+6);

	// printf("Conic Arc State = [%.6f, 0, 0, 0, %.6f, 0], Period = %.6f\n", IC[0], IC[4], orbT/pSys->getCharT());

	// waitForUser();
	// Correct to be periodic
	printf("Correcting initial DRO from conic...\n");
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(pSys, icVec, orbT/pSys->getCharT(), Mirror_tp::MIRROR_XZ, tol);

	printf("Creating Family...\n");

	if(contType == Continuation_tp::NAT_PARAM){
		pFam->setSortType(FamSort_tp::SORT_X);
		std::vector<int> indVars {0, 4};			// Vary x and vy
		std::vector<int> depVars {0, 4, 6};			// Predict x, vy, and period
		std::vector<Mirror_tp> mirrorTypes{Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		// Set the simple step size to be negative so that it moves away from P2
		if(step_simple > 0)
			step_simple *= -1;

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, 1);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);
		std::vector<int> initDir {-1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *  \brief Generate the Low Prograde Orbit (LPO) family 
 *  \details This family is initialized from a conic orbit with an orbital 
 *  radius that corresponds to the minimum flyby altitude, or, if that value is 
 *  less than 1 km, an altitude of 100 km.
 * 
 *  \param pFam pointer to a family data object
 */
void FamGenerator::cr3bp_generateLPO(Fam_cr3bp *pFam){
	SysData_cr3bp *pSys = pFam->getSysDataPtr();
	BodyData P2Data = BodyData(pSys->getPrimID(1));
	double orbR = P2Data.getBodyRad() + 
		(P2Data.getMinFlyBy() > P2Data.getBodyRad() ? P2Data.getMinFlyBy() : P2Data.getBodyRad());	// minimum acceptable orbital radius, km
	double orbV = sqrt(P2Data.getGravParam()/orbR);							// Circular velocity at orbR, km/s
	double orbT = 2*PI*sqrt(pow(orbR, 3)/P2Data.getGravParam());					// Orbital period, sec

	double IC[] {1 - pSys->getMu() - orbR/pSys->getCharL(), 0, 0,
				 0, -orbV*pSys->getCharT()/pSys->getCharL(), 0};		// IC for a DRO from the conic
	std::vector<double> icVec (IC, IC+6);

	// printf("Conic Arc State = [%.6f, 0, 0, 0, %.6f, 0], Period = %.6f\n", IC[0], IC[4], orbT/pSys->getCharT());

	// waitForUser();
	// Correct to be periodic
	printf("Correcting initial LPO from conic...\n");
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(pSys, icVec, orbT/pSys->getCharT(), Mirror_tp::MIRROR_XZ, tol);

	printf("Creating Family...\n");

	if(contType == Continuation_tp::NAT_PARAM){
		pFam->setSortType(FamSort_tp::SORT_X);
		std::vector<int> indVars {0, 4};			// Vary x and vy
		std::vector<int> depVars {0, 4, 6};			// Predict x, vy, and period
		std::vector<Mirror_tp> mirrorTypes{Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		// Set the simple step size to be negative so that it moves away from P2
		if(step_simple > 0)
			step_simple *= -1;

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, 1);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);
		std::vector<int> initDir {-1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *  \brief Generate a family of Distance Prograde Orbits (DPOs)
 * 
 *  \param pFam pointer to a family object
 */
void FamGenerator::cr3bp_generateDPO(Fam_cr3bp *pFam){
	SysData_cr3bp *pSys = pFam->getSysDataPtr();

	BodyData P2Data = BodyData(pSys->getPrimID(1));
	double L1_pos[3];
	DynamicsModel_cr3bp::getEquilibPt(pSys, 1, 1e-12, L1_pos);

	// Approximate initial state: 41.53% of the way between P2 and L1
	// double orbR_nondim = 0.415362*(1 - pSys->getMu() - L1_pos[0]);
	double orbR_nondim = 0.43*(1 - pSys->getMu() - L1_pos[0]);
	// Approximate: velocity in circular orbit at orbR
	double v_circ = sqrt(P2Data.getGravParam()/(pSys->getCharL()*orbR_nondim));
	SysData_cr3bp sysSE("sun", "earth");

	// Approximate: initial velocity is about -1.17*v_circ in SE and -1.18*v_circ in EM, so
	// craft crude interpolation between the two
	double IC[] = {1 - pSys->getMu() - orbR_nondim, 0, 0,
				   0, -(1.7 + (2.5e-6)*(pSys->getMu()/sysSE.getMu()))*v_circ*pSys->getCharT()/pSys->getCharL(), 0};

	std::vector<double> icVec (IC, IC+6);

	// getPeriodic() function will use the mirror condition to stop integration at the mirror plane, so
	// I use a large TOF to make sure it gets there.
	printf("Correcting initial DPO from approximation...\n");
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(pSys, icVec, PI, Mirror_tp::MIRROR_XZ, tol);

	perOrbit.saveToMat("initial_dpo.mat");
	waitForUser();

	printf("Creating family...\n");
	if(contType == Continuation_tp::NAT_PARAM){
		pFam->setSortType(FamSort_tp::SORT_VY);
		std::vector<int> indVars {4, 0};		// Vary vy and x
		std::vector<int> depVars {0, 4, 6};		// Predict x, vy, and period
		std::vector<Mirror_tp> mirrorTypes{Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, 1);

		// run the other direction too
		step_simple *= -1;
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, 1);

	}else if(contType == Continuation_tp:: PSEUDO_ARC){
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);
		std::vector<int> initDir {0, 0, 0, 0, 1, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);

		// run the other direction too
		initDir[5] *= -1;
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *  \brief Compute a family of p:q resonant orbits
 * 
 *  \param p Resonance ratio numerator; the orbit completes p revolutions in an inertial frame
 *  in the same amount of time as the CR3BP system completes q revolutions in an inertial frame.
 *  \param q Resonance ratio denominator
 *  \param pFam pointer to a family object
 *  
 *  \throws Exception of the resonance ratio p:q is not implemented or recognized
 */
void FamGenerator::cr3bp_generateRes(int p, int q, Fam_cr3bp *pFam){
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
				default: throw Exception("FamGenerator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 2:
		{
			switch(q){
				case 1: x = 1.1237405523; vy = -0.8620117202; T = 6.04619177; order = 3; break;
				case 3: x = 0.5259693391; vy = 1.2027315054; T = 18.53612002; order = 3; break;
				case 5: x = 0.6502418226; vy = 0.9609312003; T = 31.00065761; order = 5; break;
				default: throw Exception("FamGenerator::cr3bp_generateRes: Unsupported resonance ratio\n");
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
				default: throw Exception("FamGenerator::cr3bp_generateRes: Unsupported resonance ratio\n");
			}
			break;
		}
		case 4:
		{
			switch(q){
				case 1: x = 0.531016298978725; vy = 0.529364977382337; T = 6.20854994765688; order = 3; break;
				case 3: x = 1.13067423947448; vy = -0.646972098815793; T = 17.99449516; order = 3; break;
				case 5: x = 0.7502059802; vy = 0.6870566313; T = 30.21938914; order = 3; break;
				default: throw Exception("FamGenerator::cr3bp_generateRes: Unsupported resonance ratio\n");
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
		default: throw Exception("FamGenerator::cr3bp_generateRes: Unsupported resonance ratio\n");
	}

	std::vector<double> ic {x, 0, 0, 0, vy, 0};		// Inititial state orthogonal to XZ plane
	std::vector<int> fixed {0};						// Fix x to begin with
	numNodes = T > 4 ? floor(T/2) : 2;

	if(p == 4 && q == 3)
		numNodes *= 2;
	
	SysData_cr3bp *pSys = pFam->getSysDataPtr();
	printf("Correcting %d:%d Resonant Orbit...\n", p, q);
	// Correct to a periodic orbit
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(pSys, ic, T, numNodes, order, Mirror_tp::MIRROR_XZ, fixed, tol);
	// perOrbit.saveToMat("resOrbit_initSoln.mat");
	// waitForUser();

	printf("Creating Family...\n");
	// Initialize variables and containers for data
	pFam->setSortType(FamSort_tp::SORT_X);

	if(contType == Continuation_tp::NAT_PARAM){
		// Butterfly-specific settings
		pFam->setSortType(FamSort_tp::SORT_X);
		// std::vector<int> indVars {0,2};
		std::vector<int> indVars {0, 4};
		// std::vector<int> depVars {4,6};
		std::vector<int> depVars {4, 6};
		std::vector<Mirror_tp> mirrorTypes {Mirror_tp::MIRROR_XZ, Mirror_tp::MIRROR_XZ};

		printf("Using natural parameter continuation...\n");
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, order);
		
		// Run the other direction too
		step_simple *= -1;
		cr3bp_natParamCont(pFam, perOrbit, mirrorTypes, indVars, depVars, order);
	}else if(contType == Continuation_tp::PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		Nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);

		std::vector<int> initDir {1, 0, 0, 0, 0, 0};
		printf("Using pseudo-arclength continuation...\n");
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);

		// Run the other direction too
		initDir[0] *= -1;
		cr3bp_pseudoArcCont(pFam, initGuess, Mirror_tp::MIRROR_XZ, initDir);
	}
}//====================================================

/**
 *  \brief Compute a family of periodic orbits using pseudo arclength continuation
 *  from a nodeset
 * 
 *  \param traj An initial guess for a periodic orbit
 *  \param mirrorType Condition describing the mirror symmetry exhibited by this family of periodic orbits
 *  \param initDir 6-element vector that indicates the initial step direction along the family. For example, to
 *  step along -xdot, use {0,0,0,-1,0,0} as initDir.
 *  \param pFam pointer to a family data object in which family member data will be stored
 */
void FamGenerator::cr3bp_pacFromTraj(Traj_cr3bp traj, Mirror_tp mirrorType, std::vector<int> initDir, Fam_cr3bp *pFam){
	Nodeset_cr3bp nodes(traj, numNodes);
	cr3bp_pseudoArcCont(pFam, nodes, mirrorType, initDir);
}//====================================================

/**
 *  \brief Compute a family of periodic orbits using pseudo arclength continuation
 *  from a nodeset
 * 
 *  \param nodes An initial guess for a periodic orbit
 *  \param mirrorType Condition describing the mirror symmetry exhibited by this family of periodic orbits
 *  \param initDir 6-element vector that indicates the initial step direction along the family. For example, to
 *  step along -xdot, use {0,0,0,-1,0,0} as initDir.
 *  \param pFam pointer to a family data object in which family member data will be stored
 */
void FamGenerator::cr3bp_pacFromNodeset(Nodeset_cr3bp nodes, Mirror_tp mirrorType, std::vector<int> initDir, Fam_cr3bp *pFam){
	cr3bp_pseudoArcCont(pFam, nodes, mirrorType, initDir);
}//====================================================

/**
 *	\brief Continue a family of orbits in the CR3BP via natural parameter continuation
 *	
 * 	\param fam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	\param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	\param mirrorTypes a vector of variables that describe how the family mirrors in the rotating 
 *	reference frame. Each entry corresponds to an independent variable in <tt>indVarIx</tt>
 *	\param indVarIx a vector containing the indices of the independent variables to be used. You MUST
 *	specify at least two; currently only two can be used. The first index in the vector will be used
 *	first in the continuation (using stupid-simple continuation), and the second will be toggled
 *	on later if the slope favors it.
 *	\param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 *	\param order the multiplicity or order of the family; i.e. the number of revs around the primary
 *	or system before the orbit repeats itself. For example, a Period-3 DRO has order 3, and a butterfly
 *	has order 2
 *	\throws Exception if <tt>indVarIx</tt> has fewer than two elements
 *	\throws Exception if <tt>mirrorTypes</tt> does not have the same size as <tt>indVarIx</tt>
 *	\throws Exception if the eigenvalues of the monodromy matrix cannot be computed
 *	\throws Exception if one of the indices stored in <tt>indVarIx</tt> or <tt>depVarIx</tt> is
 *	out of range
 */
void FamGenerator::cr3bp_natParamCont(Fam_cr3bp *fam, Traj_cr3bp initialGuess,
	std::vector<Mirror_tp> mirrorTypes, std::vector<int> indVarIx, std::vector<int> depVarIx, int order){

	SysData_cr3bp sys = fam->getSysData();

	if(indVarIx.size() < 2)
		throw Exception("FamGenerator::cr3bp_natParamCont: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw Exception("FamGenerator::cr3bp_natParamCont: there must be an equal number of ind. vars and mirror types");

	int indVar1 = indVarIx[0];
	int indVar2 = indVarIx[1];
	Mirror_tp mirrorType = mirrorTypes[0];

	// Initially assume that we're fixing indVar1
	std::vector<int> fixStates;
	fixStates.push_back(indVar1);

	// Get info from the initial guess trajectory
	std::vector<double> IC = initialGuess.getStateByIx(0);
	double tof = initialGuess.getTimeByIx(-1);
	double tof0 = tof;

	// Initialize counters and storage containers
	int orbitCount = 0;
	double indVarSlope = NAN;
	double deltaVar1 = 1;
	double deltaVar2 = 1;

	std::vector<Traj_cr3bp> members;
	bool diverged = false;

	// Create a dummy nodeset and create an iteration data object on the stack
	// The cr3bp_getPeriodic() function will only pass an iteration data pointer back
	// if the one passed in is not NULL, hence we create a valid object and delete it
	// before exiting the function
	Nodeset_cr3bp tempNodes(static_cast<const SysData_cr3bp *>(initialGuess.getSysData()));
	MultShootData *pItData = new MultShootData(&tempNodes);

	while(orbitCount < numOrbits){
		Traj_cr3bp perOrbit(&sys);
		try{
			printf("Guess for IC: [%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f] %.4f\n", IC[0], IC[1], IC[2], IC[3],
				IC[4], IC[5], tof);
			printf("Fix States: ");
			for(unsigned int i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			// Simulate the orbit
			perOrbit = cr3bp_getPeriodic(&sys, IC, tof, numNodes, order, mirrorType, fixStates, tol, pItData);

			diverged = false;
			printf("Orbit %03d converged!\n", static_cast<int>(members.size()));
		}catch(DivergeException &e){
			diverged = true;
		}catch(LinAlgException &e){
			astrohelion::printErr("There was a linear algebra error during family continuation...\n");
			break;
		}

		// Check for large changes in period to detect leaving family
		if(!diverged && orbitCount > 2){
			// difference in TOF; use abs() because corrector may employ reverse time and switch to forward time
			double dTOF = std::abs(perOrbit.getTimeByIx(-1)) - std::abs(members[members.size()-1].getTimeByIx(-1));
			double percChange = std::abs(dTOF/perOrbit.getTimeByIx(-1));
			if(percChange > 0.25){
				printWarn("percChange = %.4f\n", percChange);
				astrohelion::printWarn("Period jumped (now = %.5f)! Left the family! Trying smaller step size...\n", perOrbit.getTimeByIx(-1));
				diverged = true;
			}
		}

		if(diverged && orbitCount == 0){
			astrohelion::printErr("Could not converge on the first family member; try a smaller step size\n");
			break;
		}else if(diverged && orbitCount > 0){
			perOrbit = members.back();	// Use previous solution as the converged solution to get correct next guess

			if(orbitCount <= numSimple){
				if(step_simple > minStepSize){
					step_simple = step_simple/2 > minStepSize ? step_simple/2 : minStepSize;
					printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", step_simple, minStepSize);
				}else{
					astrohelion::printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}else{
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq > minStepSize){
					dq = dq/2 > minStepSize ? dq/2 : minStepSize;
					printColor(MAGENTA, "  Decreased step size to %0.4e (min = %.4e)\n", dq, minStepSize);

					if(std::abs(indVarSlope) > slopeThresh)
						step_fitted_1 = dq;
					else
						step_fitted_2 = dq;
				}else{
					astrohelion::printErr("Minimum step size reached, could not converge... exiting\n");
					break;
				}
			}
		}else{	// Did not diverge

			// Save the computed orbit
			members.push_back(perOrbit);
			orbitCount++;

			// Check to see if we should update the step size
			if(pItData->count < 4 && orbitCount > numSimple){
				double dq = std::abs(indVarSlope) > slopeThresh ? step_fitted_1 : step_fitted_2;

				if(dq < maxStepSize){
					dq = 2*dq < maxStepSize ? 2*dq : maxStepSize;
					printColor(MAGENTA, "  Increased step size to %0.4e (max = %.4e)\n", dq, maxStepSize);
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
				astrohelion::printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);
			
			// Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
		 //    if(eigensolver.info() != Eigen::Success){
		 //    	if(pItData){
			// 		delete(pItData);
			// 		pItData = nullptr;
			// 	}
		 //        throw Exception("FamGenerator::cr3bp_natParamCont: Could not compute eigenvalues of monodromy matrix");
		 //    }

		    // Eigen::VectorXcd vals = eigensolver.eigenvalues();
		    // std::vector<cdouble> eigVals(vals.data(), vals.data()+6);

			// For debugging:
			// printf("Eigenvalues:\n");
			// for(int i = 0; i < 6; i++){ printf("  %s\n", complexToStr(eigVals[i]).c_str()); }
			// waitForUser();

			// Add orbit to family
			FamMember_cr3bp child(perOrbit);
			// child.setEigVals(eigVals);
			fam->addMember(child);
		}

		// Create next initial guess
		tof = perOrbit.getTimeByIx(-1);

		if(tof*tof0 < 0){
			printErr("Time-of-Flight changed sign: ending continuation process\n");
			break;
		}

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
			int first = static_cast<int>(members.size()) - curveFitMem < 0 ? 0 : static_cast<int>(members.size()) - curveFitMem;

			for(unsigned int n = first; n < members.size(); n++){
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
				IC.at(indVar1) = perOrbit.getStateByIx(0).at(indVar1) + astrohelion::sign(deltaVar1)*step_fitted_1;
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
				IC.at(indVar2) = perOrbit.getStateByIx(0).at(indVar2) + astrohelion::sign(deltaVar2)*step_fitted_2;
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
			for(unsigned int n = 0; n < allDepVars.size(); n++){
				int ix = allDepVars[n];
				if(ix < 6)
					IC[ix] = predictedIC[ix];
				else if(ix == 6)
					tof = predictedIC[ix];
				else{
					if(pItData){
						delete(pItData);
						pItData = nullptr;
					}
					throw Exception("FamGenerator::cr3bp_natParamCont: Cannot update independent variable; index out of range");
				}
			}

			// Update slope
			deltaVar1 = members[orbitCount-1].getStateByIx(0)[indVar1] - members[orbitCount-2].getStateByIx(0)[indVar1];
			deltaVar2 = members[orbitCount-1].getStateByIx(0)[indVar2] - members[orbitCount-2].getStateByIx(0)[indVar2];
			indVarSlope = deltaVar1/deltaVar2;
		}
	}// end of while loop

	if(pItData){
		delete(pItData);
		pItData = nullptr;
	}
}//==================================================

/**
 *	\brief Continue a family of orbits in the CR3BP
 *	
 *	This algorithm is specifically tailored to periodic orbits; see the commented out
 *	code about constraint method 1 for a more general approach that can be applied to 
 *	families of non-periodic orbits
 *	
 * 	\param fam a pointer to a family object to store family members in; the family MUST have
 *	defined its system data object
 *	\param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	\param mirrorType describes how this orbit mirrors
 *	\param initDir a vector that contains one non-zero element that indicates the sign of 
 *	the desired step for the family. The index corresponds to the state index. For example,
 *	if I wish the family to continue with a step in the negative z direction for the first node,
 *	I would input a vector of the form {0, 0, -1, 0, 0, 0, ...}. Technically, you can constrain
 * 	a step on any of the free variables, but only one step direction will be considered (the first one
 * 	as the algorithm reads through the vector)
 * 	\throws Exception if the mirrorType is not recognized
 * 	\throws Exception if the free variable vector contains fewer than six states
 * 	\throws Exception if the eigenvalues of the monodromy matrix cannot be computed
 */
void FamGenerator::cr3bp_pseudoArcCont(Fam_cr3bp *fam, Nodeset_cr3bp initialGuess,
	Mirror_tp mirrorType, std::vector<int> initDir){

	// Check inputs (Only applies to Constraint Method 1)
	// if(periodicityIgnoreIx < 0 || periodicityIgnoreIx > 5)
	// 	throw Exception("FamGenerator::cr3bp_pseudoArcCont: Periodicity Ignore Index out of range");

	// if(fixToVal_ix < 0 || fixToVal_ix > 7)
	// 	throw Exception("FamGenerator::cr3bp_pseudoArcCont: FixToVal Index out of range");

	// TODO - Make these editable?
	double stepSize = 0.001;
	double maxStepSize = 0.5;
	double minStepSize = 1e-7;

	SysData_cr3bp sys = fam->getSysData();
	Nodeset_cr3bp familyMember(initialGuess);	// Copy the input initial guess
	double tof0 = initialGuess.getTotalTOF();

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
	 *	\param periodicityIgnoreIx (int) index of a state variable to ignore when enforcing periodicity. It is best
	 *	to ignore one of the planar components (i.e. x or y, corresponding to indices 0 and 1, respectively)
	 *	\param fixToVal_ix (int) index of a variable to fix to <tt>fixToVal_val</tt> at the first node. If the 
	 *	index is between 0 and 5, it will constraint one of the usual state variables. An index of 6 will 
	 *	constrain total TOF, and an index of 7 will constrain Jacobi at the first node
	 *	\param fixToVal_val (double) the value to constrain state <tt>fixToVal_ix</tt> to
	 */
	 /*
	// Create a periodicity constraint
	double periodicConData[] = {0,0,0,0,0,0};
	periodicConData[periodicityIgnoreIx] = NAN;
	Constraint periodicCon(Constraint_tp::MATCH_CUST, familyMember.getNumNodes()-1, periodicConData, 6);

	Constraint extraCon;
	if(fixToVal_ix < 6){
		double extraCon_data[] = {NAN, NAN, NAN, NAN, NAN, NAN};
		extraCon_data[fixToVal_ix] = fixToVal_val;
		extraCon = Constraint(Constraint_tp::STATE, 0, extraCon_data, 6);
	}else if(fixToVal_ix == 6){
		double val = fixToVal_val;
		extraCon = Constraint(Constraint_tp::TOF, 0, &val, 1);
	}else if(fixToVal_ix == 7){
		double val = fixToVal_val;
		extraCon = Constraint(Constraint_tp::JC, 0, &val, 1);
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
		case Mirror_tp::MIRROR_XZ:
			perpCross_data[1] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case Mirror_tp::MIRROR_YZ:
            perpCross_data[0] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case Mirror_tp::MIRROR_XY:
            perpCross_data[2] = 0;
			perpCross_data[3] = 0;
			perpCross_data[4] = 0;
            break;
        case Mirror_tp::MIRROR_X_AX_H:
        case Mirror_tp::MIRROR_X_AX_V:
        	perpCross_data[1] = 0;
			perpCross_data[2] = 0;
			perpCross_data[3] = 0;
            break;
        default:
            throw Exception("Mirror type either not defined or not implemented");
	}
	// constrain perpendicular crossings and periodicity
	Constraint perpCross1_Con(Constraint_tp::STATE, 0, perpCross_data, 6);
	Constraint perpCross2_Con(Constraint_tp::STATE, familyMember.getNumNodes()-1, perpCross_data, 6);

	familyMember.addConstraint(perpCross1_Con);
	familyMember.addConstraint(perpCross2_Con);

	std::vector<Constraint> constraints {perpCross1_Con, perpCross2_Con};

	// Correct the nodeset to retrieve a free-variable vector for a family member
	CorrectionEngine corrector;
	corrector.setVarTime(true);			// Variable time MUST be enabled for PAC
	corrector.setEqualArcTime(true);	// MUST use equal arc time to get propper # of constraints
	corrector.setTol(tol);
	corrector.setIgnoreCrash(true);		// Ignore crashes into primary
	corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	MultShootData familyItData(&familyMember);
	Nodeset_cr3bp tempNodes(static_cast<const SysData_cr3bp *>(initialGuess.getSysData()));

	// Initialize this nodeset outside the loop because the familyItData will end up with a pointer
	// to this nodeset after the multiple shooting processs; if the declaration is in the loop,
	// the nodeset is destroyed each iteration and the pointer ceases to be useful.
	Nodeset_cr3bp perNodes(static_cast<const SysData_cr3bp *>(initialGuess.getSysData()));

	Nodeset_cr3bp newMember(static_cast<const SysData_cr3bp *>(initialGuess.getSysData()));

	try{
		familyItData = corrector.multShoot(&familyMember, &tempNodes);
	}catch(DivergeException &e){
		astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: Could not converge initial guess!\n");
	}catch(LinAlgException &e){
		astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
	}

	printf("Applying continuation to compute family...\n");
	
	// Initialize counters and storage containers
	int orbitCount = 0;
	Eigen::VectorXd convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);
	Eigen::VectorXd prevN(familyItData.totalFree, 1);
	
	std::vector<Traj_cr3bp> members;

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
			astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: Nullspace is zero-dimensional; cannot proceed...\n");
			return;
		}		

		// // For debugging, save nullspace vectors to file
		// char filename[16];
		// sprintf(filename, "N%02d.csv", orbitCount);
		// N.astrohelion::toCSV(filename);

		/**
		 *	Choose the nullspace vector that is closest to the previous one (which converged)
		 */
		printf("Choosing Nullspace Vector (%ldD, %ld elements)\n", N.cols(), N.rows());
		if(orbitCount == 0){
			if(N.cols() > 1){
				astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: Nullspace is multidimensional on first iteration; unsure how to proceed...\n");
				return;
			}

			bool sameDir = true;
			for(unsigned int i = 0; i < initDir.size(); i++){
				// Find a non-zero element
				if(initDir[i] != 0 && i < static_cast<unsigned int>(N.rows())){
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
					astrohelion::printColor(CYAN, "Angle is %.4f > pi/2; changing sign and angle\n", angle);
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

		// Nodeset_cr3bp newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData);
		newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData);
		// Reset perNodes
		perNodes = Nodeset_cr3bp(static_cast<const SysData_cr3bp *>(initialGuess.getSysData()));

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
						astrohelion::printColor(MAGENTA, "Increased Step Size to %.4e (max %.4e)!\n", stepSize, maxStepSize);
					}else if(familyItData.count > 10){
						// If convergence was slow, take smaller steps
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						astrohelion::printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);
					}

					// Exit this loop; we converged!
					break;
				}catch(DivergeException &e){
					if(stepSize > minStepSize){
						astrohelion::printWarn("FamGenerator::cr3bp_pseudoArcCont: Corrector diverged... trying smaller step size\n");

						// Decrease step size and try again
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						astrohelion::printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

						// Re-Create the initial guess using the new step size
						newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData);
					}else{
						astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: Could not converge new family member!\n");
						killLoop = true;
						break;
					}
				}
			}
		}catch(LinAlgException &e){
			astrohelion::printErr("FamGenerator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
			killLoop = true;
		}

		if(killLoop)
			break;

		printf("Orbit %03d converged!\n", static_cast<int>(members.size()));

		
		if(familyItData.totalFree < 6)
			throw Exception("FamGenerator::PAC algorithm expects at least 6 states in the free variable vector");
		// Check to see if the converged family vector is significantly different from previously computed family member
		// Note that only the first 6 states (the IC for the trajectory) is checked; differences in other nodes 
		// are assumed to be insignificant (if IC is the same, only possible change is change in TOF)
		Eigen::VectorXd newX_init = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), 6, 1);
		Eigen::VectorXd oldX_init = Eigen::Map<Eigen::VectorXd>(convergedFreeVarVec.data(), 6, 1);

		Eigen::VectorXd diff = newX_init - oldX_init;
		double diffX = diff.norm();
		astrohelion::printErr("||diff in X(1:6)|| = %.4e\n", diffX);

		if(diffX < tol){
			astrohelion::printErr("Solutions from pseudo-arc-length have ceased changing; ending continuation\n");
			break;
		}

		if(diffX > 2*maxStepSize){
			astrohelion::printErr("Solution changed by amount greater than maxStepSize; solution jumped, ending continuation\n");
			break;
		}
		
		// Save new converged family vector
		convergedFreeVarVec = Eigen::Map<Eigen::VectorXd>(&(familyItData.X[0]), familyItData.totalFree, 1);

		// Convert converged nodeset to an orbit to save; TODO - could be improved to be much faster!
		Traj_cr3bp perOrbit = Traj_cr3bp::fromNodeset(perNodes);
		// perNodes.saveToMat("temp_perNodes_pac.mat");
		// newMember.saveToMat("temp_newMember.mat");
		// perOrbit.saveToMat("temp_perOrbit_pac.mat");

		if(perOrbit.getTotalTOF()*tof0 < 0){
			printErr("Time-of-Flight changed sign; ending continuation\n");
			break;
		}

		members.push_back(perOrbit);
		orbitCount++;

		// Compute eigenvalues
		MatrixXRd mono = perOrbit.getSTMByIx(-1);
		
		double monoErr = std::abs(1.0 - mono.determinant());
		if(monoErr > 1e-5)
			astrohelion::printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);
		
		// Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
	 //    if(eigensolver.info() != Eigen::Success)
	 //        throw Exception("FamGenerator::cr3bp_pseudoArcCont: Could not compute eigenvalues of monodromy matrix");

	 //    Eigen::VectorXcd vals = eigensolver.eigenvalues();
	 //    std::vector<cdouble> eigVals(vals.data(), vals.data()+6);

		// Add orbit to family
		FamMember_cr3bp child(perOrbit);
		// child.setEigVals(eigVals);
		fam->addMember(child);
	}// end of while loop
}//==================================================

/**
 *	\brief Create a nodeset that contains an initial guess for a family member
 *	using pseudo arclength continuation
 *
 *	\param convergedFreeVarVec a matrix containing the free variable vector of the previous
 *	(nearest) converged family member
 *	\param N a 1D nullspace vector that lies tangent to the family
 *	\param stepSize scales the size of the step by scaling the nullspace vector
 *	\param pFamilyItData pointer to a MultShootData object containing corrections information about the
 *	previous (nearest) converged family member
 */
Nodeset_cr3bp FamGenerator::cr3bp_getNextPACGuess(Eigen::VectorXd convergedFreeVarVec,
	Eigen::VectorXd N, double stepSize, MultShootData pFamilyItData){

	/**
	 *	Step forwards away from previously converged solution
	 */

	Eigen::VectorXd newFreeVarVec = convergedFreeVarVec + stepSize*N;
	double *X = newFreeVarVec.data();

	// Convert into a new nodeset (TODO: Make this more flexible by putting conversion code in a model?)
	const SysData_cr3bp *sys = static_cast<const SysData_cr3bp *>(pFamilyItData.sysData);
	Nodeset_cr3bp newMember(sys);

	sys->getDynamicsModel()->multShoot_createOutput(&pFamilyItData, pFamilyItData.nodeset, false, &newMember);

	// Get rid of any pre-existing pseudo arclength constraints from previous corrections
	std::vector<Constraint> arcCons = newMember.getArcConstraints();
	newMember.clearArcConstraints();
	for(Constraint &con : arcCons){
		if(con.getType() != Constraint_tp::PSEUDOARC)
			newMember.addConstraint(con);
	}

	/* 
	 *	Form the Pseudo-Arclength Continuation constraint
	 */
	std::vector<double> pacCon_data = pFamilyItData.X;
	// Append the null vector (i.e. step direction)
	pacCon_data.insert(pacCon_data.end(), N.data(), N.data()+N.rows());
	// Append the step size
	pacCon_data.insert(pacCon_data.end(), stepSize);
	// Create the actual constraint
	Constraint pacCon(Constraint_tp::PSEUDOARC, pFamilyItData.numNodes-1, pacCon_data);
	newMember.addConstraint(pacCon);

	// Outputs for debugging and sanity checks
	printf("New IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...] tof = %.4f\n",
		X[0], X[1], X[2], X[3], X[4], X[5], X[newFreeVarVec.rows()-1]);

	return newMember;
}//====================================================

/**
 *	\brief Reset all parameters to their default values
 */
void FamGenerator::reset(){
	if(!bIsClean)
		cleanEngine();

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

/**
 *  \brief Utility function to reset any parameters and flags that are
 *  assigned during generation processes
 */
void FamGenerator::cleanEngine(){
	// Nothing stored, so nothing to do
	bIsClean = true;
}//====================================================

}// END of Astrohelion namespace