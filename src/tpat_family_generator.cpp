/**
 *	@file tpat_family_generator.cpp
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
#include "tpat.hpp"
 
#include "tpat_family_generator.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_linear_motion_engine.hpp"
#include "tpat_matrix.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <cmath>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief simple, do-nothing constructor
 */
tpat_family_generator::tpat_family_generator(){}

/** 
 *	@brief Copy constructor
 *	@param f a family generator reference
 */
tpat_family_generator::tpat_family_generator(const tpat_family_generator &f){
	copyMe(f);
}//====================================================

tpat_family_generator::~tpat_family_generator(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Copy operator
 *	@param f a family generator reference
 */
tpat_family_generator& tpat_family_generator::operator =(const tpat_family_generator &f){
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
void tpat_family_generator::setContType(cont_t t){ contType = t; }

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
void tpat_family_generator::setSlopeThresh(double d){ slopeThresh = std::abs(d); }

/** 
 *	@brief Set the step size we take when performing simple continuation
 *
 *	The default value is 0.0005
 *	@param d the step size, non-dimensional units
 */
void tpat_family_generator::setStep_simple(double d){ step_simple = d; }

/** 
 *	@brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	@param d the step size, non-dimensional units
 */
void tpat_family_generator::setStep_fitted_1(double d){ step_fitted_1 = d; }

/** 
 *	@brief Set the step size we take when performing advanced continuation
 *	using independent variable #1
 *
 *	The default value is 0.005
 *	@param d the step size, non-dimensional units
 */
void tpat_family_generator::setStep_fitted_2(double d){ step_fitted_2 = d; }

/**
 * @brief  Set the corrector tolerance for the family generator
 * @details This is the tolerance that a periodic orbit will be
 * judged by; if constraints are not met to this tolerance, no
 * periodic orbit will be returned. Note that a looser tolerance may
 * allow continuation to make more prorgress
 * 
 * @param t corrector tolerance, non-dimensional
 */
void tpat_family_generator::setTol(double t){ tol = t; }

/**
 *	@brief Set the number of nodes used for corrections processes
 *
 *	The default value is 3
 *	@param n the number of nodes
 */
void tpat_family_generator::setNumNodes(int n){ numNodes = n; }

/**
 *	@brief Set the maximum number of orbits for this family
 *
 *	The default value is 500
 *	@param n the number of orbits
 */
void tpat_family_generator::setNumOrbits(int n){ numOrbits = n; }

//-----------------------------------------------------
//      Operations and Utility
//-----------------------------------------------------

/** 
 *	@brief Copy this object from the specified guy
 *	@param f the source family generator
 */
void tpat_family_generator::copyMe(const tpat_family_generator &f){
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
tpat_family_cr3bp tpat_family_generator::cr3bp_generateAxial(const char* lyapFamFile, double initStepSize){
	tpat_family_cr3bp lyapFam(lyapFamFile);

	tpat_family_cr3bp axialFam(lyapFam.getSysData());
	axialFam.setSortType(tpat_family_cr3bp::SORT_X);

	// Set the simple step size to be negative if the user inputs a negative step-off distance
	if(step_simple > 0 && initStepSize < 0)
		step_simple *= -1;


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

	tpat_traj_cr3bp firstAxial = cr3bp_getPeriodic(axialFam.getSysDataPtr(), IC, period,
		numNodes, 1, MIRROR_X_AX_H, fixStates, tol);

	std::vector<int> indVars {5,4};	// begin stepping in z-dot, optionally use y-dot
	std::vector<int> depVars {0,6};	// Predict x and period with least squares
	std::vector<mirror_t> mirrorTypes {MIRROR_X_AX_H, MIRROR_X_AX_V};

	cr3bp_natParamCont(&axialFam, firstAxial, mirrorTypes, indVars, depVars, 1);

	return axialFam;
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
tpat_family_cr3bp tpat_family_generator::cr3bp_generateHalo(const char* lyapFamFile, double initStepSize){
	tpat_family_cr3bp lyapFam(lyapFamFile);
	
	tpat_family_cr3bp haloFam(lyapFam.getSysData());
	haloFam.setSortType(tpat_family_cr3bp::SORT_X);

	// Set the simple step size to be negative if the user inputs a negative step-off distance
	if(step_simple > 0 && initStepSize < 0)
		step_simple *= -1;

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

	tpat_traj_cr3bp firstHalo = cr3bp_getPeriodic(haloFam.getSysDataPtr(), IC, period,
		numNodes, 1, MIRROR_XZ, fixStates, tol);

	if(contType == NAT_PARAM){
		std::vector<int> indVars {2,0};	// begin stepping in z, optionally using x
		std::vector<int> depVars {4};	// Predict y-dot with least-squares
		std::vector<mirror_t> mirrorTypes {MIRROR_XZ, MIRROR_XZ};

		cr3bp_natParamCont(&haloFam, firstHalo, mirrorTypes, indVars, depVars, 1);
	}else if(contType == PSEUDO_ARC){

		// Turn trajectory object into nodeset; double number of nodes
		tpat_nodeset_cr3bp initGuess(firstHalo, 2*numNodes-1);

		int sign = initStepSize < 0 ? -1 : 1;
		std::vector<int> initDir {0, 0, sign, 0, 0, 0};
		cr3bp_pseudoArcCont(&haloFam, initGuess, MIRROR_XZ, initDir);
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
 */
tpat_family_cr3bp tpat_family_generator::cr3bp_generateLyap(tpat_sys_data_cr3bp sysData, int LPt, double x0){
	if(LPt < 1 || LPt > 3)
		throw tpat_exception("tpat_family_generator::cr3bp_generateLyap: Invalide LPt number");

	// Initialize variables and containers for data
	tpat_family_cr3bp fam(sysData);

	// Get initial guess from linearization
	double LPt_data[] = {0,0,0};
	cr3bp_getEquilibPt(sysData, LPt, 1e-14, LPt_data);

	// Begin solving - get linear approximation at ICs
	double r0[] = {x0, 0, 0};
	tpat_linear_motion_engine linEngine;
	tpat_traj_cr3bp linTraj = linEngine.getCR3BPLinear(LPt, r0,
		tpat_linear_motion_engine::ELLIP, fam.getSysDataPtr());

	fam.setSortType(tpat_family_cr3bp::SORT_X);

	if(contType == NAT_PARAM){
	
		std::vector<int> indVars;
		indVars.push_back(0);	// We're going to fix the x-coordinate in the corrector to keep it from slipping
		indVars.push_back(4);	// Optionally, allow y-dot to be an independent variable if x is changing too quickly
		std::vector<int> depVars {4}; // Predict y-dot with least squares in the algorithm
		std::vector<mirror_t> mirrorTypes {MIRROR_XZ, MIRROR_XZ};
		cr3bp_natParamCont(&fam, linTraj, mirrorTypes, indVars, depVars, 1);

	}else if(contType == PSEUDO_ARC){

		// Make a copy of sysData so we can pass in a pointer
		tpat_sys_data_cr3bp sys(sysData);

		// Get the initial state and tof from the linearization
		std::vector<double> IC = linTraj.getState(0);
		double tof = linTraj.getTime(-1);

		// Correct the initial guess to a true periodic orbit; we need a full DF matrix
		// for a CONVERGED family member to start PAC
		std::vector<int> fixStates {0};
		int order = 1;
		tpat_traj_cr3bp perOrbit = cr3bp_getPeriodic(&sys, IC, tof, numNodes, order, MIRROR_XZ, fixStates, tol);

		// Turn trajectory object into nodeset; double number of nodes
		tpat_nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);

		// Apply Pseudo Arclength Continuation: Ignore y (ix = 0) for periodicity, force y to equal 0 at node 0
		int sign = IC[0] - LPt_data[0] < 0 ? -1 : 1;	// force the first step to be away from Lagrange point
		std::vector<int> initDir {sign, 0, 0, 0, 0, 0};
		cr3bp_pseudoArcCont(&fam, initGuess, MIRROR_XZ, initDir);
	}

	return fam;
}//====================================================

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
 */
tpat_family_cr3bp tpat_family_generator::cr3bp_generateButterfly(tpat_sys_data_cr3bp *sysData, int LPt){
	if(LPt != 2)
		throw tpat_exception("tpat_family_generator::cr3bp_butterfly: LPts != 2 are not implemented");

	double LPt_data[] = {0,0,0};
	cr3bp_getEquilibPt(*sysData, LPt, 1e-14, LPt_data);

	// The butterfly orbits bifurcate from the Halo Family, but I don't have good enough data
	// and/or bifurcation detection algorithms to find the proper bifurcation. For now,
	// use this IC for the bifucating Halo from Dan Grebow's Thesis
	double ic[] = {1.0406, 0, 0.1735, 0, -0.0770, 0};
	std::vector<double> icVec (ic, ic+6);
	double tof = 2.8077;

	printf("Correcting Halo...\n");
	// Correct to a periodic orbit
	std::vector<int> fixed {4};
	tpat_traj_cr3bp perOrbit = cr3bp_getPeriodic(sysData, icVec, tof, 8, 2, MIRROR_XZ, fixed, tol);

	printf("Creating Family...\n");
	// Initialize variables and containers for data
	tpat_family_cr3bp fam(*sysData);

	if(contType == NAT_PARAM){
		// Butterfly-specific settings
		fam.setSortType(tpat_family_cr3bp::SORT_X);
		// std::vector<int> indVars {0,2};
		std::vector<int> indVars {0, 2};
		// std::vector<int> depVars {4,6};
		std::vector<int> depVars {4, 6};
		std::vector<mirror_t> mirrorTypes {MIRROR_XZ, MIRROR_XZ};

		printf("Using continuation...\n");
		cr3bp_natParamCont(&fam, perOrbit, mirrorTypes, indVars, depVars, 2);
	}else if(contType == PSEUDO_ARC){
		// Turn trajectory object into nodeset; double number of nodes
		tpat_nodeset_cr3bp initGuess(perOrbit, 2*numNodes-1);

		std::vector<int> initDir {1, 0, 0, 0, 0, 0};
		cr3bp_pseudoArcCont(&fam, initGuess, MIRROR_XZ, initDir);
	}
	return fam;
}//====================================================

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
 */
void tpat_family_generator::cr3bp_natParamCont(tpat_family_cr3bp *fam, tpat_traj_cr3bp initialGuess,
	std::vector<mirror_t> mirrorTypes, std::vector<int> indVarIx, std::vector<int> depVarIx, int order){

	tpat_sys_data_cr3bp sys = fam->getSysData();

	if(indVarIx.size() < 2)
		throw tpat_exception("tpat_family_generator::cr3bp_natParamCont: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw tpat_exception("tpat_family_generator::cr3bp_natParamCont: there must be an equal number of ind. vars and mirror types");

	int indVar1 = indVarIx[0];
	int indVar2 = indVarIx[1];
	mirror_t mirrorType = mirrorTypes[0];

	// Initially assume that we're fixing indVar1
	std::vector<int> fixStates;
	fixStates.push_back(indVar1);

	// Get info from the initial guess trajectory
	std::vector<double> IC = initialGuess.getState(0);
	double tof = initialGuess.getTime(-1);

	// Initialize counters and storage containers
	int orbitCount = 0;
	double indVarSlope = NAN;
	double deltaVar1 = 1;
	double deltaVar2 = 1;

	std::vector<tpat_traj_cr3bp> members;
	while(orbitCount < numOrbits){
		tpat_traj_cr3bp perOrbit(&sys);
		try{
			printf("Guess for IC: [%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f] %.4f\n", IC[0], IC[1], IC[2], IC[3],
				IC[4], IC[5], tof);
			printf("Fix States: ");
			for(size_t i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			// Simulate the orbit
			perOrbit = cr3bp_getPeriodic(&sys, IC, tof, numNodes, order, mirrorType, fixStates, tol);
		}catch(tpat_diverge &e){
			break;
		}catch(tpat_linalg_err &e){
			printErr("There was a linear algebra error during family continuation...\n");
			break;
		}

		printf("Orbit %03d converged!\n", ((int)members.size()));

		bool leftFamily = false;

		// Check for large changes in period to detect leaving family
		if(orbitCount > 2){
			double dTOF = perOrbit.getTime(-1) - members[members.size()-1].getTime(-1);
			double percChange = std::abs(dTOF/perOrbit.getTime(-1));
			if(percChange > 0.25){
				leftFamily = true;
				printf("percChange = %.4f\n", percChange);
				printWarn("Period jumped (now = %.5f)! Left the family! Exiting...\n", perOrbit.getTime(-1));
				break;
			}
		}

		if(!leftFamily){
			members.push_back(perOrbit);
			orbitCount++;

			tof = perOrbit.getTime(-1);

			if(orbitCount < numSimple){
				// Use simple continuation; copy the converged IC, step forward in the independent variable
				IC = perOrbit.getState(0);
				IC.at(indVar1) += step_simple;
			}else{

				// Compute the slope for the first time
				if(orbitCount == numSimple){
					deltaVar1 = members[orbitCount-1].getState(0)[indVar1] - members[orbitCount-2].getState(0)[indVar1];
					deltaVar2 = members[orbitCount-1].getState(0)[indVar2] - members[orbitCount-2].getState(0)[indVar2];
					indVarSlope = deltaVar1/deltaVar2;
				}

				// Use least-squares fitting to predict the value for the independent and dependent variables
				std::vector<double> prevStates;
				int first = ((int)members.size()) - curveFitMem < 0 ? 0 : ((int)members.size()) - curveFitMem;

				for(size_t n = first; n < members.size(); n++){
					std::vector<double> ic = members[n].getState(0);
					prevStates.insert(prevStates.end(), ic.begin(), ic.begin()+6);
					prevStates.push_back(members[n].getTime(-1));
					prevStates.push_back(members[n].getJacobi(0));
				}

				// This will hold the input depVars plus the unused independent variable
				std::vector<int> allDepVars;
				std::vector<double> predictedIC;
				if(std::abs(indVarSlope) > slopeThresh){
					mirrorType = mirrorTypes[0];
					// Use continuation in indVar1
					IC.at(indVar1) = perOrbit.getState(0).at(indVar1) + tpat_util::sign(deltaVar1)*step_fitted_1;
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
					IC.at(indVar2) = perOrbit.getState(0).at(indVar2) + tpat_util::sign(deltaVar2)*step_fitted_2;
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
						throw tpat_exception("tpat_family_generator::cr3bp_natParamCont: Cannot update independent variable; index out of range");
				}

				// Update slope
				deltaVar1 = members[orbitCount-1].getState(0)[indVar1] - members[orbitCount-2].getState(0)[indVar1];
				deltaVar2 = members[orbitCount-1].getState(0)[indVar2] - members[orbitCount-2].getState(0)[indVar2];
				indVarSlope = deltaVar1/deltaVar2;
			}

			// Compute eigenvalues
			tpat_matrix mono = perOrbit.getSTM(-1);
			double monoErr = std::abs(1.0 - det(mono));
			if(monoErr > 1e-5)
				printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);
			
			std::vector< std::vector<cdouble> > eigData = eig(mono);
			std::vector<cdouble> eigVals = eigData[0];

			// Add orbit to family
			tpat_family_member_cr3bp child(perOrbit);
			child.setEigVals(eigVals);
			fam->addMember(child);
		}// End of leftFamily?
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
 */
void tpat_family_generator::cr3bp_pseudoArcCont(tpat_family_cr3bp *fam, tpat_nodeset_cr3bp initialGuess,
	mirror_t mirrorType, std::vector<int> initDir){

	// Check inputs (Only applies to Constraint Method 1)
	// if(periodicityIgnoreIx < 0 || periodicityIgnoreIx > 5)
	// 	throw tpat_exception("tpat_family_generator::cr3bp_pseudoArcCont: Periodicity Ignore Index out of range");

	// if(fixToVal_ix < 0 || fixToVal_ix > 7)
	// 	throw tpat_exception("tpat_family_generator::cr3bp_pseudoArcCont: FixToVal Index out of range");

	// TODO - Make these editable?
	double stepSize = 0.001;
	double maxStepSize = 0.01;
	double minStepSize = 1e-7;

	tpat_sys_data_cr3bp sys = fam->getSysData();
	tpat_nodeset_cr3bp familyMember(initialGuess);	// Copy the input initial guess

	printf("Correcting Initial Guess...\n");
	
	// Clear constraints and add new ones
	familyMember.clearConstraints();

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
	tpat_constraint periodicCon(tpat_constraint::MATCH_CUST, familyMember.getNumNodes()-1, periodicConData, 6);

	tpat_constraint extraCon;
	if(fixToVal_ix < 6){
		double extraCon_data[] = {NAN, NAN, NAN, NAN, NAN, NAN};
		extraCon_data[fixToVal_ix] = fixToVal_val;
		extraCon = tpat_constraint(tpat_constraint::STATE, 0, extraCon_data, 6);
	}else if(fixToVal_ix == 6){
		double val = fixToVal_val;
		extraCon = tpat_constraint(tpat_constraint::TOF, 0, &val, 1);
	}else if(fixToVal_ix == 7){
		double val = fixToVal_val;
		extraCon = tpat_constraint(tpat_constraint::JC, 0, &val, 1);
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
		case MIRROR_XZ:
			perpCross_data[1] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case MIRROR_YZ:
            perpCross_data[0] = 0;
			perpCross_data[3] = 0;
			perpCross_data[5] = 0;
            break;
        case MIRROR_XY:
            perpCross_data[2] = 0;
			perpCross_data[3] = 0;
			perpCross_data[4] = 0;
            break;
        case MIRROR_X_AX_H:
        case MIRROR_X_AX_V:
        	perpCross_data[1] = 0;
			perpCross_data[2] = 0;
			perpCross_data[3] = 0;
            break;
        default:
            throw tpat_exception("Mirror type either not defined or not implemented");
	}
	// Just for kicks: try constraining perpendicular crossings and periodicity
	tpat_constraint perpCross1_Con(tpat_constraint::STATE, 0, perpCross_data, 6);
	tpat_constraint perpCross2_Con(tpat_constraint::STATE, familyMember.getNumNodes()-1, perpCross_data, 6);

	familyMember.addConstraint(perpCross1_Con);
	familyMember.addConstraint(perpCross2_Con);

	std::vector<tpat_constraint> constraints {perpCross1_Con, perpCross2_Con};

	// Correct the nodeset to retrieve a free-variable vector for a family member
	tpat_correction_engine corrector;
	corrector.setVarTime(true);			// Variable time MUST be enabled for PAC
	corrector.setEqualArcTime(true);	// MUST use equal arc time to get propper # of constraints
	corrector.setTol(tol);
	corrector.setIgnoreCrash(true);		// Ignore crashes into primary
	iterationData familyItData;
	try{
		familyItData = corrector.correct(&familyMember);
	}catch(tpat_diverge &e){
		printErr("tpat_family_generator::cr3bp_pseudoArcCont: Could not converge initial guess!\n");
	}catch(tpat_linalg_err &e){
		printErr("tpat_family_generator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
	}

	printf("Applying continuation to compute family...\n");
	
	// Initialize counters and storage containers
	int orbitCount = 0;
	tpat_matrix convergedFreeVarVec(familyItData.totalFree, 1, familyItData.X);
	tpat_matrix prevN(familyItData.totalFree, 1);
	std::vector<tpat_traj_cr3bp> members;

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
		tpat_matrix DF(familyItData.totalFree-1, familyItData.totalFree, DF_data);
		tpat_matrix N = null_qr(DF);

		printf("DF has dimensions %d x %d\n", DF.getRows(), DF.getCols());
		// Check to make sure the IS a nullspace
		if(N.getRows() == 1){
			printErr("tpat_family_generator::cr3bp_pseudoArcCont: Nullspace is zero-dimensional; cannot proceed...\n");
			return;
		}		

		// // For debugging, save nullspace vectors to file
		// char filename[16];
		// sprintf(filename, "N%02d.csv", orbitCount);
		// N.toCSV(filename);

		/**
		 *	Choose the nullspace vector that is closest to the previous one (which converged)
		 */
		printf("Choosing Nullspace Vector (%dD, %d elements)\n", N.getCols(), N.getRows());
		if(orbitCount == 0){
			if(N.getCols() > 1){
				printErr("tpat_family_generator::cr3bp_pseudoArcCont: Nullspace is multidimensional on first iteration; unsure how to proceed...\n");
				return;
			}

			bool sameDir = true;
			for(size_t i = 0; i < initDir.size(); i++){
				// Find a non-zero element
				if(initDir[i] != 0 && i < (size_t)(N.getCols())){
					// If signs are different, assume direction is different
					if(N.at(i)*initDir[i] < 0){
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
			for(int i = 0; i < N.getCols(); i++){
				// Compute angle from dot product
				tpat_matrix col_i = N.getCols() > 1 ? N.getCol(i) : N;
				tpat_matrix dotProd = trans(prevN)*col_i / (norm(prevN)*norm(N));
				double angle = std::acos(dotProd.at(0));
				int sign = 1;

				// Flip the sign if the angle is greater than 90 and change the value to 180 - angle
				printf("dot product = %.4f\n", dotProd.at(0));
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
			tpat_matrix temp = N.getCols() > 1 ? N.getCol(best_ix) : N;
			N = best_sign*temp;	// Apply sign change, if needed
		}

		prevN = N;	// Update memory
		printf("Chose N with first elements = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...]\n",
			N.at(0), N.at(1), N.at(2), N.at(3), N.at(4), N.at(5));

		tpat_nodeset_cr3bp newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData, constraints);

		/*
		 *	Apply multiple shooting to converge the new guess to be a member of the family
		 */
		bool killLoop = false;
		try{
			while(stepSize >= minStepSize){
				try{
					familyItData = corrector.correct(&newMember);

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
				}catch(tpat_diverge &e){
					if(stepSize > minStepSize){
						printWarn("tpat_family_generator::cr3bp_pseudoArcCont: Corrector diverged... trying smaller step size\n");

						// Decrease step size and try again
						stepSize = stepSize/2 > minStepSize ? stepSize/2 : minStepSize;
						printColor(MAGENTA, "Decreased Step Size to %.4e (min %.4e)!\n", stepSize, minStepSize);

						// Re-Create the initial guess using the new step size
						newMember = cr3bp_getNextPACGuess(convergedFreeVarVec, N, stepSize, familyItData, constraints);
					}else{
						printErr("tpat_family_generator::cr3bp_pseudoArcCont: Could not converge new family member!\n");
						killLoop = true;
						break;
					}
				}
			}
		}catch(tpat_linalg_err &e){
			printErr("tpat_family_generator::cr3bp_pseudoArcCont: There was a linear algebra error...\n");
			killLoop = true;
		}

		if(killLoop)
			break;

		printf("Orbit %03d converged!\n", ((int)members.size()));

		// Save new converged family vector
		convergedFreeVarVec = tpat_matrix(familyItData.totalFree, 1, familyItData.X);

		// Convert converged nodeset to an orbit to save; TODO - could be improved to be much faster!
		tpat_nodeset_cr3bp perNodes = corrector.getCR3BP_Output();

		tpat_traj_cr3bp perOrbit = tpat_traj_cr3bp::fromNodeset(perNodes);

		members.push_back(perOrbit);
		orbitCount++;

		// Compute eigenvalues
		tpat_matrix mono = perOrbit.getSTM(-1);
		double monoErr = std::abs(1.0 - det(mono));
		// if(monoErr > 1e-5)
		printColor(BOLDRED, "Monodromy Matrix error = %.4e; This will affect eigenvalue accuracy!\n", monoErr);

		std::vector< std::vector<cdouble> > eigData = eig(mono);
		std::vector<cdouble> eigVals = eigData[0];

		// Add orbit to family
		tpat_family_member_cr3bp child(perOrbit);
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
 *	@param familyItData an iterationData object containing corrections information about the
 *	previous (nearest) converged family member
 *	@param cons a vector of constraints to place on the nodeset
 */
tpat_nodeset_cr3bp tpat_family_generator::cr3bp_getNextPACGuess(tpat_matrix convergedFreeVarVec,
	tpat_matrix N, double stepSize, iterationData familyItData, std::vector<tpat_constraint> cons){

	/**
	 *	Step forwards away from previously converged solution
	 */
	tpat_matrix newFreeVarVec = convergedFreeVarVec + stepSize*N;
	double *X = newFreeVarVec.getDataPtr();

	// Convert into a new nodeset (TODO: Make this more flexible by putting conversion code in a model?)
	tpat_sys_data_cr3bp *sys = static_cast<tpat_sys_data_cr3bp *>(familyItData.sysData);
	tpat_nodeset_cr3bp newMember(sys);

	for(int n = 0; n < familyItData.numNodes; n++){
		// tof stored in the element after all the nodes; NAN for last node
		double tof = n < familyItData.numNodes - 1 ? X[6*familyItData.numNodes]/(familyItData.numNodes - 1) : NAN;
		tpat_node node(X+6*n, tof);
		newMember.appendNode(node);
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
	pacCon_data.insert(pacCon_data.end(), N.getDataPtr(), N.getDataPtr()+N.getRows());
	// Append the step size
	pacCon_data.insert(pacCon_data.end(), stepSize);
	// Create the actual constraint
	tpat_constraint pacCon(tpat_constraint::PSEUDOARC, familyItData.numNodes-1, pacCon_data);
	newMember.addConstraint(pacCon);

	// Outputs for debugging and sanity checks
	printColor(RED, "||delta-X|| = %.4e\n", norm(newFreeVarVec - convergedFreeVarVec));
	printf("New IC = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, ...] tof = %.4f\n",
		X[0], X[1], X[2], X[3], X[4], X[5], X[newFreeVarVec.getRows()-1]);

	return newMember;
}

/**
 *	@brief Reset all parameters to their default values
 */
void tpat_family_generator::reset(){
	numOrbits = 500;
	numSimple = 3;
	step_simple = 0.0005;
	step_fitted_1 = 0.005;
	step_fitted_2 = 0.005;
	curveFitMem = 5;
	numNodes = 3;
	slopeThresh = 1;
	tol = 1e-12;
}//======================================

//