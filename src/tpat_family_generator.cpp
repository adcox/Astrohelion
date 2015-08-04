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

#include "tpat_correction_engine.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_constraint.hpp"
#include "tpat_linear_motion_engine.hpp"
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
		printErr("Could not locate any bifurcations in the Lyapunov family; extiting...");
		return axialFam;
	}

	if(bifs.size() != 3)
		printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!");

	std::vector<double> IC = lyapFam.getMember(bifs[1]).getIC();
	double period = lyapFam.getMember(bifs[1]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {5};	// force z-dot to be non-zero
	IC[5] += initStepSize;

	tpat_traj_cr3bp firstAxial = cr3bp_getPeriodic(axialFam.getSysData(), IC, period,
		numNodes, MIRROR_X_AX_H, fixStates);

	std::vector<int> indVars {5,4};	// begin stepping in z-dot, optionally use y-dot
	std::vector<int> depVars {0,6};	// Predict x and period with least squares
	std::vector<mirror_t> mirrorTypes {MIRROR_X_AX_H, MIRROR_X_AX_V};

	cr3bp_continueFamily(&axialFam, firstAxial, mirrorTypes, indVars, depVars);

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
		printErr("Could not locate any bifurcations in the Lyapunov family; extiting...");
		return haloFam;
	}

	if(bifs.size() != 3)
		printWarn("The # of bifurcations in the Lyap family != 3... something may be wrong!");

	std::vector<double> IC = lyapFam.getMember(bifs[0]).getIC();
	double period = lyapFam.getMember(bifs[0]).getTOF();
	int numNodes = 3;
	std::vector<int> fixStates {2};	// force z to be out of plane
	IC[2] += initStepSize;

	tpat_traj_cr3bp firstHalo = cr3bp_getPeriodic(haloFam.getSysData(), IC, period,
		numNodes, MIRROR_XZ, fixStates);

	std::vector<int> indVars {2,0};	// begin stepping in z, optionally using x
	std::vector<int> depVars {4};	// Predict y-dot with least-squares
	std::vector<mirror_t> mirrorTypes {MIRROR_XZ, MIRROR_XZ};

	cr3bp_continueFamily(&haloFam, firstHalo, mirrorTypes, indVars, depVars);

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

	double LPt_data[] = {0,0,0};
	cr3bp_getEquilibPt(sysData, LPt, 1e-14, LPt_data);

	// Begin solving - get linear approximation at ICs
	double r0[] = {x0, 0, 0};
	tpat_linear_motion_engine linEngine;
	tpat_traj_cr3bp linTraj = linEngine.getCR3BPLinear(LPt, r0,
		tpat_linear_motion_engine::ELLIP, sysData);

	// Initialize variables and containers for data
	tpat_family_cr3bp fam(sysData);

	// Lyapunov-specific settings
	fam.setSortType(tpat_family_cr3bp::SORT_X);
	std::vector<int> indVars;
	indVars.push_back(0);	// We're going to fix the x-coordinate in the corrector to keep it from slipping
	indVars.push_back(4);	// Optionally, allow y-dot to be an independent variable if x is changing too quickly
	std::vector<int> depVars {4}; // Predict y-dot with least squares in the algorithm
	std::vector<mirror_t> mirrorTypes {MIRROR_XZ, MIRROR_XZ};

	cr3bp_continueFamily(&fam, linTraj, mirrorTypes, indVars, depVars);

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
 */
void tpat_family_generator::cr3bp_continueFamily(tpat_family_cr3bp *fam,
	tpat_traj_cr3bp initialGuess, std::vector<mirror_t> mirrorTypes, std::vector<int> indVarIx, std::vector<int> depVarIx){

	tpat_sys_data_cr3bp sys = fam->getSysData();

	if(indVarIx.size() < 2)
		throw tpat_exception("tpat_family_generator::cr3bp_continueFamily: Must specify two independent variables");

	if(mirrorTypes.size() != indVarIx.size())
		throw tpat_exception("tpat_family_generator::cr3bp_continueFamily: there must be an equal number of ind. vars and mirror types");

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
	double indVarSlope = 0;
	double deltaVar1 = 1;
	double deltaVar2 = 1;

	std::vector<tpat_traj_cr3bp> members;
	while(orbitCount < numOrbits){
		tpat_traj_cr3bp perOrbit;
		try{
			printf("IC: [%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f] %.4f\n", IC[0], IC[1], IC[2], IC[3],
				IC[4], IC[5], tof);
			printf("Fix States: ");
			for(size_t i = 0; i < fixStates.size(); i++){ printf("%d, ", fixStates[i]); }
			printf("\n");
			printf("Slope = %.3f\n", indVarSlope);

			// Simulate the orbit
			perOrbit = cr3bp_getPeriodic(sys, IC, tof, numNodes, mirrorType, fixStates);
		}catch(tpat_diverge &e){
			break;
		}catch(tpat_linalg_err &e){
			printErr("There was a linear algebra error during family continuation...");
			break;
		}

		printf("Orbit %03d converged!\n", ((int)members.size()));

		bool leftFamily = false;

		// Check for large changes in period to detect leaving family
		if(orbitCount > 2){
			if(perOrbit.getTime(-1) > 1.75*members[members.size()-1].getTime(-1) ||
				perOrbit.getTime(-1) < 0.4*members[members.size()-1].getTime(-1)){
				
				leftFamily = true;
				printWarn("Period jumped! Left the family! Exiting...");
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
						throw tpat_exception("tpat_family_generator::cr3bp_continueFamily: Cannot update independent variable; index out of range");
				}

				// Update slope
				deltaVar1 = members[orbitCount-1].getState(0)[indVar1] - members[orbitCount-2].getState(0)[indVar1];
				deltaVar2 = members[orbitCount-1].getState(0)[indVar2] - members[orbitCount-2].getState(0)[indVar2];
				indVarSlope = deltaVar1/deltaVar2;
			}

			// Compute eigenvalues
			tpat_matrix mono = perOrbit.getSTM(-1);
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
}//======================================

//