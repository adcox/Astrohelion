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
#include "tpat_cr3bp_family.hpp"
#include "tpat_cr3bp_family_member.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_cr3bp_traj.hpp"
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
 */
tpat_family_generator& tpat_family_generator::operator =(const tpat_family_generator &f){
	copyMe(f);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Methods
//-----------------------------------------------------

/** 
 *	@brief Set the step size we take when performing simple continuation
 *
 *	The default value is 0.0005
 *	@param d the step size, non-dimensional units
 */
void tpat_family_generator::setStep_simple(double d){ step_simple = d; }

/** 
 *	@brief Set the step size we take when performing simple continuation
 *
 *	The default value is 0.005
 *	@param d the step size, non-dimensional units
 */
void tpat_family_generator::setStep_fitted(double d){ step_fitted = d; }

/**
 *	@brief Set the number of nodes used for corrections processes
 *
 *	The default value is 3
 *	@param n the number of nodes
 */
void tpat_family_generator::setNumNodes(int n){ numNodes = n; }

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
	step_fitted = f.step_fitted;
	curveFitMem = f.curveFitMem;
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
 *	@param r0 the initial displacement from the Lagrange point along the x-axis.
 *
 *	@return a family of orbits
 */
tpat_cr3bp_family tpat_family_generator::cr3bp_generateLyap(tpat_cr3bp_sys_data sysData, int LPt, double x0){
	if(LPt < 1 || LPt > 3)
		throw tpat_exception("tpat_family_generator::cr3bp_generateLyap: Invalide LPt number");

	double LPt_data[] = {0,0,0};
	cr3bp_getEquilibPt(sysData, LPt, 1e-14, LPt_data);

	// Begin solving - get linear approximation at ICs
	double r0[] = {x0, 0, 0};
	tpat_linear_motion_engine linEngine;
	tpat_cr3bp_traj linTraj = linEngine.getCR3BPLinear(LPt, r0,
		tpat_linear_motion_engine::ELLIP, sysData);

	// Initialize variables and containers for data
	tpat_cr3bp_family fam(sysData);

	// Lyapunov-specific settings
	fam.setSortType(tpat_cr3bp_family::SORT_X);
	std::vector<int> fixStates;
	fixStates.push_back(0);	// We're going to fix the x-coordinate in the corrector to keep it from slipping
	std::vector<int> depVars;
	depVars.push_back(4);	// Predict y-dot with least squares in the algorithm

	cr3bp_continueFamily(&fam, sysData, linTraj, MIRROR_XZ, fixStates, depVars);

	return fam;
}//====================================================

/**
 *	@brief Continue a family of orbits in the CR3BP
 *	
 * 	@param fam a pointer to a family object to store family members in
 *	@param sys the system data for the system we're integrating in
 *	@param initialGuess a trajectory that is a good initial guess for the "first" member of the family
 *	@param mirrorType describes how the family mirrors in the rotating reference frame
 *	@param fixStates a list of state indices telling the algorithm which states should be
 *	held constant in the initial state. It is usually useful to fix the independent variable
 *	and let the others vary
 *	@param depVarIx a list of state indices telling the algorithm which states should be predicted
 *	by a 2nd-order least squares approximation. If left empty, the continuation scheme will use
 *	simple techniques that don't perform very well.
 */
void tpat_family_generator::cr3bp_continueFamily(tpat_cr3bp_family *fam, tpat_cr3bp_sys_data sys,
	tpat_cr3bp_traj initialGuess, mirror_t mirrorType, std::vector<int> fixStates,
	std::vector<int> depVarIx){

	// the sort type enum is set up to map to the index of the independent variable
	int indVarIx = static_cast<int>(fam->getSortType());

	// Get info from the initial guess trajectory
	std::vector<double> IC = initialGuess.getState(0);
	double tof = initialGuess.getTime(-1);

	// Initialize counters and storage containers
	int orbitCount = 0;
	std::vector<tpat_cr3bp_traj> members;
	while(orbitCount < numOrbits){
		tpat_cr3bp_traj perOrbit;
		try{
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
			if(perOrbit.getTime(-1) > 1.75*members[members.size()-1].getTime(-1)){
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
				IC.at(indVarIx) += step_simple;
			}else{
				// Use least-squares fitting to predict the value for the independent and dependent variables
				std::vector<double> prevStates;
				int first = ((int)members.size()) - curveFitMem < 0 ? 0 : ((int)members.size()) - curveFitMem;

				for(size_t n = first; n < members.size(); n++){
					std::vector<double> ic = members[n].getState(0);
					prevStates.insert(prevStates.end(), ic.begin(), ic.begin()+6);
					prevStates.push_back(members[n].getTime(-1));
					prevStates.push_back(members[n].getJC(0));
				}

				// Update independent variable
				IC.at(indVarIx) = perOrbit.getState(0).at(indVarIx) + step_fitted;

				// Use least squares to predict new values for the dependent variables
				std::vector<double> predictedIC = familyCont_LS(indVarIx, IC.at(indVarIx),
					depVarIx, prevStates);

				// Update IC with predicted variables
				for(size_t n = 0; n < depVarIx.size(); n++){
					int ix = depVarIx[n];
					if(ix < 6)
						IC[ix] = predictedIC[ix];
					else if(ix == 6)
						tof = predictedIC[ix];
					else
						throw tpat_exception("tpat_family_generator::cr3bp_continueFamily: Cannot update independent variable; index out of range");
				}
			}

			// Compute eigenvalues
			tpat_matrix mono = perOrbit.getSTM(-1);
			std::vector<cdouble> eigVals = eig(mono);

			// Add orbit to family
			tpat_cr3bp_family_member child(perOrbit);
			child.setEigVals(eigVals);
			fam->addMember(child);
		}// End of leftFamily?
	}// end of while loop
}//==================================================



//