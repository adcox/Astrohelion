/**
 *  @file FamGenerator.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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
#pragma once

#include "Core.hpp"
#include "Engine.hpp"

#include "Calculations.hpp"
#include "CorrectionEngine.hpp"
#include "EigenDefs.hpp"
 
#include <vector>


namespace astrohelion{

// forward declarations
class Fam_cr3bp;
class Nodeset_cr3bp;
class SysData_cr3bp;
class Traj_cr3bp;

/**
 *	@brief Type of continuation to use when generating a family
 */
enum class Continuation_tp{
	NAT_PARAM,	//!< Use natural parameter continuation
	PSEUDO_ARC	//!> Use pseudo arclength continuation
};
		
/**
 *	@ingroup engine
 *	@brief An object to generate families of orbits
 */
class FamGenerator : public Core, public Engine{
	public:
		/**
		 *  @name *structors
		 *  @{
		 */
		FamGenerator();
		FamGenerator(const FamGenerator&);
		~FamGenerator();
		//@}

		// operators
		FamGenerator& operator =(const FamGenerator&);

		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		void setContType(Continuation_tp);
		void setMaxStepSize(double);
		void setMinStepSize(double);
		void setNumNodes(int);
		void setNumOrbits(int);
		void setSlopeThresh(double);
		void setStep_simple(double);
		void setStep_fitted_1(double);
		void setStep_fitted_2(double);
		void setTol(double);
		//@}

		/**
		 *  @name Family Generation
		 *  @{
		 */
		void cr3bp_generateAxial(const char*, double, Fam_cr3bp*);
		void cr3bp_generateButterfly(int, Fam_cr3bp*);
		void cr3bp_generateDPO(Fam_cr3bp*);
		void cr3bp_generateDRO(Fam_cr3bp*);
		void cr3bp_generateHalo(const char*, double, Fam_cr3bp*);
		void cr3bp_generateLPO(Fam_cr3bp*);
		void cr3bp_generateLyap(int, double, Fam_cr3bp*);
		void cr3bp_generateRes(int, int, Fam_cr3bp*);
		void cr3bp_generateVertical(const char*, double, Fam_cr3bp*);

		void cr3bp_pacFromTraj(Traj_cr3bp, Mirror_tp, std::vector<int>, Fam_cr3bp*);
		void cr3bp_pacFromNodeset(Nodeset_cr3bp, Mirror_tp, std::vector<int>, Fam_cr3bp*);
		//@}
		
		void reset();

	private:
		Continuation_tp contType = Continuation_tp::NAT_PARAM;	//!< Type of continuation to use
		double tol = 1e-12;				//!< Tolerance for corrections

		// Settings for Natural Parameter Continuation
		int numOrbits = 500;			//!< Maximum number of family members to generate
		int numSimple = 3;				//!< Number of simply-continued family members
		double step_simple = 0.0005;	//!< Step size in the independent variable when using simple continuation
		double step_fitted_1 = 0.005;	//!< Step size in first ind. var. when using advanced continuation
		double step_fitted_2 = 0.005;	//!< Step size in second ind. var. when using advanced continuation
		double minStepSize = 1e-6;		//!< Minimum allowable step size
		double maxStepSize = 0.05;		//!< Maximum allowable step size
		int curveFitMem = 5;			//!< Number of points to use with Least-Squares algorithm
		int numNodes = 3;				//!< Number of nodes to use when correcting HALF a periodic orbit
		double slopeThresh = 1;			//!< Minimum slope for stepping in indVar1; else step in indVar2

		void copyMe(const FamGenerator&);
		void cleanEngine();
		void cr3bp_natParamCont(Fam_cr3bp*, Traj_cr3bp, std::vector<Mirror_tp>, std::vector<int>, std::vector<int>, int order = 1);
		void cr3bp_pseudoArcCont(Fam_cr3bp*, Nodeset_cr3bp, Mirror_tp, std::vector<int>);
		Nodeset_cr3bp cr3bp_getNextPACGuess(Eigen::VectorXd, Eigen::VectorXd, double, MultShootData);
};

}