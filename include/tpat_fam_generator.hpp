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

#ifndef H_FAMILY_GENERATOR
#define H_FAMILY_GENERATOR

#include "tpat.hpp"
 
#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_eigen_defs.hpp"
 
#include <vector>

// forward declarations
class TPAT_Fam_CR3BP;
class TPAT_Nodeset_CR3BP;
class TPAT_Sys_Data_CR3BP;
class TPAT_Traj_CR3BP;

/**
 *	@brief Type of continuation to use when generating a family
 */
enum class TPAT_Continuation_Tp{
	NAT_PARAM,	//!< Use natural parameter continuation
	PSEUDO_ARC	//!> Use pseudo arclength continuation
};
		
/**
 *	@brief An object to generate families of orbits
 */
class TPAT_Fam_Generator : public TPAT{
	public:
		// *structors
		TPAT_Fam_Generator();
		TPAT_Fam_Generator(const TPAT_Fam_Generator&);
		~TPAT_Fam_Generator();

		// operators
		TPAT_Fam_Generator& operator =(const TPAT_Fam_Generator&);

		// Set and Get
		void setContType(TPAT_Continuation_Tp);
		void setMaxStepSize(double);
		void setMinStepSize(double);
		void setNumNodes(int);
		void setNumOrbits(int);
		void setSlopeThresh(double);
		void setStep_simple(double);
		void setStep_fitted_1(double);
		void setStep_fitted_2(double);
		void setTol(double);
		
		// Operations & Utility
		TPAT_Fam_CR3BP cr3bp_generateAxial(const char*, double);
		TPAT_Fam_CR3BP cr3bp_generateButterfly(TPAT_Sys_Data_CR3BP*, int);
		TPAT_Fam_CR3BP cr3bp_generateDRO(TPAT_Sys_Data_CR3BP*);
		TPAT_Fam_CR3BP cr3bp_generateHalo(const char*, double);
		TPAT_Fam_CR3BP cr3bp_generateLPO(TPAT_Sys_Data_CR3BP*);
		TPAT_Fam_CR3BP cr3bp_generateLyap(TPAT_Sys_Data_CR3BP, int, double);
		TPAT_Fam_CR3BP cr3bp_generateRes(TPAT_Sys_Data_CR3BP*, int, int);
		TPAT_Fam_CR3BP cr3bp_generateVertical(const char*, double);
		
		void reset();

	private:
		TPAT_Continuation_Tp contType = TPAT_Continuation_Tp::NAT_PARAM;	//!< Type of continuation to use
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

		void copyMe(const TPAT_Fam_Generator&);
		void cr3bp_natParamCont(TPAT_Fam_CR3BP*, TPAT_Traj_CR3BP, std::vector<TPAT_Mirror_Tp>, std::vector<int>, std::vector<int>, int);
		void cr3bp_pseudoArcCont(TPAT_Fam_CR3BP*, TPAT_Nodeset_CR3BP, TPAT_Mirror_Tp, std::vector<int>);
		TPAT_Nodeset_CR3BP cr3bp_getNextPACGuess(Eigen::VectorXd, Eigen::VectorXd, double, TPAT_MultShoot_Data, std::vector<TPAT_Constraint>);
};

#endif