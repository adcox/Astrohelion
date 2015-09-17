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

#include "tpat_calculations.hpp"

#include <vector>

// forward declarations
class tpat_family_cr3bp;
class tpat_sys_data_cr3bp;
class tpat_traj_cr3bp;

/**
 *	@brief An object to generate families of orbits
 */
class tpat_family_generator{
	public:
		// *structors
		tpat_family_generator();
		tpat_family_generator(const tpat_family_generator&);
		~tpat_family_generator();

		// operators
		tpat_family_generator& operator =(const tpat_family_generator&);

		// Set and Get
		void setNumNodes(int);
		void setNumOrbits(int);
		void setSlopeThresh(double);
		void setStep_simple(double);
		void setStep_fitted_1(double);
		void setStep_fitted_2(double);

		// Operations & Utility
		tpat_family_cr3bp cr3bp_generateAxial(const char*, double);
		tpat_family_cr3bp cr3bp_generateButterfly(tpat_sys_data_cr3bp, int, double);
		tpat_family_cr3bp cr3bp_generateHalo(const char*, double);
		tpat_family_cr3bp cr3bp_generateLyap(tpat_sys_data_cr3bp, int, double);
		void reset();

	private:
		int numOrbits = 500;			//!< Maximum number of family members to generate
		int numSimple = 3;				//!< Number of simply-continued family members
		double step_simple = 0.0005;	//!< Step size in the independent variable when using simple continuation
		double step_fitted_1 = 0.005;	//!< Step size in first ind. var. when using advanced continuation
		double step_fitted_2 = 0.005;	//!< Step size in second ind. var. when using advanced continuation
		int curveFitMem = 5;			//!< Number of points to use with Least-Squares algorithm
		int numNodes = 3;				//!< Number of nodes to use when correcting HALF a periodic orbit
		double slopeThresh = 1;			//!< Minimum slope for stepping in indVar1; else step in indVar2
		void copyMe(const tpat_family_generator&);
		void cr3bp_continueFamily(tpat_family_cr3bp*, tpat_traj_cr3bp, std::vector<mirror_t>, std::vector<int>, std::vector<int>);
};

#endif