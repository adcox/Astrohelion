/**
 * \file ManifoldEngine.hpp
 * \brief Computes manifolds
 * 
 * \author Andrew Cox
 * \version April 19, 2017
 * \copyright GNU GPL v3.0
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
#pragma once

#include <vector>

#include "Core.hpp"
#include "Engine.hpp"

#include "EigenDefs.hpp"

namespace astrohelion{

// Forward Declarations
class Arcset_cr3bp;
class SysData_cr3bp;

/**
 *	\brief Describes the type of manifold, both stability and direction
 *	\details The value of the type species the number, stability, and direction of the manifold.
 *	For example, all unstable manifold options are positive values, all stable manifold options are
 *	negative (think reverse time propagation), and the MAN_ALL option is zero. Similarly, selecting 
 *	all unstable manifolds has a value of 1 (all stable is a value of -1)
 */
enum class Manifold_tp : int{
	MAN_ALL = 0,		//!< Compute all types of manifolds, arriving/departing from all directions
	MAN_U = 1,			//!< Unstable, compute both arcs arriving from +x and -x directions
	MAN_U_RIGHT = 2,	//!< Unstable, departing towards +x direction
	MAN_U_LEFT = 3,		//!< Unstable, departing towards -x direction
	MAN_S = -1,			//!< Stable, compute both arcs arriving from +x and -x directions
	MAN_S_RIGHT = -2,	//!< Stable, arriving from +x direction
	MAN_S_LEFT = -3		//!< Stable, arriving from -x direction
};

enum class Manifold_StepOff_tp{
	STEP_MATCH_JC,		//!< Step along eigenvector as much as possible but adjust to match the Jacobi constant value of the origin
	STEP_VEC_NORMPOS,	//!< Step along eigenvector by prescribed distance; eigenvector is normalized by magnitude of position components
	STEP_VEC_NORMFULL	//!< Step along eigenvector by prescribed distance; eigenvector is normalized by magnitude of all (position and velocity) components
};

/**
 *  \ingroup engine
 *  \brief Compute manifolds
 * 
 * 	\author Andrew Cox
 * 	\version April 19, 2017
 * 	\copyright GNU GPL v3.0
 */
class ManifoldEngine : public Core, public Engine{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	ManifoldEngine();
	ManifoldEngine(const ManifoldEngine&);
	//\}

	/**
	 * \name Set and Get Functions
	 * \{
	 */
	double getStepOffDist() const;
	void setStepOffDist(double);
	//\}

	/**
	 *  \name Manifold Generation Algorithms
	 *  \{
	 */
	
	std::vector<Arcset_cr3bp> computeSetFromPeriodic(Manifold_tp, const Arcset_cr3bp*, unsigned int, double, Manifold_StepOff_tp stepType = Manifold_StepOff_tp::STEP_MATCH_JC);
	std::vector<Arcset_cr3bp> computeSingleFromPeriodic(Manifold_tp, const Arcset_cr3bp*, double, double, Manifold_StepOff_tp stepType = Manifold_StepOff_tp::STEP_MATCH_JC);
	std::vector<Arcset_cr3bp> manifoldsFromPOPoint(Manifold_tp, std::vector<double>, MatrixXRd, std::vector<cdouble>, MatrixXRd, double, const SysData_cr3bp*, Manifold_StepOff_tp stepType = Manifold_StepOff_tp::STEP_MATCH_JC);
	
	// Add a function for manifolds from map fixed points; use Wayne's method?

	//\}

	/**
	 * \name Utility Functions
	 * \{
	 */

	MatrixXRd eigVecValFromPeriodic(Manifold_tp, const Arcset_cr3bp*, std::vector<cdouble> *eigVal_final = nullptr);
	void reset();
	//\}
private:

	double tol_eigVal = 1e-5;	//!< Tolerance with which eigenvalues are evaluated to determine if they are on the unit circle
	double stepOffDist = 20;	//!< Dimensional distance to step along eigenvector, km

	void copyMe(const ManifoldEngine&);
	void cleanEngine(){}
	static bool compareCDouble(cdouble, cdouble);
};// End of ManifoldEngine class

}// END of astrohelion namespace