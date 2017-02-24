/**
 *  @file tpat_linMotion_Engine.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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

#include "Core.hpp"
#include "Engine.hpp"

namespace astrohelion{

// Forward delcarations
class Traj_cr3bp;
class SysData_cr3bp;

/**
 *	@ingroup engine
 *	@brief type of linear motion
 *
 *	This type tells the engine what kind of linearization to produce
 */
enum class LinMotion_tp {
	NONE,	//!< No motion specified
	HYP,	//!< Hyperbolic motion, applies only to linearizations near collinear points
	ELLIP,	//!< Elliptical motion, applies only to linearizations near collinear points
	SPO,	//!< Short-Period Orbit, applies only to Case I linearizations near triangular points
	LPO,	//!< Long-Period Orbit, applies only to Case I linearizations near triangular points
	MPO,	//!< Mixed-Period Orbit, applies only to Case I linearizations near triangular points
	CONVERGE,	//!< Convergent motion, applies only to Case III linearizations near triangular points
	DIVERGE		//!< Divergent motion, applies only to Case III linearizations near triangular points
};	

/**
 *	@brief An engine that will generate a trajectory from the linearized CR3BP EOMs
 */
class LinMotionEngine : public Core, public Engine{
	public:
		// *structors
		LinMotionEngine();

		// Operators

		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		double getMPORatio() const;
		int getNumRevs() const;
		double getTimeStep() const;
		double getTol() const;

		void setMPORatio(double);
		void setNumRevs(int);
		void setTimeStep(double);
		void setTol(double);
		const char* getTypeStr(LinMotion_tp) const;
		//@}

		/**
		 *  @name Orbit Generation
		 *  @{
		 */
		Traj_cr3bp getCR3BPLinear(int, double[3], SysData_cr3bp*);
		Traj_cr3bp getCR3BPLinear(int, double[3], LinMotion_tp, SysData_cr3bp*);
		Traj_cr3bp getCR3BPLinear(int, double[3], double, double, LinMotion_tp, SysData_cr3bp*);
		Traj_cr3bp getCR3BPLiss(int, double, bool, double, double, double, SysData_cr3bp*);
		//@}
		
	private:
		/** @brief step size between points on linear motion trajectory */
		double t_step = 0.001;

		/** @brief Number of rotations to propagate */
		double rots = 1;

		/** @brief tolerance for numerical methods, like locating Lagrange points */
		double tol = 1e-14;

		/** Ratio between SPO and LPO behavior when constructing an MPO */
		double nu = 1;

		void cleanEngine();
		void reset();
};

}// END of Astrohelion namespace