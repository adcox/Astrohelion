/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#ifndef H_LIN_MOTION
#define H_LIN_MOTION

#include "tpat.hpp"
 
// Forward delcarations
class TPAT_Traj_CR3BP;
class TPAT_Sys_Data_CR3BP;

/**
 *	@brief type of linear motion
 *
 *	This type tells the engine what kind of linearization to produce
 */
enum class TPAT_LinMotion_Tp {
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
class TPAT_LinMotion_Engine : public TPAT{
	public:
		// *structors
		TPAT_LinMotion_Engine();

		// Operators

		// Set and Get
		double getMPORatio() const;
		int getNumRevs() const;
		double getTimeStep() const;
		double getTol() const;

		void setMPORatio(double);
		void setNumRevs(int);
		void setTimeStep(double);
		void setTol(double);
		const char* getTypeStr(TPAT_LinMotion_Tp) const;

		// Misc
		TPAT_Traj_CR3BP getCR3BPLinear(int, double[3], TPAT_Sys_Data_CR3BP*);
		TPAT_Traj_CR3BP getCR3BPLinear(int, double[3], TPAT_LinMotion_Tp, TPAT_Sys_Data_CR3BP*);
		TPAT_Traj_CR3BP getCR3BPLinear(int, double[3], double, double, TPAT_LinMotion_Tp, TPAT_Sys_Data_CR3BP*);
		TPAT_Traj_CR3BP getCR3BPLiss(int, double, bool, double, double, double, TPAT_Sys_Data_CR3BP*);
	private:
		/** @brief step size between points on linear motion trajectory */
		double t_step = 0.001;

		/** @brief Number of rotations to propagate */
		double rots = 1;

		/** @brief tolerance for numerical methods, like locating Lagrange points */
		double tol = 1e-14;

		/** Ratio between SPO and LPO behavior when constructing an MPO */
		double nu = 1;
};

#endif