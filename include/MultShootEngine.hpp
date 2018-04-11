/**
 *  @file MultShootEngine.hpp
 *	@brief Performs corrections on a nodeset subject to a set of constraints
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <vector>

#include "Core.hpp"
#include "Engine.hpp"

#include "Arcset.hpp"
#include "Common.hpp"
#include "Constraint.hpp"
#include "EigenDefs.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "SysData.hpp"



namespace astrohelion{

// Forward declarations
class Constraint;
class MultShootData;
class Arcset_bc4bp;
class Arcset_cr3bp;
class SimEngine;
class SysData_bc4bp;

/**
 *	\ingroup engine
 *	@brief An engine object to perform corrections, such as multiple shooting.
 *
 *	### Multiple Shooting Algorithm
 *	The multiple shooting algorithm is initiated by calling the multShoot() 
 *	function. Typical use calls the version which requires a nodeset pointer 
 *	as the only input.
 *	
 *	@author Andrew Cox
 *	@version Mar 26, 2018
 *	@copyright GNU GPL v3.0
 */
class MultShootEngine : public Core, public Engine{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		MultShootEngine();
		MultShootEngine(const MultShootEngine&);
		~MultShootEngine();
		//\}

		/**
		 * \name Operators
		 * \{
		 */
		MultShootEngine& operator =(const MultShootEngine &e);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		double getMaxErr() const;
		int getMaxIts() const;
		double getTol() const;
		bool doesFullFinalProp() const;
		bool doesLineSearch() const;
		bool isFindingEvent() const;
		MSTOF_tp getTOFType() const;

		void setDoLineSearch(bool);
		void setFullFinalProp(bool);
		void setIgnoreCrash(bool);
		void setIgnoreDiverge(bool);
		void setFindEvent(bool);
		void setMaxErr(double);
		void setMaxIts(int);
		void setSaveEachIt(bool);
		void setTol(double);
		void setTOFType(MSTOF_tp);
		
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		void multShoot(const Arcset*, Arcset*,
			MultShootData *pData = nullptr);
		void multShoot(MultShootData*);
		//\}

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void reset();
		static double getTotalDV(const MultShootData&);
		static bool finiteDiff_checkMultShoot(const Arcset*, 
			Verbosity_tp verbosity = Verbosity_tp::SOME_MSG, 
			bool writeToFile = false);
		static bool finiteDiff_checkMultShoot(const Arcset*, MultShootEngine, 
			Verbosity_tp verbosity = Verbosity_tp::SOME_MSG, 
			bool writeToFile = false);
		static void propSegsFromFreeVars(MultShootData&, SimEngine&, 
			Verbosity_tp verbosity = Verbosity_tp::NO_MSG);
		//\}
		
	private:

		/** Describe the way that time is parameterized in the design variable 
		vector */
		MSTOF_tp tofTp = MSTOF_tp::VAR_FREE;

		/** 
		 * Flag that indicates whether or not the LU solver has been 
		 * initialized for the current corrections process
		 */
		bool bInitLUSolver = false;

		/** 
		 * Flag that indicates whether or not the QR solver has been 
		 * initialized for the current corrections process
		 */
		bool bInitQRSolver = false;

		/** Whether or not to conduct final round of propagations with an 
			integrator that leverages variable step time. When set to TRUE,
			the output arcset will have Segments with many states. If set to 
			FALSE, the final arcset will have Segments with only the first
			and final state of the propagation. This latter option is best for
			applications that require speed and efficiency */
		bool bFullFinalProp = true;

		/** Flag to turn on when the algorithm is used to locate an event */
		bool bFindEvent = false;

		/** Flag to turn off crash detection in the simulation engine */
		bool bIgnoreCrash = false;

		/** 
		 * 	Flag to ignore diverge (i.e. don't throw an exception) and return 
		 * 	the partially converged iteration data instead
		 */
		bool bIgnoreDiverge = false;

		/** 
		 * Whether or not to use a rough line search to choose the size of the 
		 * Newton step. Default is false
		 */
		bool bLineSearchAttenFactor = false;

		/** Whether or not the LU solver failed */
		bool bLUFailed = false;

		/** Whether or not to save the solution after each update; 
		useful for debugging */
		bool bSaveEachIt = false;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Maximum error value permitted; if the error rises above this value, 
		the corrections are considered diverged */
		double maxErr = 1e3;

		/** 
		 * 	Maximum permissible constraint vector magnitude for convergence. 
		 *	This tolerance also influences the simulation tolerance.
		 */
		double tolF = 1e-12;

		/**
		 * Maximum permissible newton step (\f$ \delta \vec{X} \f$) magnitude 
		 * for convergence. I.e., if the design variable changes by less than 
		 * this amount, consider the process converged
		 */
		double tolX = 1e-14;

		/**
		 * Tolerance that sets the criterion to decide whether a spurious 
		 * convergence has occured during the line search
		 */
		double tolA = 1e-12;

		/**
		 * 	Line Search: Determines how large the average rate of decrease of 
		 * 	the constraint vector, F, can be realtive to the initial rate of 
		 * 	decrease, i.e.,
		 * 	\f[
		 * 		f(\vec{X}_{n+1}) \leq f(\vec{X}_n) + 
		 * 		\alpha \vec{\nabla} f(\vec{X}_n) \cdot (\vec{X}_{n+1} - 
		 * 		\vec{X}_n)\,,
		 * 	\f]
		 * 	where \f$ f = \frac{1}{2} \vec{F}^T \vec{F} \f$. For more details, 
		 * 	see the MathSpec document.
		 */
		double ls_alpha = 1e-4;

		/** 
		 * Line Search: Maximum permissible magnitude of the Newton update. 
		 * Any \f$ \vec{X} \f$ vectors with magnitudes larger than this are 
		 * scaled down
		 */
		double ls_maxStepSize = 1e2;

		/**
		 * Sparse LU Solver to factorize the Jacobian matrix when solving
		 * the update equation. 
		 */
		Eigen::SparseLU<SparseMatXCd, Eigen::COLAMDOrdering<int> > luSolver {};

		/**
		 * Sparse QR Solver to factorize the Jacobian matrix when solving
		 * the update equation.
		 */
		Eigen::SparseQR<SparseMatXCd, Eigen::COLAMDOrdering<int> > qrSolver {};

		/**
		 * \name Analysis Functions
		 * \{
		 */
		void checkDFSingularities(MatrixXRd);
		void chooseStep_LineSearch(MultShootData&, const Eigen::VectorXd&, 
			const Eigen::VectorXd&, const SparseMatXCd&, const Eigen::VectorXd&,
			Eigen::VectorXd&, bool&);
		void factorizeJacobian(const SparseMatXCd&, const Eigen::VectorXd&, 
			Eigen::VectorXd&, bool);
		void reportConMags(const MultShootData&);
		void solveUpdateEq(MultShootData&, const Eigen::VectorXd&, 
			const Eigen::VectorXd&, Eigen::VectorXd&);
		Eigen::ComputationInfo QR(const SparseMatXCd&, const Eigen::VectorXd&, 
			Eigen::VectorXd&);
		Eigen::ComputationInfo LU(const SparseMatXCd&, const Eigen::VectorXd&, 
			Eigen::VectorXd&, bool);
		//\}

		/**
		 * \name Utility Functions
		 * \{
		 */
		void cleanEngine();
		void copyMe(const MultShootEngine&);
		//\}
};

}// END of Astrohelion namespace