/**
 *  \file MultShootEngine.hpp
 *	\brief Performs corrections on a nodeset subject to a set of constraints
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
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
 *	\brief An engine object to perform corrections, such as multiple shooting.
 *
 *	### Multiple Shooting Algorithm
 *	The multiple shooting algorithm is initiated by calling the multShoot() function. 
 *	Typical use calls the version which requires a nodeset pointer as the only input.
 *	
 *	Best Practices:
 *	
 *	- Variable scaling may be turned on to (hopefully) improve numerical performance when
 *	some variables have different orders of magnitude. Basic testing has revealed that
 *	this scaling has a positive effect when the scaling constants are relatively mild
 *	(e.g., with magnitudes between 10^-2 and 10). If larger scaling factors are required,
 *	it may be disadvantageous to scale the variables. In this case, it may be necessary to 
 *	reformulate the problem to shift large variables.
 *	
 *	\author Andrew Cox
 *	\version August 3, 2015
 *	\copyright GNU GPL v3.0
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

		// Operators
		MultShootEngine& operator =(const MultShootEngine &e);

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		int getMaxIts() const;
		double getTol() const;
		bool doesFullFinalProp() const;
		bool doesLineSearch() const;
		bool isFindingEvent() const;
		MSTOF_tp getTOFType() const;

		void setAttenuation(double, double limit = 1e-8);
		void setDoLineSearch(bool);
		void setFullFinalProp(bool);
		void setIgnoreCrash(bool);
		void setIgnoreDiverge(bool);
		void setFindEvent(bool);
		void setMaxIts(int);
		void setTol(double);
		void setTOFType(MSTOF_tp);
		
		//\}

		/**
		 *  \name Correction Algorithms
		 *  \{
		 */
		MultShootData multShoot(const Arcset*, Arcset*);
		MultShootData multShoot(MultShootData);
		//\}

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void reset();
		static double getTotalDV(const MultShootData*);
		static bool finiteDiff_checkMultShoot(const Arcset*, Verbosity_tp verbosity = Verbosity_tp::SOME_MSG, bool writeToFile = false);
		static bool finiteDiff_checkMultShoot(const Arcset*, MultShootEngine, Verbosity_tp verbosity = Verbosity_tp::SOME_MSG, bool writeToFile = false);
		//\}
		
	private:

		MSTOF_tp tofTp = MSTOF_tp::VAR_FREE;

		/**
		 * Whether or not to do a line search to choose the best step size; this slows the engine down
		 * but can increase accuracy!
		 */
		bool bDoLineSearch = false;

		/** Whether or not to conduct final round of propagations with an 
			integrator that leverages variable step time. When set to TRUE,
			the output arcset will have Segments with many states. If set to 
			FALSE, the final arcset will have Segments with only the first
			and final state of the propagation. This latter option is best for
			applications that require speed and efficiency */
		bool bFullFinalProp = true;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Scale step size by this value */
		double attenuation = 1;

		/** Do not scale steps size if error is below this tolerance */
		double attenuationLimitTol = 1e-8;

		/** Maximum acceptable error value, non-dimensional units.
			This tolerance also influences the simulation tolerance.
		 */
		double tol = 1e-12;

		/** Flag to turn on when this algorithm is being used to locate an event */
		bool bFindEvent = false;

		/** Flag to turn off crash detection in the simulation engine */
		bool bIgnoreCrash = false;

		/** Flag to ignore diverge (i.e. don't throw an exception) and return the partially converged iteration data instead */
		bool bIgnoreDiverge = false;

		void checkDFSingularities(MatrixXRd);
		void chooseStep_LineSearch(MultShootData*, const Eigen::VectorXd*, const Eigen::VectorXd*, const Eigen::VectorXd*, Eigen::VectorXd*);
		void cleanEngine();
		void copyMe(const MultShootEngine&);
		void propSegsFromFreeVars(MultShootData*, SimEngine*);
		void reportConMags(const MultShootData*);
		void solveUpdateEq(MultShootData*, const Eigen::VectorXd*, const Eigen::VectorXd*, Eigen::VectorXd*);
};

}// END of Astrohelion namespace