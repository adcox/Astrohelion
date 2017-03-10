/**
 *  \file CorrectionEngine.hpp
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

#include "Core.hpp"
#include "Engine.hpp"

#include "Common.hpp"
#include "Constraint.hpp"
#include "EigenDefs.hpp"
#include "Node.hpp"
#include "Nodeset.hpp"
#include "SysData.hpp"
#include "Traj.hpp"

#include <vector>


namespace astrohelion{

// Forward declarations
class Constraint;
class MultShootData;
class Nodeset_bc4bp;
class Nodeset_cr3bp;
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
class CorrectionEngine : public Core, public Engine{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		/** Default, do-nothing constructor */
		CorrectionEngine(){}
		CorrectionEngine(const CorrectionEngine&);
		~CorrectionEngine();
		//\}

		// Operators
		CorrectionEngine& operator =(const CorrectionEngine &e);

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		int getMaxIts() const;
		double getTol() const;
		bool isFindingEvent() const;
		bool usesEqualArcTime() const;
		bool usesScaledVars() const;
		bool usesVarTime() const;

		void setEqualArcTime(bool);
		void setIgnoreCrash(bool);
		void setIgnoreDiverge(bool);
		void setFindEvent(bool);
		void setMaxIts(int);
		void setScaleVars(bool);
		void setStepScale(double, double limit = 1e-8);
		void setTol(double);
		void setVarTime(bool);
		//\}

		/**
		 *  \name Correction Algorithms
		 *  \{
		 */
		MultShootData multShoot(const Nodeset*, Nodeset*);
		MultShootData multShoot(MultShootData, Nodeset*);
		//\}

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void reset();
		//\}
		
	private:

		/** Whether or not to use variable time in the corrections process */
		bool bVarTime = true;

		/** Whether or not to force all arcs to have the same length (in time);
		 	only applies if variable time is enabled */
		bool bEqualArcTime = false;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Scale step size by this value */
		double stepScale = 1;

		/** Do not scale steps size if error is below this tolerance */
		double stepScaleLimitTol = 1e-8;

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

		/** Flag to apply scaling to variables, constraint values, and partial derivatives to ease numerical processes */
		bool bScaleVars = false;

		void cleanEngine();
		void copyMe(const CorrectionEngine&);
		void reportConMags(const MultShootData*);
		Eigen::VectorXd solveUpdateEq(MultShootData*);
};

}// END of Astrohelion namespace