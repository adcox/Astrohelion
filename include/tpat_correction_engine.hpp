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

#ifndef H_CORRECTIONS
#define H_CORRECTIONS

#include "tpat.hpp"
 
#include "tpat_constants.hpp"
#include "tpat_constraint.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_traj.hpp"

#include <vector>

// Forward declarations
class TPAT_Constraint;
class TPAT_MultShoot_Data;
class TPAT_Nodeset_BC4BP;
class TPAT_Nodeset_CR3BP;
class TPAT_Sim_Engine;
class TPAT_Sys_Data_BC4BP;



/**
 *	@brief An engine object to perform corrections, such as multiple shooting.
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
 *	@author Andrew Cox
 *	@version August 3, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Correction_Engine : public TPAT{
	public:
		// *structors
		/** Default, do-nothing constructor */
		TPAT_Correction_Engine(){}
		TPAT_Correction_Engine(const TPAT_Correction_Engine&);
		~TPAT_Correction_Engine();

		// Operators
		TPAT_Correction_Engine& operator =(const TPAT_Correction_Engine &e);

		// Set and get functions
		int getMaxIts() const;
		double getTol() const;
		TPAT_Verbosity_Tp isVerbose() const;
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
		void setTol(double);
		void setVarTime(bool);
		void setVerbose(TPAT_Verbosity_Tp);
		
		// Utility/Action functions
		TPAT_MultShoot_Data multShoot(const TPAT_Nodeset*, TPAT_Nodeset*);
		TPAT_MultShoot_Data multShoot(TPAT_MultShoot_Data, TPAT_Nodeset*);

	private:
		/** Describes how many messages to spit out */
		TPAT_Verbosity_Tp verbose = TPAT_Verbosity_Tp::SOME_MSG;

		/** Whether or not to use variable time in the corrections process */
		bool varTime = true;

		/** Whether or not to force all arcs to have the same length (in time);
		 	only applies if variable time is enabled */
		bool equalArcTime = false;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Maximum acceptable error value, non-dimensional units.
			This tolerance also influences the simulation tolerance.
		 */
		double tol = 1e-12;

		/** Flag to turn on when this algorithm is being used to locate an event */
		bool findEvent = false;

		/** Flag to turn off crash detection in the simulation engine */
		bool ignoreCrash = false;

		/** Flag to ignore diverge (i.e. don't throw an exception) and return the partially converged iteration data instead */
		bool ignoreDiverge = false;

		/** Flag to apply scaling to variables, constraint values, and partial derivatives to ease numerical processes */
		bool scaleVars = false;

		/** Whether or not the engine is ready to be cleaned and/or deconstructed */
		bool isClean = true;

		void cleanEngine();
		void copyEngine(const TPAT_Correction_Engine&);
		void reportConMags(const TPAT_MultShoot_Data*);
		Eigen::VectorXd solveUpdateEq(TPAT_MultShoot_Data*);
};

#endif