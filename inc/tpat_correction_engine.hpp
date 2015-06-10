/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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

#ifndef __H_CORRECTIONS_
#define __H_CORRECTIONS_

#include "tpat_nodeset.hpp"
#include "tpat_sys_data.hpp"

// Forward declarations
class tpat_bcr4bpr_nodeset;
class tpat_bcr4bpr_sys_data;
class tpat_constraint;
class tpat_cr3bp_nodeset;
class tpat_matrix;
class tpat_simulation_engine;
class iterationData;

/**
 *	@brief An engine object to perform corrections, such as single-shooting
 *	and multiple shooting algorithms
 *
 *	@author Andrew Cox
 *	@version May 25, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_correction_engine{
	public:
		// *structors
		/** Default, do-nothing constructor */
		tpat_correction_engine(){}
		tpat_correction_engine(const tpat_correction_engine&);
		~tpat_correction_engine();

		// Operators
		tpat_correction_engine& operator =(const tpat_correction_engine &e);

		// Set and get functions
		bool usesVarTime() const;
		int getMaxIts() const;
		double getTol() const;
		bool isVerbose() const;
		bool isFindingEvent() const;
		tpat_cr3bp_nodeset getCR3BPOutput();
		tpat_bcr4bpr_nodeset getBCR4BPROutput();

		void setVarTime(bool);
		void setMaxIts(int);
		void setTol(double);
		void setVerbose(bool);
		void setFindEvent(bool);

		// Utility/Action functions
		void correct_cr3bp(tpat_cr3bp_nodeset*);
		void correct_bcr4bpr(tpat_bcr4bpr_nodeset*);

	private:
		/** Whether or not to spit out lots of messages */
		bool verbose = true;

		/** Whether or not to use variable time in the corrections process */
		bool varTime = true;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Maximum acceptable error value, non-dimensional units */
		double tol = 1e-12;

		/** Whether or not an input nodeset was supplied and corrected */
		bool receivedNodesetIn = false;

		/** Whether or not space has been dynamically allocated for nodeset_out */
		bool createdNodesetOut = false;

		/** Flag to turn on when this algorithm is being used to locate an event */
		bool findEvent = false;

		/** The input nodeset, not modified*/
		tpat_nodeset *nodeset_in = 0;

		/** The output nodeset, constructed from the corrected arcs */
		tpat_nodeset *nodeset_out = 0;

		void copyEngine(const tpat_correction_engine&);
		void correct(tpat_nodeset*);
		void createOutput(iterationData*);

		void createPosVelCons(iterationData*, tpat_sys_data::system_t, int);
		void targetState(iterationData*, tpat_constraint, int);
		void targetMatchAll(iterationData*, tpat_constraint, int);
		void targetMatchCust(iterationData*, tpat_constraint, int);
		void targetDist(iterationData*, tpat_constraint, tpat_sys_data*, int, int);
		void targetSP(iterationData*, tpat_bcr4bpr_sys_data*, int);
		void updateDeltaVCon(iterationData*, tpat_sys_data*, int);
		void updatePrimPos(iterationData*, tpat_sys_data*, double);
		void updatePrimVel(iterationData*, tpat_sys_data*, double);

		tpat_matrix solveUpdateEq(iterationData*);
};

#endif