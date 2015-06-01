/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __H_CORRECTIONS_
#define __H_CORRECTIONS_

#include "adtk_nodeset.hpp"

// Forward declarations
class adtk_cr3bp_nodeset;
class adtk_bcr4bpr_nodeset;
class adtk_simulation_engine;

/**
 *	@brief An engine object to perform corrections, such as single-shooting
 *	and multiple shooting algorithms
 *
 *	@author Andrew Cox
 *	@version May 25, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_correction_engine{
	public:
		// *structors
		/** Default, do-nothing constructor */
		adtk_correction_engine(){}
		adtk_correction_engine(const adtk_correction_engine&);
		~adtk_correction_engine();

		// Operators
		adtk_correction_engine& operator =(const adtk_correction_engine &e);

		// Set and get functions
		bool usesVarTime() const;
		int getMaxIts() const;
		double getTol() const;
		bool isVerbose() const;
		bool isFindingEvent() const;
		adtk_cr3bp_nodeset getCR3BPOutput();
		adtk_bcr4bpr_nodeset getBCR4BPROutput();

		void setVarTime(bool);
		void setMaxIts(int);
		void setTol(double);
		void setVerbose(bool);
		void setFindEvent(bool);

		// Utility/Action functions
		void correct_cr3bp(adtk_cr3bp_nodeset*);
		void correct_bcr4bpr(adtk_bcr4bpr_nodeset*);

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
		adtk_nodeset *nodeset_in = 0;

		/** The output nodeset, constructed from the corrected arcs */
		adtk_nodeset *nodeset_out = 0;

		void copyEngine(const adtk_correction_engine&);
		void correct(adtk_nodeset*);
		void createOutput(std::vector<double>, adtk_simulation_engine);
};

#endif