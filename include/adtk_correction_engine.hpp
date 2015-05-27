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

/**
 *	An engine objecto to perform corrections, such as single-shooting
 *	and multiple shooting algorithms
 *
 *	Author: Andrew Cox
 *	Version: May 25, 2015
 */
class adtk_correction_engine{
	public:
		
		void correct_cr3bp(adtk_cr3bp_nodeset*);
		void correct_bcr4bpr(adtk_bcr4bpr_nodeset*);

		// Set and get functions
		bool usesVarTime() const;
		int getMaxIts() const;
		double getTol() const;
		bool isVerbose() const;

		void setVarTime(bool);
		void setMaxIts(int);
		void setTol(double);
		void setVerbose(bool);

	private:
		/** Whether or not to spit out lots of messages */
		bool verbose = true;

		/** Whether or not to use variable time in the corrections process */
		bool varTime = true;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Maximum acceptable error value, non-dimensional units */
		double tol = 1e-12;

		/** */
		adtk_nodeset *nodeset;

		void correct(adtk_nodeset*);
};

#endif