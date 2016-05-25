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

#ifndef H_BCR4BPR_NODESET
#define H_BCR4BPR_NODESET

#include "tpat_nodeset.hpp"

#include "matio.h"

// Forward Declarations
class TPAT_Event;
class TPAT_Sys_Data_BC4BP;

/**
 *	@brief his derivative of the TPAT_Nodeset object contains additional information
 *	for the BCR4BP
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes and epoch times at each node are also stored
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Nodeset_BC4BP : public TPAT_Nodeset{

public:
	// *structors
	TPAT_Nodeset_BC4BP(const TPAT_Sys_Data_BC4BP*);
	TPAT_Nodeset_BC4BP(const double[6], const TPAT_Sys_Data_BC4BP*, double, double, int);
	TPAT_Nodeset_BC4BP(std::vector<double>, const TPAT_Sys_Data_BC4BP*, double, double, int);
	TPAT_Nodeset_BC4BP(const double[6], const TPAT_Sys_Data_BC4BP*, double, double, int,
		tpat_nodeDistro_tp);
	TPAT_Nodeset_BC4BP(std::vector<double>, const TPAT_Sys_Data_BC4BP*, double, double, int,
		tpat_nodeDistro_tp);
	// TPAT_Nodeset_BC4BP(const TPAT_Nodeset_BC4BP&, int, int);
	TPAT_Nodeset_BC4BP(const TPAT_Nodeset_BC4BP&);
	TPAT_Nodeset_BC4BP(const TPAT_Base_Arcset&);
	baseArcsetPtr create(const TPAT_Sys_Data*) const;
	baseArcsetPtr clone() const;

	// Operators

	// Set and Get Functions

	// Utility Functions
	
private:
	void initExtraParam();
};

#endif