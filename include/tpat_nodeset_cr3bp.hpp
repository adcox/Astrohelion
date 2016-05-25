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

#ifndef H_CR3BP_NODESET
#define H_CR3BP_NODESET

#include "tpat_nodeset.hpp"

// Forward Declarations
class TPAT_Sys_Data_CR3BP;
class TPAT_Traj_CR3BP;

/**
 *	@brief This derivative of the TPAT_Nodeset contains additional information for 
 *	the BCR4BPR
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes are recorded in a separate vector.
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Nodeset_CR3BP : public TPAT_Nodeset{

public:
	// *structors
	TPAT_Nodeset_CR3BP(const TPAT_Sys_Data_CR3BP*);
	TPAT_Nodeset_CR3BP(const double[6], const TPAT_Sys_Data_CR3BP*, double, int);
	TPAT_Nodeset_CR3BP(std::vector<double>, const TPAT_Sys_Data_CR3BP*, double, int);
	TPAT_Nodeset_CR3BP(const double[6], const TPAT_Sys_Data_CR3BP*, double, int, tpat_nodeDistro_tp);
	TPAT_Nodeset_CR3BP(std::vector<double>, const TPAT_Sys_Data_CR3BP*, double, int, tpat_nodeDistro_tp);
	TPAT_Nodeset_CR3BP(TPAT_Traj_CR3BP, int);
	TPAT_Nodeset_CR3BP(TPAT_Traj_CR3BP, int, tpat_nodeDistro_tp);
	TPAT_Nodeset_CR3BP(const TPAT_Nodeset_CR3BP&, int, int);
	TPAT_Nodeset_CR3BP(const TPAT_Nodeset_CR3BP&);
	TPAT_Nodeset_CR3BP(const TPAT_Base_Arcset&);
	baseArcsetPtr create(const TPAT_Sys_Data*) const;
	baseArcsetPtr clone() const;

	// Operators

	// Set and Get 
	double getJacobi(int) const;
	double getJacobiByIx(int) const;
	void setJacobi(int, double);
	void setJacobiByIx(int, double);

protected:
	// Utility
	void initExtraParam();
};

#endif