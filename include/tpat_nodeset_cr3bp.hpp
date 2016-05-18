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
class tpat_sys_data_cr3bp;
class tpat_traj_cr3bp;

/**
 *	@brief This derivative of the tpat_nodeset contains additional information for 
 *	the BCR4BPR
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes are recorded in a separate vector.
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_nodeset_cr3bp : public tpat_nodeset{

public:
	// *structors
	tpat_nodeset_cr3bp(const tpat_sys_data_cr3bp*);
	tpat_nodeset_cr3bp(const double[6], const tpat_sys_data_cr3bp*, double, int);
	tpat_nodeset_cr3bp(std::vector<double>, const tpat_sys_data_cr3bp*, double, int);
	tpat_nodeset_cr3bp(const double[6], const tpat_sys_data_cr3bp*, double, int, tpat_nodeDistro_tp);
	tpat_nodeset_cr3bp(std::vector<double>, const tpat_sys_data_cr3bp*, double, int, tpat_nodeDistro_tp);
	tpat_nodeset_cr3bp(tpat_traj_cr3bp, int);
	tpat_nodeset_cr3bp(tpat_traj_cr3bp, int, tpat_nodeDistro_tp);
	tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp&, int, int);
	tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp&);
	tpat_nodeset_cr3bp(const tpat_base_arcset&);
	baseArcsetPtr create(const tpat_sys_data*) const;
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