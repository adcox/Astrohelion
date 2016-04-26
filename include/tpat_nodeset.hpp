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

#ifndef H_NODESET
#define H_NODESET

#include "tpat_arc_data.hpp"
#include "matio.h"

// Forward Declarations
class tpat_event;
class tpat_node;
class tpat_traj;

/**
 *	@brief Similar to tpat_traj, but only holds state data at specific "nodes"
 * 
 *	The nodeset object is similar to a trajectory object, but a nodeset only contains a few
 *	distinct states, or "nodes" and is used in corrections processes to break a trajectory
 *	into smaller pieces, which can improve the corrector's performance.
 *
 *	In addition to nodes, a nodeset stores information about the constraints that should
 *	be applied when the nodeset is passed through a corrections algorithm
 *
 *	@see tpat_node
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_nodeset : public tpat_arc_data{

public:
	/**
		 *	@brief Node distribution type
		 *
		 *	Specified how nodes are distributed along an integrated trajectory
		 */
		enum tpat_nodeDistro_tp {
			DISTRO_NONE, 	//!< There is no organizational method; nodes may be input by user.
			DISTRO_TIME,	//!< Nodes spread evenly in time
			DISTRO_ARCLENGTH};	//!< Nodes spread evenly along trajectory by arclength (approx.)

	// *structors
	tpat_nodeset(const tpat_sys_data*);
	tpat_nodeset(const tpat_nodeset&);
	tpat_nodeset(const tpat_arc_data&);
	tpat_nodeset(const tpat_nodeset&, int, int);
	tpat_nodeset(const char*);
	// Operators

	// Set and Get Functions
	std::vector<tpat_constraint> getNodeCons(int) const;
	tpat_node getNode(int) const;

	int getNumCons() const;
	int getNumNodes() const;
	double getTOF(int) const;
	double getTotalTOF() const;
	
	void addConstraint(tpat_constraint);
	virtual void appendNode(tpat_node);
	int createNodesAtEvent(int, tpat_event);
	virtual int createNodesAtEvents(int, std::vector<tpat_event>);
	virtual void deleteNode(int);
	virtual void insertNode(int, tpat_node);

	void allowDV_at(std::vector<int>);
	void allowDV_all();
	void allowDV_none();

	// Utility Functions
	void clearConstraints();
	void readFromMat(const char*);
	void print() const;
	void reverseOrder();
	void saveToMat(const char*) const;
	void initExtraParam();

protected:

	void initSetFromICs(const double[6], double, double, int, tpat_nodeDistro_tp);
	void initSetFromICs_time(const double[6], double, double, int);
	void initSetFromICs_arclength(const double[6], double, double, int);
	void initSetFromTraj(tpat_traj, int, tpat_nodeDistro_tp);
	void initStepVectorFromMat(mat_t *, const char*);
	void saveTOFs(mat_t*) const;

};

#endif