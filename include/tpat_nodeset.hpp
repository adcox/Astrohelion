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

#include "tpat_constraint.hpp"
#include "tpat_node.hpp"

#include "matio.h"

#include <vector>

// Forward Declarations
class tpat_sys_data;
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
 *	@version May 21, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_nodeset{
	public:
		/**
		 *	@brief Node distribution type
		 *
		 *	Specified how nodes are distributed along an integrated trajectory
		 */
		enum node_distro_t {
			DISTRO_NONE, 	//!< There is no organizational method; nodes may be input by user.
			DISTRO_TIME,	//!< Nodes spread evenly in time
			DISTRO_ARCLENGTH};	//!< Nodes spread evenly along trajectory by arclength (approx.)

		// *structors
		tpat_nodeset();
		tpat_nodeset(const tpat_nodeset&);
		tpat_nodeset(const tpat_nodeset&, int, int);

		virtual ~tpat_nodeset();
		
		// Operators
		tpat_nodeset& operator =(const tpat_nodeset&);

		// Set and Get functions
		tpat_node getNode(int) const;
		double getTOF(int) const;
		double getTotalTOF() const;
		int getNumCons() const;
		int getNumNodes() const;
		node_distro_t getNodeDistro() const;
		tpat_constraint getConstraint(int) const;

		/**
		 *	Extend this function to return a system data object from derivative classes
		 *	@return a pointer to the system data object
		 */
		virtual tpat_sys_data* getSysData() = 0;

		void addConstraint(tpat_constraint);
		void appendNode(tpat_node);
		void deleteNode(int);
		void insertNode(int, tpat_node);

		void setNodeDistro(node_distro_t);
		void setVelConNodes_allBut(std::vector<int>);

		// Utility Functions
		void clearConstraints();
		void copyMe(const tpat_nodeset&);
		virtual void print() const = 0;		//!< @brief Output a human-readable description of the nodeset
		void reverseOrder();
		void saveToMat(const char*);

	protected:

		/** How nodes are distributed */
		node_distro_t nodeDistro = DISTRO_NONE;

		/** A vector of node objects */
		std::vector<tpat_node> nodes;

		/** Vector of constraints to be applied to this nodeset*/
		std::vector<tpat_constraint> constraints;

		static void basicConcat(const tpat_nodeset&, const tpat_nodeset&, tpat_nodeset*);
		void initSetFromICs(double[6], tpat_sys_data*, double, double, int, node_distro_t);
		void initSetFromTraj(tpat_traj, tpat_sys_data*, int, node_distro_t);
		void saveNodes(mat_t*);
		void saveTOFs(mat_t*);
};

#endif