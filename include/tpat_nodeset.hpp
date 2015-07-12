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

#ifndef __H_NODESET__
#define __H_NODESET__

#include "tpat_constraint.hpp"
#include "matio.h"

#include <vector>

// Forward Declarations
class tpat_sys_data;

/**
 *	@brief Similar to tpat_trajectory, but only holds state data at specific "nodes"
 * 
 *	The nodeset object is similar to a trajectory object, but a nodeset only contains a few
 *	distinct states, or "nodes" and is used in corrections processes to break a trajectory
 *	into smaller pieces, which can improve the corrector's performance.
 *
 *	Each node has 6 states: three position, and three velocity. Times-of-flight between nodes
 *	are stored in non-dimensional time units. In addition to nodes and times-of-flight, 
 *	a nodeset stores information about velocity continuity and a list of constraints.
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
			NONE, 	//!< There is no organizational method; nodes may be input by user.
			TIME,	//!< Nodes spread evenly in time
			ARCLENGTH};	//!< Nodes spread evenly along trajectory by arclength (approx.)

		// *structors
		tpat_nodeset(const int);
		tpat_nodeset(const tpat_nodeset&);
		tpat_nodeset(const tpat_nodeset&, int, int);

		virtual ~tpat_nodeset();
		
		// Operators
		tpat_nodeset& operator =(const tpat_nodeset&);

		// Set and Get functions
		std::vector<double>* getNodes();
		std::vector<double>* getTOFs();
		double getTotalTOF() const;
		std::vector<double> getNode(int) const;
		double getTOF(int) const;
		int getNumNodes() const;
		int getNodeSize() const;
		node_distro_t getNodeDistro() const;
		std::vector<int> getVelConNodes();
		int getNumCons() const;
		tpat_constraint getConstraint(int) const;

		/**
		 *	Extend this function to return a system data object from derivative classes
		 *	@return a pointer to the system data object
		 */
		virtual tpat_sys_data* getSysData() = 0;

		void addConstraint(tpat_constraint);
		void appendNode(std::vector<double>);
		void appendNode(double*);
		void appendTOF(double);
		void deleteTOF(int);
		void insertNode(int, std::vector<double>);
		void insertNode(int, double*);
		void insertTOF(int, double);

		void setNodeDistro(node_distro_t);
		void setTOF(int, double);
		void setVelConNodes(std::vector<int>);
		void setVelConNodes_allBut(std::vector<int>);

		// Utility Functions
		void reverseOrder();
		void clearConstraints();
		virtual void print() const = 0;		//!< @brief Output a human-readable description of the nodeset
		void saveToMat(const char*);

	protected:
		/** The number of states in one node */
		const int nodeSize;

		/** How nodes are distributed */
		node_distro_t nodeDistro = NONE;

		/** A vector of nodes; organized in row-major order */
		std::vector<double> nodes;

		/** A vector of TOFs between nodes */
		std::vector<double> tofs;

		/** Vector of constraints to be applied to this nodeset*/
		std::vector<tpat_constraint> constraints;

		/** List of node indices; the nodes included in this list should have continuous velocity */
		std::vector<int>velConNodes;

		/** Whether or not velocity continuity nodes have been initiailized or set by user */
		bool velConSet = false;

		static void basicConcat(const tpat_nodeset&, const tpat_nodeset&, tpat_nodeset*);
		void initSetFromICs(double[], tpat_sys_data*, double, double, int, node_distro_t);
		void saveNodes(mat_t*);
		void saveTOFs(mat_t*);
		// std::vector<double> concatNodes(std::vector<double>, std::vector<double>);
};

#endif