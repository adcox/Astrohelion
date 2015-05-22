/**
 *	The nodeset object is similar to a trajectory object, but a nodeset only contains a few
 *	distinct states, or "nodes" and is used in corrections processes to break a trajectory
 *	into smaller pieces, which can improve the corrector's performance.
 *
 *	Author: Andrew Cox
 *	Version: May 21, 2015
 */

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

#ifndef __H_NODESET__
#define __H_NODESET__

#include "adtk_constraint.hpp"

#include <vector>

// Forward Declarations
class adtk_sys_data;

class adtk_nodeset{
	public:
		/**
		 *	Node Distribution Type
		 *
		 *	This type describes how nodes are distributed, or how they are chosen
		 *	from a trajectory. Values are:
		 *
		 *	NONE 		- 	There is no organizational method; nodes may be input by
		 *					user.
		 *
		 *	TIME 		- 	nodes spread evenly in time
		 *
		 *	ARCLENGTH 	- 	nodes spread evenly along trajectory by arclength; this
		 *					method is approximate, so not all legs will be exactly 
		 *					the same length, but should be close.
		 */
		enum node_distro_t {NONE, TIME, ARCLENGTH};

		adtk_nodeset(const int);
		adtk_nodeset(const adtk_nodeset&);

		adtk_nodeset& operator =(const adtk_nodeset&);

		std::vector<double> getNode(int) const;
		double getTOF(int) const;
		int getNumNodes() const;
		int getNodeSize() const;
		node_distro_t getNodeDistro() const;
		adtk_sys_data* getSysData() const;

		void appendNode(std::vector<double>);
		void appendTOF(double);
		void setNodeDistro(node_distro_t);

	protected:
		/** The number of states in one node */
		const int nodeSize;

		/** How nodes are distributed */
		node_distro_t nodeDistro = NONE;

		/** A vector of nodes; organized in row-major order */
		std::vector<double> nodes;

		/** A vector of TOFs between nodes */
		std::vector<double> tofs;

		void initSetFromICs(double[], adtk_sys_data*, double, double, int, node_distro_t);
};

#endif