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
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __H_CONSTRAINT__
#define __H_CONSTRAINT__
 
#include <string>
#include <vector>

/**
 *	@brief Contains information about how a particular node should be 
 *	constrained during a corrections process
 *
 *	@author Andrew Cox
 *	@version May 21, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_constraint{
	public:
		/**
		 *	Constraint Types
		 *
		 *	STATE 		-	Constrain specific states to specified numeric values.
		 *					The <tt>node</tt> value represents the node to constraint,
		 *					and the <tt>data</tt> vector contains values for each state.
		 *					to leave a state unconstrained, place a NAN value in the
		 *					place of that state.
		 *
		 *	MATCH_ALL 	- 	Make one node match another in all states.
		 *					The <tt>node</tt> value represents the source node, and the
		 *					first value in <tt>data</tt> represents the node that will
		 *					be matched to the first. Note that, for non-autonomous systems,
		 *					the corrector will not attempt to match epoch time because it
		 *					doesn't make sense to make two nodes occur at the same time.
		 *
		 *	MATCH_CUST 	-	Make one node match another in the specified states.
		 *					The <tt>node</tt> value represents the source node; place the
		 *					index of the constrained node in the place of each state you
		 *					wish to constrain. For example, to constrain node 5 to match 
		 *					node 2 in x and z position, the data vector would contain
		 *					the following values: {5, NAN, 5, NAN, NAN, NAN}
		 *
		 *	DIST 		- 	Constrain a node to be at a specific distance from a primmary.
		 *					The <tt>node</tt> field identifies the constrained node, the 
		 *					first value in the <tt>data</tt> field (i.e. <tt>data[0]</tt>) 
		 *					identifies the primary number (index 0) and the second value in 
		 *					<tt>data</tt> gives the desired distance from the primary's center
		 *					in non-dimensional units
		 *
		 *	MIN_DIST 	- 	Constrain a node to be at or above a specific distance relative
		 *					to a specific primary. Data is input in the same way as the DIST
		 *					consraint described above.
		 *
		 *	MAX_DIST 	- 	Constrain a node to be at or below a specific distance relative to 
		 *					specific primary. Data is input in the same way as the DIST constraint
		 *					described above.
		 *
		 *	MAX_DELTA_V	- 	Constrain total delta V to be less than the specified amount
		 *					The <tt>node</tt> value is unused in this constraint; make the
		 *					first value in the data vector be the maximum delta-V in 
		 *					non-dimensional velocity units
		 *
		 *	DELTA_V 	- 	Constrain delta V to be exactly the specified amount
		 *					Data is input in the same format as MAX_DELTA_V
		 *
		 *	SP 			- 	Constrain the node to intersect the saddle point (BCR4BPR only)
		 *					Place the index of the node you want to constrain in the <tt>node</tt>
		 *					variable; <tt>data</tt> is unused in this constraint.
		 */
		enum constraint_t {NONE, STATE, MATCH_ALL, MATCH_CUST, DIST, MIN_DIST, MAX_DIST, 
			DELTA_V, MAX_DELTA_V, SP};
		
		adtk_constraint(int);
		adtk_constraint(int, constraint_t);
		adtk_constraint(int, constraint_t, int, std::vector<double>);
		adtk_constraint(int, constraint_t, int, double*);
		adtk_constraint(const adtk_constraint&);

		adtk_constraint& operator =(const adtk_constraint&);

		constraint_t getType() const;
		int getNode() const;
		int getNodeSize() const;
		std::vector<double> getData() const;
		int countConstrainedStates() const;

		void setType(constraint_t);
		void setNode(int);
		void setData(std::vector<double>);

		void print() const;
		const char* getTypeStr(constraint_t) const;
	private:
		const int nodeSize;

		constraint_t type = NONE;
		int node = 0;
		std::vector<double> data;
};

#endif