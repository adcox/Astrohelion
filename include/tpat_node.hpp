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

#ifndef H_TPAT_NODE
#define H_TPAT_NODE

#include "tpat_arc_step.hpp"

/**
 *	@brief Derived from tpat_arc_step with specific calls for nodes
 *
 *	Values Stored in ExtraParam:
 *	* 0 	- 	Node time-of-flight
 *	* 1 	- 	Epoch (BCR4BP)
 */
class tpat_node : public tpat_arc_step{
public:
	// *structors
	tpat_node(double*, double);
	tpat_node(std::vector<double>, double);
	tpat_node(const tpat_node&);
	tpat_node(const tpat_arc_step&);
	
	// Set and Get Functions
	double getTOF() const;
	std::vector<bool> getVelCon() const;

	void setTOF(double);
	void setVel_AllCon();
	void setVel_AllDiscon();
	void setVelCon(bool[3]);
	void setVelCon(std::vector<bool>);
	void setVelCon(bool, bool, bool);

private:
	void initArrays();
	
};

#endif