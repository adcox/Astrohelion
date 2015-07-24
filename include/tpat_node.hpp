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

#include <cmath>
#include <vector>

/**
 *	@brief A basic structure to contain info about a node in a nodeset
 */
class tpat_node{
public:
	// *structors
	tpat_node();
	tpat_node(double*, double);
	tpat_node(std::vector<double>, double);
	tpat_node(const tpat_node&);
	~tpat_node();

	// Operators
	tpat_node& operator =(const tpat_node&);
	friend bool operator ==(const tpat_node&, const tpat_node&);

	// Set and Get functions
	std::vector<double> getPosVelState() const;
	double getExtraParam(int) const;
	std::vector<double> getExtraParams() const;
	double getTOF() const;
	std::vector<bool> getVelCon() const;

	void setExtraParam(int, double);
	void setExtraParams(std::vector<double>);
	void setPosVelState(double*);
	void setPosVelState(std::vector<double>);
	void setTOF(double);
	void setVel_AllCon();
	void setVel_AllDiscon();
	void setVelCon(bool[3]);
	void setVelCon(bool, bool, bool);

	// Utility Functions
	void copyMe(const tpat_node&);
	void initArrays();
	
private:

	/** Stores the 6 position and velocity states */
	double posVelState[6];

	/** Determines the length of time (non-dimensional) the corrector will integrate
	 *	past this node */
	double tof = 0;

	/** Stores any extra parameters, like epoch, initial mass, thrust level, etc. */
	std::vector<double> extraParam;

	/** Determines which velocity states are continuous. By default, 
	 *	all of them are. 1 = continuous, 0 = discontinuous */
	bool velCon[3];
};

#endif