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

#include "tpat_linkable.hpp"

#include "tpat_constraint.hpp"

#include <vector>
// Forward Declarations


/**
 *	@brief A brief description
 *
 *	@author Andrew Cox
 *	@version 
 *	@copyright GNU GPL v3.0
 */
class tpat_node : public tpat_linkable{

public:
	// *structors
	tpat_node();
	tpat_node(const double[6], double);
	tpat_node(std::vector<double>, double);
	tpat_node(const double[6], const double[3], double);
	tpat_node(std::vector<double>, std::vector<double>, double);
	tpat_node(const tpat_node&);
	virtual ~tpat_node();

	// Operators
	tpat_node& operator =(const tpat_node&);
	friend bool operator ==(const tpat_node&, const tpat_node&);
	friend bool operator !=(const tpat_node&, const tpat_node&);

	// Set and Get functions
	void addConstraint(tpat_constraint);
	void clearConstraints();
	std::vector<double> getAccel() const;
	std::vector<tpat_constraint> getConstraints() const;
	double getEpoch() const;
	double getExtraParam(int) const;
	std::vector<double> getExtraParams() const;
	int getNumCons() const;
	std::vector<double> getState() const;
	std::vector<bool> getVelCon() const;
	void removeConstraint(int);
	void setAccel(const double*);
	void setAccel(std::vector<double>);
	void setConstraints(std::vector<tpat_constraint>);
	void setConstraintNodeNum(int);
	void setEpoch(double);
	void setExtraParam(int, double);
	void setExtraParams(std::vector<double>);
	void setState(const double*);
	void setState(std::vector<double>);
	void setVel_AllCon();
	void setVel_AllDiscon();
	void setVelCon(const bool[3]);
	void setVelCon(std::vector<bool>);
	void setVelCon(bool, bool, bool);

protected:
	virtual void copyMe(const tpat_node&);
	virtual void initArrays();

	double state[6];	//!< Stores 3 position and 3 velocity states
	double accel[3];	//!< Stores 3 acceleration states
	double epoch = 0;	//!< The epoch associated with this node, relative to some base epoch

	/** Stores extra parameters like mass, costates, etc. */
	std::vector<double> extraParam;

	/** Stores flags, which may be interpreted by derived classes */
	std::vector<bool> flags;

	/** Stores constraints on this node (especially usefull in nodesets) */
	std::vector<tpat_constraint> cons;
};

#endif