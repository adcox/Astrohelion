/**
 *  @file Node.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "Linkable.hpp"

#include "Constraint.hpp"

#include <cmath>
#include <vector>
#include <map>
 
// Forward Declarations
// (none)
 
namespace astrohelion{
/**
 *	@brief A brief description
 *
 *	@author Andrew Cox
 *	@version 
 *	@copyright GNU GPL v3.0
 */
class Node : public Linkable{

public:
	// *structors
	Node();
	Node(const double[6], double);
	Node(std::vector<double>, double);
	Node(const double[6], const double[3], double);
	Node(std::vector<double>, std::vector<double>, double);
	Node(const Node&);
	// ~Node();

	// Operators
	Node& operator =(const Node&);
	friend bool operator ==(const Node&, const Node&);
	friend bool operator !=(const Node&, const Node&);

	// Set and Get functions
	void addConstraint(Constraint);
	void clearConstraints();
	std::vector<double> getAccel() const;
	std::vector<Constraint> getConstraints() const;
	double getEpoch() const;
	double getExtraParam(std::string) const;
	std::vector<double> getExtraParamVec(std::string) const;
	std::map<std::string, double> getExtraParams() const;
	std::map<std::string, std::vector<double> > getExtraParamVec() const;
	int getNumCons() const;
	std::vector<double> getState() const;

	void print() const;

	void removeConstraint(int);
	void setAccel(const double*);
	void setAccel(std::vector<double>);
	void setConstraints(std::vector<Constraint>);
	void setConstraintID(int);
	void setEpoch(double);
	void setExtraParam(std::string, double);
	void setExtraParamVec(std::string, std::vector<double>);
	void setExtraParams(std::map<std::string, double>);
	void setExtraParamVec(std::map<std::string, std::vector<double> >);
	void setState(const double*);
	void setState(std::vector<double>);

protected:
	virtual void copyMe(const Node&);

	double state[6] = {NAN, NAN, NAN, NAN, NAN, NAN};	//!< Stores 3 position and 3 velocity states
	double accel[3] = {NAN, NAN, NAN};					//!< Stores 3 acceleration states
	double epoch = 0;	//!< The epoch associated with this node, relative to some base epoch

	/** Stores extra parameters (scalars) like mass, costates, etc. */
	// std::vector<double> extraParam {};
	std::map<std::string, double> extraParam {};

	/** Stores extra parameter vectors, like dqdt, etc. */
	std::map<std::string, std::vector<double> > extraParamVecs {};

	/** Stores constraints on this node (especially usefull in nodesets) */
	std::vector<Constraint> cons {};
};

}// END of Astrohelion namespace