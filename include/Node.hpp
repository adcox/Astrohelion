/**
 *  \file Node.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <cmath>
#include <vector>
#include <map>

#include "Linkable.hpp"

#include "Constraint.hpp"
#include "Event.hpp"

 
// Forward Declarations
// (none)
 
namespace astrohelion{

/**
 *	\ingroup traj
 *	\brief A single point or state on an arc
 *
 *	\author Andrew Cox
 *	\version September 30, 2016
 *	\copyright GNU GPL v3.0
 */
class Node : public Linkable{

public:
	/**
	 *  \name *structors
	 *  \{
	 */
	Node();
	Node(const double*, unsigned int, double);
	Node(std::vector<double>, double);
	Node(const Node&);
	// ~Node();
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	Node& operator =(const Node&);
	friend bool operator ==(const Node&, const Node&);
	friend bool operator !=(const Node&, const Node&);
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	std::vector<Constraint> getConstraints() const;
	double getEpoch() const;
	double getExtraParam(std::string) const;
	std::vector<double> getExtraParamVec(std::string) const;
	std::map<std::string, double> getExtraParams() const;
	std::map<std::string, std::vector<double> > getExtraParamVec() const;
	unsigned int getNumCons() const;
	std::vector<double> getState() const;
	const std::vector<double>& getStateRef_const() const;
	Event_tp getTriggerEvent() const;
	
	void setConstraints(std::vector<Constraint>);
	void setEpoch(double);
	void setExtraParam(std::string, double);
	void setExtraParamVec(std::string, std::vector<double>);
	void setExtraParams(std::map<std::string, double>);
	void setExtraParamVec(std::map<std::string, std::vector<double> >);
	void setID(int) override;
	void setState(const double*, unsigned int);
	void setState(std::vector<double>);
	void setTriggerEvent(Event_tp);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void addConstraint(Constraint);
	void clearConstraints();
	void removeConstraint(int);
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	void print() const;
	//\}

protected:
	std::vector<double> state = {};		//!< Stores states, e.g., position, velocity, mass, thrust
	double epoch = 0;	//!< The epoch associated with this node, relative to some base epoch

	/** Stores extra parameters (scalars) like mass, costates, etc. */
	std::map<std::string, double> extraParam {};

	/** Stores extra parameter vectors, like dqdt, etc. */
	std::map<std::string, std::vector<double> > extraParamVecs {};

	/** Stores constraints on this node (especially usefull in nodesets) */
	std::vector<Constraint> cons {};

	/** The event type that triggered the creation of this node */
	Event_tp triggerEventTp = Event_tp::NONE;

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	virtual void copyMe(const Node&);
	//\}
};

}// END of Astrohelion namespace