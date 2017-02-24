/**
 *  @file Nodeset.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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

#include "BaseArcset.hpp"
#include "matio.h"

namespace astrohelion{

// Forward Declarations
class Event;
class Node;
class Traj;

/**
 *	@ingroup traj
 *	@brief Similar to Traj, but only holds state data at specific "nodes"
 * 
 *	The nodeset object is similar to a trajectory object, but a nodeset only contains a few
 *	distinct states, or "nodes" and is used in corrections processes to break a trajectory
 *	into smaller pieces, which can improve the corrector's performance.
 *
 *	In addition to nodes, a nodeset stores information about the constraints that should
 *	be applied when the nodeset is passed through a corrections algorithm
 *
 *	@see Node
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class Nodeset : public BaseArcset{

public:
	/**
	 *	@brief Node distribution type
	 *
	 *	Specified how nodes are distributed along an integrated trajectory
	 */
	enum NodeDistro_tp {
		NONE, 		//!< There is no organizational method; nodes may be input by user.
		TIME,		//!< Nodes spread evenly in time
		ARCLENGTH	//!< Nodes spread evenly along trajectory by arclength (approx.)
	};

	/**
	 *  @name *structors
	 *  @{
	 */
	Nodeset(const SysData*);
	Nodeset(const Nodeset&);
	Nodeset(const BaseArcset&);
	Nodeset(const Nodeset&, int, int);
	Nodeset(const char*);
	virtual ~Nodeset();
	virtual baseArcsetPtr create(const SysData*) const;
	virtual baseArcsetPtr clone() const;
	//@}

	// Operators
	friend Nodeset operator +(const Nodeset&, const Nodeset&);
	virtual Nodeset& operator +=(const Nodeset&);

	/**
	 *  @name Set and Get Functions
	 *  @{
	 */
	void allowDV_at(std::vector<int>);
	void allowDV_all();
	void allowDV_none();
	int createNodesAtEvent(int, Event, double minTimeDiff = 1e-2);
	virtual int createNodesAtEvents(int, std::vector<Event>, double minTimeDiff = 1e-2);
	//@}
	
	// Utility Functions
	virtual void readFromMat(const char*);
	virtual void print() const;
	void reverseOrder();
	virtual void saveToMat(const char*) const;
	virtual void initExtraParam();

protected:

	void initFromICs(std::vector<double>, double, double, int, NodeDistro_tp);
	void initFromICs_time(std::vector<double>, double, double, int);
	void initFromICs_arclength(std::vector<double>, double, double, int);
	void initFromTraj(Traj, int, NodeDistro_tp);

};


}