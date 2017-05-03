/**
 *  \file Arcset.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version April 28, 2017
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

#include <vector>

#include "BaseArcset.hpp"

#include "matio.h"

namespace astrohelion{

// Forward Declarations
class Event;
class SysData;
class Arcset;

class Arcset : public BaseArcset{

public:
	/**
	 *	\brief Node distribution type
	 *
	 *	Specified how nodes are distributed along an integrated trajectory
	 */
	enum NodeDistro_tp {
		NONE, 		//!< There is no organizational method; nodes may be input by user.
		TIME,		//!< Nodes spread evenly in time
		ARCLENGTH	//!< Nodes spread evenly along trajectory by arclength (approx.)
	};

	/**
	 *  \name *structors
	 *  \{
	 */
	Arcset(const SysData*);
	Arcset(const Arcset&);
	Arcset(const BaseArcset&);
	Arcset(const char*);

	virtual ~Arcset();
	virtual baseArcsetPtr create(const SysData*) const;
	virtual baseArcsetPtr clone() const;
	//\}

	/**
	 * \name Operators
	 * \{
	 */
	friend Arcset operator +(const Arcset&, const Arcset&);
	virtual Arcset& operator +=(const Arcset&);
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	void allowDV_at(std::vector<int>);
	void allowDV_all();
	void allowDV_none();
	int createNodesAtEvent(int, Event, double minTimeDiff = 1e-2);
	virtual int createNodesAtEvents(int, std::vector<Event>, double minTimeDiff = 1e-2);

	double getTimeByIx(int) const;

	void setTimeByIx(int, double);
	void shiftAllTimes(double);
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	virtual void print() const;
	void readFromMat(const char*);
	void saveToMat(const char*) const;
	//\}

protected:
	virtual void saveCmds(mat_t*) const;
	virtual void readCmds(mat_t*);
	
	void initFromICs(std::vector<double>, double, double, int, NodeDistro_tp type = NodeDistro_tp::TIME, unsigned int ctrlLawID = ControlLaw::NO_CTRL);
	void initFromICs_time(std::vector<double>, double, double, int, unsigned int ctrlLawID = ControlLaw::NO_CTRL);
	void initFromICs_arclength(std::vector<double>, double, double, int, unsigned int ctrlLawID = ControlLaw::NO_CTRL);
	void initFromTraj(Arcset, int, NodeDistro_tp);
};

}// End of astrohelion namespace