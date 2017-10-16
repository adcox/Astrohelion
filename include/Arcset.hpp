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

/**
 *	\ingroup traj
 *	\brief Extends BaseArcset in a few areas.
 *	\details This class is largely unneeded; may replace BaseArcset with Arcset
 *	or visa-versa; there is no need for an extra inheritance layer here
 *
 *	\author Andrew Cox
 *	\version June 9, 2017
 *	\copyright GNU GPL v3.0
 */
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
	double getTimeByIx(int) const;

	void setSTMs_cumulative();
	void setSTMs_sequence();
	void setTimeByIx(int, double);
	
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void allowDV_at(std::vector<int>);
	void allowDV_all();
	void allowDV_none();
	int createNodesAtEvent(int, Event, double minTimeDiff = 1e-2);
	virtual int createNodesAtEvents(int, std::vector<Event>, double minTimeDiff = 1e-2);
	void shiftAllTimes(double);
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	virtual void print() const;
	//\}

	/**
	 *  \name File I/O
	 *  \{
	 */
	void readFromMat(const char*, std::vector<ControlLaw*>&);
	void readFromStruct(matvar_t*, unsigned int, std::vector<ControlLaw*>&);
	void saveToMat(const char*, Save_tp saveTp = Save_tp::SAVE_ALL) const;
	void saveToStruct(matvar_t*, unsigned int, Save_tp saveTp = Save_tp::SAVE_ALL) const;
	//\}

protected:
	/**
	 *  \name File I/O
	 *  \{
	 */
	virtual void saveCmds_toFile(mat_t*, Save_tp saveTp = Save_tp::SAVE_ALL) const;
	virtual void saveCmds_toStruct(matvar_t*, unsigned int, Save_tp saveTp = Save_tp::SAVE_ALL) const;
	virtual void readCmds_fromFile(mat_t*, std::vector<ControlLaw*>&);
	virtual void readCmds_fromStruct(matvar_t*, unsigned int, std::vector<ControlLaw*>&);
	//\}
};

}// End of astrohelion namespace