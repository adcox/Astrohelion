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

#ifndef __H_EVENT_
#define __H_EVENT_

#include "adtk_constraint.hpp"
#include "adtk_sys_data.hpp"

#include <vector>

/**
 *	@brief A data object containing information about an event that may
 *	occur during simulation
 *
 *	TODO: The time variable isn't used. should it be? How?
 *
 *	@author Andrew Cox
 *	@version May 29, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_event{
	public:
		/**
		 *	Event Type
		 *
		 *	NONE 		- 	No type has been specified; cannot be used in integration
		 *
		 *	YZ_Plane 	- 	Event occurs when trajectory crosses the YZ-plane (x = 0)
		 *					In 2D, this is equivalent to a y-axis crossing
		 *
		 *	XZ_PLANE 	- 	Event occurs when trajectory crosses the XZ-plane (y = 0)
		 *					In 2D, this is equivalent to an x-axis crossing
		 *
		 *	XY_PLANE	- 	Event occurs when trajectory crosses the XY-plane (z = 0)
		 *
		 *	CRASH		- 	Event occurs when trajectory falls below minimum acceptable
		 *					altitude or the surface of one of the system primaries.
		 */
		enum event_t {NONE, YZ_PLANE, XZ_PLANE, XY_PLANE, CRASH};

		// *structors
		adtk_event();
		adtk_event(adtk_sys_data*, event_t, int, bool);
		adtk_event(adtk_sys_data*, event_t, int, bool, double*);
		adtk_event(const adtk_event&);

		// Operators
		adtk_event& operator =(const adtk_event&);

		// Get and Set Functions
		int getDir() const;
		event_t getType() const;
		const char* getTypeStr() const;
		double getTime() const;
		std::vector<double>* getState();
		bool stopOnEvent() const;

		adtk_constraint::constraint_t getConType() const;
		std::vector<double> getConData() const;
		int getConNode() const;

		void setDir(int);
		void setSysData(adtk_sys_data*);

		// Computations, etc.
		bool crossedEvent(double[6], double) const;
		void updateDist(double[6], double);

		void printStatus() const;
	private:
		event_t type = NONE; //!< The type of event this is

		/** Direction of desired event crossing: +1 for positive, -1 for negative, 0 for both */
		int triggerDir = 0;

		/** Whether or not to stop integration when this event occurs */
		bool stop = true;

		double dist = 100000;		//!< Distance to the event; must be able to change sign
		double lastDist = 100000; 	//!< distance to event at previous iteration
		double theTime = 0;			//!< Time at which the even occurs

		/** State at which the event occurs; also used to store last state whenver
		updateDist() is called so we can determine the direction */
		std::vector<double> state;

		/** Type of constraint used by the shooting algorithm to locate this event */
		adtk_constraint::constraint_t conType;

		/** Data for the constraint used by the shooting algorithm to locate this event */
		std::vector<double> conData;

		adtk_sys_data* sysData; 	//!< Copy of the system data pointer

		void copyEvent(const adtk_event&);
		void initEvent(event_t, int, bool, double*);
		double getDist(double[6], double) const;
		int getDir(double[6], double) const;
};

#endif