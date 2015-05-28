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
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __H_EVENT_
#define __H_EVENT_

#include "adtk_constraint.hpp"

#include <vector>

class adtk_event{
	public:
		/**
		 *	Event Type
		 *
		 *	NONE 		- 	No type has been specified; cannot be used in integration
		 *
		 *	YZ-Plane 	- 	Event occurs when trajectory crosses the YZ-plane (x = 0)
		 *					In 2D, this is equivalent to a y-axis crossing
		 *
		 *	XZ-PLANE 	- 	Event occurs when trajectory crosses the XZ-plane (y = 0)
		 *					In 2D, this is equivalent to an x-axis crossing
		 *
		 *	XY-Plane	- 	Event occurs when trajectory crosses the XY-plane (z = 0)
		 */
		enum event_t {NONE, YZ_PLANE, XZ_PLANE, XY_PLANE};

		// *structors
		adtk_event();
		adtk_event(event_t, int, bool, int);

		// Get and Set Functions
		int getDir() const;
		event_t getType() const;
		const char* getTypeStr();
		double getTime() const;
		std::vector<double>* getState();
		int getIndex() const;
		bool stopOnEvent() const;

		adtk_constraint::constraint_t getConType() const;
		std::vector<double> getConData() const;
		int getConNode() const;

		void setDir(int);
		void setType(event_t);
		void setIndex(int);

		// Computations, etc.
		bool crossedEvent(double[6]);
		void updateDist(double[6]);
	private:
		/** The type of event this is */
		event_t type = NONE;

		/** Direction of desired event crossing: +1 for positive, -1 for negative, 0 for both */
		int direction = 0;

		/** Whether or not to stop integration when this event occurs */
		bool stop = true;

		/** Distance to the event; must be able to change sign */
		double dist = 100;

		/** Time at which the even occurs */
		double theTime = 0;

		/** State at which the event occurs; also used to store last state whenver
		updateDist() is called so we can determine the direction */
		std::vector<double> state;

		/** A unique index assigned to this event by the integrator (in case
		 there are multiple similar events) */
		int index = -1;

		/** Type of constraint used by the shooting algorithm to locate this event */
		adtk_constraint::constraint_t conType;

		/** Data for the constraint used by the shooting algorithm to locate this event */
		std::vector<double> conData;

		double getDist(double[6]);
		int getDir(double[6]);
};

#endif