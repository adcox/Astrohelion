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

#ifndef H_EVENT
#define H_EVENT

#include "tpat_constraint.hpp"
#include "tpat_sys_data.hpp"

#include <vector>

/**
 *	@brief A data object containing information about an event that may
 *	occur during simulation
 *
 *	**Adding a New Event**
 *	* Add a new enumerated type value with detailed documentation of how it will
 *		be parsed.
 * 	* Update the initEvent() function
 *	* Update the getConTypeStr() function
 *	* Define or identify a constraint type that can locate this event to a high
 *		degree of accuracy.
 *	* Update the getDist() function
 * 	* Update the getDir() function
 *
 *	@author Andrew Cox
 *	@version August 3, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_event{
	public:
		/**
		 *	@brief The type of event
		 *
		 *	This tells the simulation and correction engines how to interpret the
		 *	data stored in this object.
		 *
		 *	For the _PLANE crossing events, pass a single double into the <tt>params</tt>
		 *	field to specify the location of the plane. For example, choosing an event
		 *	with type YZ_PLANE and passing a <tt>param</tt> of 0.5 will create an event
		 *	that fires when the trajectory passes through the <tt>x=0.5</tt> plane. By 
		 * 	default, 0 is used, so the planes include the axes.
		 */
		enum event_t {
			NONE,		//!< No type has been specified; cannot be used in integration
			YZ_PLANE,	/*!< Event occurs when trajectory crosses an YZ-plane. By default, 
						 *	this plane occurs at x = 0, but a custom x-coordinate may be 
						 *	specified by placing the value of that coordinate in the 
						 *	<tt>param</tt> array.
						 */
			XZ_PLANE,	/*!< Event occurs when trajectory crosses an XZ-plane. By default, 
						 *	this plane occurs at y = 0, but a custom y-coordinate may be 
						 *	specified by placing the value of that coordinate in the 
						 *	<tt>param</tt> array.
						 */
			XY_PLANE,	/*!< Event occurs when trajectory crosses an XY-plane. By default, 
						 *	this plane occurs at z = 0, but a custom z-coordinate may be 
						 *	specified by placing the value of that coordinate in the 
						 *	<tt>param</tt> array.
						 */
			CRASH,		/*!< Event occurs when trajectory falls below minimum acceptable
		 				 * 	altitude or the surface of one of the system primaries.
		 				 *	The <tt>param</tt> array should have the first element specifying the 
						 *	primary index (0 for P1, 1 for P2, etc.) The minimum acceptable radius
						 *	will be the radius of the primary plus the minimum acceptable fly-by distance
						 *	specified in the tpat_body_data class.
		 				 */
		 	JC, 		/*!< Event occurs when the Jacobi value reaches the specified value
		 				 * 	of Jacobi Constant. Place this JC value in the first element of
		 				 * 	the <tt>params</tt> vector present in the 
		 				 * 	tpat_event(tpat_sys_data*, event_t, int, bool, double*) constructor.
		 				 * 	This event can only be supported by dynamic models that have associated
		 				 * 	system data objects that can be cast to cr3bp system data objects.
		 				 */
		 	APSE,		/*!< Event occurs when an apse is reached. The <tt>param</tt> array should have
		 				 * 	the first element specifiying the primary index (0 for P1, 1 for P2, etc.).
		 				 * 	The <tt>direction</tt> of the event identifies what type of apse, i.e.
		 				 *	0 will catch all apsides, -1 will catch only apopases, and +1 will catch 
		 				 *	only periapses
		 				 */
		 	DIST 		/*!< Event occurs when a distance from a primary is reached. The <tt>param</tt>
		 				 * 	array should have two elements: element 0 specifies the primary index, and
		 				 * 	element 1 specifies the distance from the center of the primary in
		 				 *	non-dimensional units. The <tt>direction</tt> of the event identifies the
		 				 *	direction of travel relative to the primary; +1 is away from the primary,
		 				 *	-1 is towards the primary, and 0 triggers for either.
		 				 */
		};

		// *structors
		tpat_event(const tpat_sys_data*);
		tpat_event(const tpat_sys_data*, event_t, int, bool);
		tpat_event(const tpat_sys_data*, event_t, int, bool, double*);
		tpat_event(const tpat_event&);
		void createEvent(event_t, int, bool);
		void createEvent(event_t, int, bool, double*);
		~tpat_event();
		
		// Operators
		tpat_event& operator =(const tpat_event&);
		friend bool operator ==(const tpat_event&, const tpat_event&);
		friend bool operator !=(const tpat_event&, const tpat_event&);
		
		// Get and Set Functions
		std::vector<double> getConData() const;
		tpat_constraint::constraint_t getConType() const;
		int getDir() const;
		double getTime() const;
		int getStopCount() const;
		int getTriggerCount() const;
		event_t getType() const;
		const char* getTypeStr() const;
		
		std::vector<double>* getState();
		void incrementCount();
		bool stopOnEvent() const;

		void setDir(int);
		void setType(event_t);
		void setStopCount(int);
		void setSysData(tpat_sys_data*);

		// Computations, etc.
		bool crossedEvent(const double[6], double) const;
		void updateDist(const double[6], double);

		void printStatus() const;
	private:

		event_t type = NONE; //!< The type of event this is

		/** Direction of desired event crossing: +1 for positive, -1 for negative, 0 for both */
		int triggerDir = 0;

		/** Number of times this event has been triggered */
		int triggerCount = 0;

		/** The number of times this event can be triggered before the simulation
		 	is stopped (if applicable, i.e. stop = true)*/
		int stopCount = 1;

		/** Whether or not to stop integration when this event occurs */
		bool stop = true;

		double dist = 100000;		//!< Distance to the event; must be able to change sign
		double lastDist = 100000; 	//!< distance to event at previous iteration
		double theTime = 0;			//!< Time at which the even occurs

		/** State at which the event occurs; also used to store last state whenver
		updateDist() is called so we can determine the direction */
		std::vector<double> state;

		/** Type of constraint used by the shooting algorithm to locate this event */
		tpat_constraint::constraint_t conType;

		/** Data for the constraint used by the shooting algorithm to locate this event */
		std::vector<double> conData;

		const tpat_sys_data* sysData; 	//!< Copy of the system data pointer

		void copyEvent(const tpat_event&);
		void initEvent(event_t, int, bool, double*);
		double getDist(const double[6], double) const;
		int getDir(const double[6], double) const;
};

#endif