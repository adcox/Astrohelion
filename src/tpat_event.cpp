/**
 *	@file tpat_event.cpp
 *	@brief Data object that stores information about a simulation event
 */
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

#include "tpat_event.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_body_data.hpp"
#include "tpat_calculations.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Basic constructor
 *  @details Use one of the two createEvent() functions to initialize the rest
 *  of the event object
 * 
 *  @param data a system data object that describes the system this event will occur in
 */
TPAT_Event::TPAT_Event(const TPAT_Sys_Data *data) : sysData(data) {}

/**
 *	@brief Create an event
 *
 *	Note that creating a CRASH event using this constructor will default to a crash
 *	with Primary #0 and a minimum acceptable distance of zero; to specify a different 
 *	primary and miss distance, use the customizable constructor.
 *
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 *	@throws TPAT_Exception if this constructor is called for an event type that requires data
 *	@throws TPAT_Exception if the event type is not recognized
 */
void TPAT_Event::createEvent(TPAT_Event_Tp t, int dir, bool willStop){
	switch(t){
		case TPAT_Event_Tp::YZ_PLANE:
		case TPAT_Event_Tp::XZ_PLANE:
		case TPAT_Event_Tp::XY_PLANE:
		case TPAT_Event_Tp::CRASH:
		{
			double params[] = {0};
			initEvent(t, dir, willStop, params);
			break;
		}
		case TPAT_Event_Tp::JC:
		case TPAT_Event_Tp::APSE:
		case TPAT_Event_Tp::DIST:
			throw TPAT_Exception("TPAT_Event_Tp::TPAT_Event: Cannot create this type of event without parameter data...");
		default: 
			throw TPAT_Exception("TPAT_Event_Tp::TPAT_Event: Creating event with no type");
	}
}//===================================================

/**
 *	@brief Create an event with custom specifications
 *	
 *	Rather than using the default parameters, this constructor allows you to
 *	create more specialized events.
 *
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 *	@param params an array of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect (otherwise it will read uninitialized memory).
 *
 *	@see TPAT_Event_Tp::TPAT_Event_Tp
 */
void TPAT_Event::createEvent(TPAT_Event_Tp t, int dir, bool willStop, double *params){
	initEvent(t, dir, willStop, params);
}//====================================================

/**
 *	@brief Create an event
 *
 *	Note that creating a CRASH event using this constructor will default to a crash
 *	with Primary #0 and a minimum acceptable distance of zero; to specify a different 
 *	primary and miss distance, use the customizable constructor.
 *
 *	@param data a system data object that describes the system this event will occur in
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 */
TPAT_Event::TPAT_Event(const TPAT_Sys_Data *data, TPAT_Event_Tp t, int dir, bool willStop) : sysData(data){
	createEvent(t, dir, willStop);
}//================================================


/**
 *	@brief Create an event with custom specifications
 *	
 *	Rather than using the default parameters, this constructor allows you to
 *	create more specialized events.
 *
 *	@param data a system data object that describes the system this event will occur in
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 *	@param params an array of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect (otherwise it will read uninitialized memory).
 *
 *	@see TPAT_Event_Tp::TPAT_Event_Tp
 *	@throws TPAT_Exception if the dynamic model does not support this event type
 *	@throws TPAT_Exception if the event type is not recognized
 *	@throws TPAT_Exception if data values refer to invalid indices
 */
TPAT_Event::TPAT_Event(const TPAT_Sys_Data *data, TPAT_Event_Tp t, int dir , bool willStop, double* params) : sysData(data){
	initEvent(t, dir, willStop, params);
}//==========================================

/**
 *	@see TPAT_Event(data, t, dir, willStop, params)
 */
void TPAT_Event::initEvent(TPAT_Event_Tp t, int dir, bool willStop, double* params){
	type = t;
	triggerDir = dir;
	stop = willStop;

	if(! sysData->getModel()->supportsEvent(type)){
		throw TPAT_Exception("TPAT_Event_Tp::initEvent: The current dynamic model does not support this event type");
	}

	// Create constraint data based on the type
	switch(type){
		case TPAT_Event_Tp::YZ_PLANE:
		case TPAT_Event_Tp::XZ_PLANE:
		case TPAT_Event_Tp::XY_PLANE:
			conType = TPAT_Constraint_Tp::STATE;
			break;
		case TPAT_Event_Tp::CRASH:
			conType = TPAT_Constraint_Tp::MAX_DIST;
			break;
		case TPAT_Event_Tp::JC:
			conType = TPAT_Constraint_Tp::JC;
			break;
		case TPAT_Event_Tp::APSE:
			conType = TPAT_Constraint_Tp::APSE;
			break;
		case TPAT_Event_Tp::DIST:
			conType = TPAT_Constraint_Tp::DIST;
			break;
		default: 
			throw TPAT_Exception("TPAT_Event_Tp::initEvent: Creating event with no type");
	}

	double data[] = {NAN, NAN, NAN, NAN, NAN, NAN};	// six empty elements
	switch(type){
		case TPAT_Event_Tp::YZ_PLANE: data[0] = params[0]; break;	// x = specified value
		case TPAT_Event_Tp::XZ_PLANE: data[1] = params[0];	break;	// y = specified value
		case TPAT_Event_Tp::XY_PLANE: data[2] = params[0]; break;	// z = specified value
		case TPAT_Event_Tp::CRASH:
		{
			data[0] = params[0];	// Index of primary
			if(data[0] < sysData->getNumPrimaries()){
				// Get body data, compute crash distance
			    TPAT_Body_Data primData(sysData->getPrimary((int)(data[0])));
			    data[1] = (primData.getRadius() + primData.getMinFlyBy())/sysData->getCharL();
			}else{
				throw TPAT_Exception("Cannot access primary for crash event");
			}
			break;
		}
		case TPAT_Event_Tp::JC: data[0] = params[0]; break;	// JC = specified value
		case TPAT_Event_Tp::APSE: data[0] = params[0]; break; 	// primary index = specified value
		case TPAT_Event_Tp::DIST:
			data[0] = params[0];
			data[1] = params[1];
			break;
		default: break;	// Do nothing
	}
	conData.insert(conData.begin(), data, data+6);
}//==========================================

/**
 *	@brief copy constructor
 */
TPAT_Event::TPAT_Event(const TPAT_Event &ev) : sysData(ev.sysData){
	copyEvent(ev);
}//==========================================

/**
 *	@brief copy the event
 *	@param ev an event
 */
void TPAT_Event::copyEvent(const TPAT_Event &ev){
	type = ev.type;
	triggerDir = ev.triggerDir;
	stop = ev.stop;
	dist = ev.dist;
	theTime = ev.theTime;
	state = ev.state;
	conType = ev.conType;
	conData = ev.conData;
	sysData = ev.sysData;	// COPY ADDRESS (ptr) of SYS DATA
	triggerCount = ev.triggerCount;
	stopCount = ev.stopCount;
}//=============================================

/**
 *	@brief Destructor
 */
TPAT_Event::~TPAT_Event(){}

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator
 */
TPAT_Event& TPAT_Event::operator =(const TPAT_Event &ev){
	copyEvent(ev);
	return *this;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are identical
 */
bool operator ==(const TPAT_Event &lhs, const TPAT_Event &rhs){
	return lhs.type == rhs.type && lhs.triggerDir == rhs.triggerDir &&
		lhs.stop == rhs.stop && lhs.sysData == rhs.sysData;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are not identical
 */
bool operator !=(const TPAT_Event &lhs, const TPAT_Event &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the trigger direction for this event; -1 for negative, +1
 *	for positive, 0 for both/either
 */
int TPAT_Event::getDir() const { return triggerDir; }

/**
 *	@return the event type
 */
TPAT_Event_Tp TPAT_Event::getType() const { return type; }

/**
 *	@return a human-readable string representing the event type
 */
const char* TPAT_Event::getTypeStr() const{
	switch(type){
		case TPAT_Event_Tp::NONE: return "NONE"; break;
		case TPAT_Event_Tp::YZ_PLANE: return "yz-plane"; break;
		case TPAT_Event_Tp::XZ_PLANE: return "xz-plane"; break;
		case TPAT_Event_Tp::XY_PLANE: return "xy-plane"; break;
		case TPAT_Event_Tp::CRASH: return "crash"; break;
		case TPAT_Event_Tp::JC: return "jacobi constant"; break;
		case TPAT_Event_Tp::APSE: return "apse"; break;
		case TPAT_Event_Tp::DIST: return "distance"; break;
		default: return "UNDEFINED!"; break;
	}
}//========================================

/**
 *	@return the time associated with this event
 */
double TPAT_Event::getTime() const { return theTime; }

/**
 *	@return a pointer to the state vector object; useful for in-place reading or writing
 */
std::vector<double>* TPAT_Event::getState() { return &state; }

/**
 *	@return whether or not this event will stop the integration
 */
bool TPAT_Event::stopOnEvent() const { return stop; }

/**
 *	@return the type of constraint this event will use to target the exact event occurence
 */
TPAT_Constraint_Tp TPAT_Event::getConType() const { return conType; }

/**
 *	@return the constraint data used to target this exact event
 */
std::vector<double> TPAT_Event::getConData() const { return conData; }

/**
 *	@brief Retrieve the current trigger count, or the number of times
 *	this event has been triggered during the current simulation
 *	@return the trigger count
 */
int TPAT_Event::getTriggerCount() const { return triggerCount; }

/**
 *	@brief Retrieve the number of triggers this event can have before 
 *	the simulation will be stopped (if applicable)
 *	@return the stopping trigger count
 */
int TPAT_Event::getStopCount() const { return stopCount; }

/**
 *	@brief Increment the trigger counter by +1
 */
void TPAT_Event::incrementCount(){ triggerCount++; }

/**
 *	@brief Set the trigger direction for this event
 *	@param d the direction: +1 for positive, -1 for negative, 0 for both/either
 */
void TPAT_Event::setDir(int d){ triggerDir = d; }

/**
 *	@brief Set the system data object for this event
 *	@param data a system data object
 */
void TPAT_Event::setSysData(TPAT_Sys_Data* data){ sysData = data; }

/**
 *	@brief Set the number of triggers this event can endure before the simulation
 *	is forced to stop (if applicable, i.e. if stopOnEvent() = true)
 *	@param c the maximum number of triggers; simulation will be stopped when this
 *	number of triggers occurs (not after)
 */
void TPAT_Event::setStopCount(int c){ stopCount = c; }

/**
 *  @brief Set the flag that determines whether a simulation ends when the event occurs
 * 
 *  @param s Whether or not the simulation should stop when this event is triggered
 */
void TPAT_Event::setStopOnEvent(bool s){ stop = s; }

//-----------------------------------------------------
//      Computations, etc.
//-----------------------------------------------------

/**
 *	@brief Determine (roughly) whether or not this event has occured between the
 *	previous trajectory state and the current one.
 *
 *	@param y the current integrated state (6 elements)
 *	@param t the current time
 *	@return whether or not the trajectory has passed through this event
 */
bool TPAT_Event::crossedEvent(const double y[6], double t) const{
	double newDist = getDist(y, t);

	// See if we have crossed (in pos. or neg. direction)
	if(newDist*dist < 0){ // have different signs
		if(triggerDir == 0){
			return true;
		}else{
			return triggerDir == getDir(y, t);
		}
	}
	return false;
}//============================================

/**
 *	@brief Update the distance variable, which will later be compared in 
 *	the <tt>crossedEvent()</tt> function to determine whether or not 
 *	the integration has crossed the event
 *
 *	@param y a 6-element state vector; if y is larger than 42 elements, 
 *	only the first 6 will be copied
 *	@param t non-dimensional time associated with state <tt>y</tt>
 */
void TPAT_Event::updateDist(const double y[6], double t){	
	// update the dist variable using information from y
	dist = getDist(y, t);

	// Save the state from y for later comparison
	state.clear();
	state.insert(state.begin(), y, y+6);
	theTime = t;
}//======================================

/**
 *	@brief Compute the distance from the input state to the event
 *	@param y a 6-element state vector representing the current integration state
 *	@param t non-dimensional time associated with state <tt>y</tt>
 *	@return the distance
 *	@throws TPAT_Exception if the event type associated with this event is not implemented
 */
double TPAT_Event::getDist(const double y[6], double t) const{
	double d = 0;
	switch(type){
		case TPAT_Event_Tp::YZ_PLANE: d = conData[0] - y[0]; break;
		case TPAT_Event_Tp::XZ_PLANE: d = conData[1] - y[1]; break;
		case TPAT_Event_Tp::XY_PLANE: d = conData[2] - y[2]; break;
		case TPAT_Event_Tp::CRASH:
		{
			std::vector<double> primPos = sysData->getModel()->getPrimPos(t, sysData);

			int Pix = (int)(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		case TPAT_Event_Tp::JC:
		{
			const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP *>(sysData);
			d = conData[0] - TPAT_Model_CR3BP::getJacobi(y, crSys->getMu());
			break;
		}
		case TPAT_Event_Tp::APSE:
		{
			int Pix = (int)(conData[0]);
			d = sysData->getModel()->getRDot(Pix, t, y, sysData);
			break;
		}
		case TPAT_Event_Tp::DIST:
		{
			std::vector<double> primPos = sysData->getModel()->getPrimPos(t, sysData);
			int Pix = (int)(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		default:
			throw TPAT_Exception("Event type not implemented");
	}

	return d;
}//====================================================

/**
 *	@brief Get the direction of propagation for the event by comparing an input state <tt>y</tt>
 *	to the state stored in the <tt>state</tt> variable (which was updated last iteration)
 *
 *	@param y a 6-element state vector
 *	@param t non-dimensional time associated with state <tt>y</tt>
 *	@return positive or negative one to correspond with the sign
 *	@throws TPAT_Exception if the event type associated with this event is not implemented
 */
int TPAT_Event::getDir(const double y[6], double t) const{
	double d = 0;
	double dt = t - theTime;

	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case TPAT_Event_Tp::YZ_PLANE: d = y[0] - state[0]; break;
		case TPAT_Event_Tp::XZ_PLANE: d = y[1] - state[1]; break;
		case TPAT_Event_Tp::XY_PLANE: d = y[2] - state[2]; break;
		case TPAT_Event_Tp::CRASH:
		case TPAT_Event_Tp::JC:
		case TPAT_Event_Tp::APSE:
		case TPAT_Event_Tp::DIST:
			d = dist - lastDist;
			break;
		default: 
			throw TPAT_Exception("Event type not implemented");
	}

	return (int)(d*dt/std::abs(d*dt));
}//====================================================

/**
 *	@brief Print out a discription of the event
 */
void TPAT_Event::printStatus() const{
	printf("Event: Type = %s, Trigger Dir = %d, KillSim = %s\n", getTypeStr(), triggerDir, 
		stop ? "YES" : "NO");
	printf("  Dist: %e Last Dist: %e\n", dist, lastDist);
}//====================================================

/**
 *  @brief Reset the event to avoid any confusion when a simulation is rerun with the same event
 */
void TPAT_Event::reset(){
	triggerCount = 0;
	dist = 100000;
	lastDist = 100000;
	theTime = 0;
	state.clear();
}//====================================================






