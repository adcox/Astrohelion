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
#include "tpat.hpp"

#include "tpat_event.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
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

/** @brief Default, do-nothing constructor */
tpat_event::tpat_event(){}

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
tpat_event::tpat_event(tpat_sys_data *data, event_t t, int dir, bool willStop){
	// Create constraint data based on the type by calling the more detailed constructor
	sysData = data;

	switch(t){
		case YZ_PLANE:
		case XZ_PLANE:
		case XY_PLANE:
		case CRASH:
		{
			double params[] = {0};
			initEvent(t, dir, willStop, params);
			break;
		}
		case JC:
		case APSE:
		case DIST:
			throw tpat_exception("tpat_event::tpat_event: Cannot create this type of event without parameter data...");
		default: 
			throw tpat_exception("tpat_event::tpat_event: Creating event with no type");
	}
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
 *	@see tpat_event::event_t
 */
tpat_event::tpat_event(tpat_sys_data *data, event_t t, int dir , bool willStop, double* params){
	sysData = data;
	initEvent(t, dir, willStop, params);
}//==========================================

/**
 *	@see tpat_event(data, t, dir, willStop, params)
 */
void tpat_event::initEvent(event_t t, int dir, bool willStop, double* params){
	type = t;
	triggerDir = dir;
	stop = willStop;

	if(! sysData->getModel()->supportsEvent(type)){
		throw tpat_exception("tpat_event::initEvent: The current dynamic model does not support this event type");
	}

	// Create constraint data based on the type
	switch(type){
		case YZ_PLANE:
		case XZ_PLANE:
		case XY_PLANE:
			conType = tpat_constraint::STATE;
			break;
		case CRASH:
			conType = tpat_constraint::MAX_DIST;
			break;
		case JC:
			conType = tpat_constraint::JC;
			break;
		case APSE:
			conType = tpat_constraint::APSE;
			break;
		case DIST:
			conType = tpat_constraint::DIST;
			break;
		default: 
			throw tpat_exception("tpat_event::initEvent: Creating event with no type");
	}

	double data[] = {NAN, NAN, NAN, NAN, NAN, NAN};	// six empty elements
	switch(type){
		case YZ_PLANE: data[0] = params[0]; break;	// x = specified value
		case XZ_PLANE: data[1] = params[0];	break;	// y = specified value
		case XY_PLANE: data[2] = params[0]; break;	// z = specified value
		case CRASH:
		{
			data[0] = params[0];	// Index of primary
			if(data[0] < sysData->getNumPrimaries()){
				// Get body data, compute crash distance
			    tpat_body_data primData(sysData->getPrimary((int)(data[0])));
			    data[1] = (primData.getRadius() + primData.getMinFlyBy())/sysData->getCharL();
			}else{
				throw tpat_exception("Cannot access primary for crash event");
			}
			break;
		}
		case JC: data[0] = params[0]; break;	// JC = specified value
		case APSE: data[0] = params[0]; break; 	// primary index = specified value
		case DIST:
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
tpat_event::tpat_event(const tpat_event &ev){
	copyEvent(ev);
}//==========================================

/**
 *	@brief copy the event
 *	@param ev an event
 */
void tpat_event::copyEvent(const tpat_event &ev){
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
tpat_event::~tpat_event(){
	state.clear();
	conData.clear();
}//=============================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator
 */
tpat_event& tpat_event::operator =(const tpat_event &ev){
	copyEvent(ev);
	return *this;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are identical
 */
bool operator ==(const tpat_event &lhs, const tpat_event &rhs){
	return lhs.type == rhs.type && lhs.triggerDir == rhs.triggerDir &&
		lhs.stop == rhs.stop && lhs.sysData == rhs.sysData;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are not identical
 */
bool operator !=(const tpat_event &lhs, const tpat_event &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the trigger direction for this event; -1 for negative, +1
 *	for positive, 0 for both/either
 */
int tpat_event::getDir() const { return triggerDir; }

/**
 *	@return the event type
 */
tpat_event::event_t tpat_event::getType() const { return type; }

/**
 *	@return a human-readable string representing the event type
 */
const char* tpat_event::getTypeStr() const{
	switch(type){
		case NONE: return "NONE"; break;
		case YZ_PLANE: return "yz-plane"; break;
		case XZ_PLANE: return "xz-plane"; break;
		case XY_PLANE: return "xy-plane"; break;
		case CRASH: return "crash"; break;
		case JC: return "jacobi constant"; break;
		case APSE: return "apse"; break;
		case DIST: return "distance"; break;
		default: return "UNDEFINED!"; break;
	}
}//========================================

/**
 *	@return the time associated with this event
 */
double tpat_event::getTime() const { return theTime; }

/**
 *	@return a pointer to the state vector object; useful for in-place reading or writing
 */
std::vector<double>* tpat_event::getState() { return &state; }

/**
 *	@return whether or not this event will stop the integration
 */
bool tpat_event::stopOnEvent() const { return stop; }

/**
 *	@return the type of constraint this event will use to target the exact event occurence
 */
tpat_constraint::constraint_t tpat_event::getConType() const { return conType; }

/**
 *	@return the constraint data used to target this exact event
 */
std::vector<double> tpat_event::getConData() const { return conData; }

/**
 *	@brief Retrieve the current trigger count, or the number of times
 *	this event has been triggered during the current simulation
 *	@return the trigger count
 */
int tpat_event::getTriggerCount() const { return triggerCount; }

/**
 *	@brief Retrieve the number of triggers this event can have before 
 *	the simulation will be stopped (if applicable)
 *	@return the stopping trigger count
 */
int tpat_event::getStopCount() const { return stopCount; }

/**
 *	@brief Increment the trigger counter by +1
 */
void tpat_event::incrementCount(){ triggerCount++; }

/**
 *	@brief Set the trigger direction for this event
 *	@param d the direction: +1 for positive, -1 for negative, 0 for both/either
 */
void tpat_event::setDir(int d){ triggerDir = d; }

/**
 *	@brief Set the system data object for this event
 *	@param data a system data object
 */
void tpat_event::setSysData(tpat_sys_data* data){ sysData = data; }

/**
 *	@brief Set the number of triggers this event can endure before the simulation
 *	is forced to stop (if applicable, i.e. if stop = true)
 *	@param c the maximum number of triggers; simulation will be stopped when this
 *	number of triggers occurs (not after)
 */
void tpat_event::setStopCount(int c){ stopCount = c; }

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
bool tpat_event::crossedEvent(const double y[6], double t) const{
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
void tpat_event::updateDist(const double y[6], double t){	
	// update the dist variable using information from y
	lastDist = dist;
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
 */
double tpat_event::getDist(const double y[6], double t) const{
	double d = 0;
	switch(type){
		case YZ_PLANE: d = conData[0] - y[0]; break;
		case XZ_PLANE: d = conData[1] - y[1]; break;
		case XY_PLANE: d = conData[2] - y[2]; break;
		case CRASH:
		{
			std::vector<double> primPos = sysData->getModel()->getPrimPos(t, sysData);

			int Pix = (int)(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		case JC:
		{
			tpat_sys_data_cr3bp *crSys = static_cast<tpat_sys_data_cr3bp *>(sysData);
			d = conData[0] - cr3bp_getJacobi(y, crSys->getMu());
			break;
		}
		case APSE:
		{
			std::vector<double> primPos = sysData->getModel()->getPrimPos(t, sysData);
			int Pix = (int)(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = dx*y[3] + dy*y[4] + dz*y[5];
			break;
		}
		case DIST:
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
			throw tpat_exception("Event type not implemented");
	}

	return d;
}//=====================================

/**
 *	@brief Get the direction of propagation for the event by comparing an input state <tt>y</tt>
 *	to the state stored in the <tt>state</tt> variable (which was updated last iteration)
 *
 *	@param y a 6-element state vector
 *	@param t non-dimensional time associated with state <tt>y</tt>
 *	@return positive or negative one to correspond with the sign
 */
int tpat_event::getDir(const double y[6], double t) const{
	double d = 0;
	double dt = t - theTime;

	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case YZ_PLANE: d = y[0] - state[0]; break;
		case XZ_PLANE: d = y[1] - state[1]; break;
		case XY_PLANE: d = y[2] - state[2]; break;
		case CRASH:
		case JC:
		case APSE:
		case DIST:
			d = dist - lastDist;
			break;
		default: 
			throw tpat_exception("Event type not implemented");
	}

	return (int)(d*dt/std::abs(d*dt));
}//==============================================

/**
 *	@brief Print out a discription of the event
 */
void tpat_event::printStatus() const{
	printf("Event: Type = %s, Trigger Dir = %d, KillSim = %s\n", getTypeStr(), triggerDir, 
		stop ? "YES" : "NO");
	printf("  Dist: %e Last Dist: %e\n", dist, lastDist);
}//======================================


