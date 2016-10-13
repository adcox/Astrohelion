/**
 *	@file Event.cpp
 *	@brief Data object that stores information about a simulation event
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

#include "Event.hpp"

#include "AsciiOutput.hpp"
#include "SysData_bc4bp.hpp"
#include "BodyData.hpp"
#include "Calculations.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_ltvp.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"
 
#include <cmath>

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Basic constructor
 *  @details Use one of the two createEvent() functions to initialize the rest
 *  of the event object
 */
Event::Event(){}

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
 *	@throws Exception if this constructor is called for an event type that requires data
 *	@throws Exception if the event type is not recognized
 */
void Event::createEvent(Event_tp t, int dir, bool willStop){
	switch(t){
		case Event_tp::YZ_PLANE:
		case Event_tp::XZ_PLANE:
		case Event_tp::XY_PLANE:
		case Event_tp::CRASH:
		{
			std::vector<double> params {0};
			initEvent(t, dir, willStop, params);
			break;
		}
		case Event_tp::JC:
		case Event_tp::APSE:
		case Event_tp::DIST:
			throw Exception("Event_tp::Event: Cannot create this type of event without parameter data...");
		default: 
			throw Exception("Event_tp::Event: Creating event with no type");
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
 *	@param params a vector of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect (otherwise it will read uninitialized memory).
 *
 *	@see Event_tp::Event_tp
 */
void Event::createEvent(Event_tp t, int dir, bool willStop, std::vector<double> params){
	initEvent(t, dir, willStop, params);
}//====================================================

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
 */
Event::Event(Event_tp t, int dir, bool willStop){
	createEvent(t, dir, willStop);
}//================================================


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
 *	@param params a vector of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect.
 *
 *	@see Event_tp::Event_tp
 *	@throws Exception if the dynamic model does not support this event type
 *	@throws Exception if the event type is not recognized
 *	@throws Exception if data values refer to invalid indices
 */
Event::Event(Event_tp t, int dir , bool willStop, std::vector<double> params){
	initEvent(t, dir, willStop, params);
}//==========================================

/**
 *	@see Event(data, t, dir, willStop, params)
 */
void Event::initEvent(Event_tp t, int dir, bool willStop, std::vector<double> params){
	// Save the parameters
	type = t;
	triggerDir = dir;
	bStop = willStop;
	paramsIn = params;	// Just temporary; initialize() function will reorder and add relevant constants

	// Determine the type of constraint that will be needed to detect this event
	switch(type){
		case Event_tp::YZ_PLANE:
		case Event_tp::XZ_PLANE:
		case Event_tp::XY_PLANE:
			conType = Constraint_tp::STATE;
			break;
		case Event_tp::CRASH:
			conType = Constraint_tp::MAX_DIST;
			break;
		case Event_tp::JC:
			conType = Constraint_tp::JC;
			break;
		case Event_tp::APSE:
			conType = Constraint_tp::APSE;
			break;
		case Event_tp::DIST:
			conType = Constraint_tp::DIST;
			break;
		default: 
			throw Exception("Event::initEvent: Creating event with no type");
	}
}//====================================================

void Event::initialize(const SysData* pSys){
	pSysData = const_cast<SysData*>(pSys);

	if(!pSysData)
		throw Exception("Event::initialize: System data pointer has not been set!");

	if(! pSys->getDynamicsModel()->supportsEvent(type)){
		throw Exception("Event::initialize: The current dynamic model does not support this event type");
	}

	double data[] = {NAN, NAN, NAN, NAN, NAN, NAN};	// six empty elements
	switch(type){
		case Event_tp::YZ_PLANE: data[0] = paramsIn[0]; break;	// x = specified value
		case Event_tp::XZ_PLANE: data[1] = paramsIn[0];	break;	// y = specified value
		case Event_tp::XY_PLANE: data[2] = paramsIn[0]; break;	// z = specified value
		case Event_tp::CRASH:
		{
			data[0] = paramsIn[0];	// Index of primary
			if(data[0] < pSysData->getNumPrimaries()){
				// Get body data, compute crash distance
			    BodyData primData(pSysData->getPrimary(static_cast<int>(data[0])));
			    data[1] = (primData.getRadius() + primData.getMinFlyBy())/pSysData->getCharL();
			}else{
				throw Exception("Cannot access primary for crash event");
			}
			break;
		}
		case Event_tp::JC: data[0] = paramsIn[0]; break;	// JC = specified value
		case Event_tp::APSE: data[0] = paramsIn[0]; break; 	// primary index = specified value
		case Event_tp::DIST:
			data[0] = paramsIn[0];
			data[1] = paramsIn[1];
			break;
		default: break;	// Do nothing
	}

	conData.clear();
	conData.insert(conData.begin(), data, data+6);
}//====================================================

/**
 *	@brief copy constructor
 */
Event::Event(const Event &ev){
	copyEvent(ev);
}//====================================================

/**
 *	@brief copy the event
 *	@param ev an event
 */
void Event::copyEvent(const Event &ev){
	type = ev.type;
	triggerDir = ev.triggerDir;
	triggerCount = ev.triggerCount;
	stopCount = ev.stopCount;
	bStop = ev.bStop;
	dist = ev.dist;
	lastDist = ev.lastDist;
	theTime = ev.theTime;
	state = ev.state;
	conType = ev.conType;
	conData = ev.conData;
	paramsIn = ev.paramsIn;
	pSysData = ev.pSysData;	// COPY ADDRESS (ptr) of SYS DATA
}//=============================================

/**
 *	@brief Destructor
 */
Event::~Event(){}

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator
 */
Event& Event::operator =(const Event &ev){
	copyEvent(ev);
	return *this;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are identical
 */
bool operator ==(const Event &lhs, const Event &rhs){
	bool same = lhs.type == rhs.type && lhs.triggerDir == rhs.triggerDir &&
		lhs.bStop == rhs.bStop;

	if(lhs.paramsIn.size() == rhs.paramsIn.size()){
		for(unsigned int i = 0; i < lhs.paramsIn.size(); i++){
			same = same && lhs.paramsIn[i] == rhs.paramsIn[i];
		}	
	}else{
		same = false;
	}

	return same;
}//====================================================

/**
 *	@brief Comparison operator
 *	@param lhs
 *	@param rhs
 *	@return true if the two events are not identical
 */
bool operator !=(const Event &lhs, const Event &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the trigger direction for this event; -1 for negative, +1
 *	for positive, 0 for both/either
 */
int Event::getDir() const { return triggerDir; }

/**
 *	@return the event type
 */
Event_tp Event::getType() const { return type; }

/**
 *	@return a human-readable string representing the event type
 */
const char* Event::getTypeStr() const{
	switch(type){
		case Event_tp::NONE: return "NONE"; break;
		case Event_tp::YZ_PLANE: return "yz-plane"; break;
		case Event_tp::XZ_PLANE: return "xz-plane"; break;
		case Event_tp::XY_PLANE: return "xy-plane"; break;
		case Event_tp::CRASH: return "crash"; break;
		case Event_tp::JC: return "jacobi constant"; break;
		case Event_tp::APSE: return "apse"; break;
		case Event_tp::DIST: return "distance"; break;
		default: return "UNDEFINED!"; break;
	}
}//========================================

/**
 *	@return the time associated with this event
 */
double Event::getTime() const { return theTime; }

/**
 *	@return a pointer to the state vector object; useful for in-place reading or writing
 */
std::vector<double>* Event::getState() { return &state; }

/**
 *	@return whether or not this event will stop the integration
 */
bool Event::stopOnEvent() const { return bStop; }

/**
 *	@return the type of constraint this event will use to target the exact event occurence
 */
Constraint_tp Event::getConType() const { return conType; }

/**
 *	@return the constraint data used to target this exact event
 */
std::vector<double> Event::getConData() const { return conData; }

/**
 *  @brief Return the system data pointer
 *  @details This pointer can only be set by runnining the initialize() function
 *  on the event or by copying another event
 *  @return the system data pointer
 */
const SysData* Event::getSysData() { return pSysData; }

/**
 *	@brief Retrieve the current trigger count, or the number of times
 *	this event has been triggered during the current simulation
 *	@return the trigger count
 */
int Event::getTriggerCount() const { return triggerCount; }

/**
 *	@brief Retrieve the number of triggers this event can have before 
 *	the simulation will be stopped (if applicable)
 *	@return the stopping trigger count
 */
int Event::getStopCount() const { return stopCount; }

/**
 *	@brief Increment the trigger counter by +1
 */
void Event::incrementCount(){ triggerCount++; }

/**
 *	@brief Set the trigger direction for this event
 *	@param d the direction: +1 for positive, -1 for negative, 0 for both/either
 */
void Event::setDir(int d){ triggerDir = d; }

/**
 *	@brief Set the number of triggers this event can endure before the simulation
 *	is forced to stop (if applicable, i.e. if stopOnEvent() = true)
 *	@param c the maximum number of triggers; simulation will be stopped when this
 *	number of triggers occurs (not after)
 */
void Event::setStopCount(int c){ stopCount = c; }

/**
 *  @brief Set the flag that determines whether a simulation ends when the event occurs
 * 
 *  @param s Whether or not the simulation should stop when this event is triggered
 */
void Event::setStopOnEvent(bool s){ bStop = s; }

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
bool Event::crossedEvent(const double y[6], double t) const{
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
void Event::updateDist(const double y[6], double t){	
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
 *	@throws Exception if the event type associated with this event is not implemented
 *	@throws Exception if the system data pointer has not been initialized via the initialize() function
 */
double Event::getDist(const double y[6], double t) const{
	if(!pSysData)
		throw Exception("Event::getDist: SysData pointer has not been initialized; please call initialize() function");

	double d = 0;
	switch(type){
		case Event_tp::YZ_PLANE: d = conData[0] - y[0]; break;
		case Event_tp::XZ_PLANE: d = conData[1] - y[1]; break;
		case Event_tp::XY_PLANE: d = conData[2] - y[2]; break;
		case Event_tp::CRASH:
		{
			std::vector<double> primPos = pSysData->getDynamicsModel()->getPrimPos(t, pSysData);

			int Pix = static_cast<int>(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		case Event_tp::JC:
		{
			const SysData_cr3bp *crSys = static_cast<const SysData_cr3bp *>(pSysData);
			d = conData[0] - DynamicsModel_cr3bp::getJacobi(y, crSys->getMu());
			break;
		}
		case Event_tp::APSE:
		{
			int Pix = static_cast<int>(conData[0]);
			d = pSysData->getDynamicsModel()->getRDot(Pix, t, y, pSysData);
			break;
		}
		case Event_tp::DIST:
		{
			std::vector<double> primPos = pSysData->getDynamicsModel()->getPrimPos(t, pSysData);
			int Pix = static_cast<int>(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		default:
			throw Exception("Event::getDist: Event type not implemented");
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
 *	@throws Exception if the event type associated with this event is not implemented
 */
int Event::getDir(const double y[6], double t) const{
	double d = 0;
	double dt = t - theTime;

	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case Event_tp::YZ_PLANE: d = y[0] - state[0]; break;
		case Event_tp::XZ_PLANE: d = y[1] - state[1]; break;
		case Event_tp::XY_PLANE: d = y[2] - state[2]; break;
		case Event_tp::CRASH:
		case Event_tp::JC:
		case Event_tp::APSE:
		case Event_tp::DIST:
			d = dist - lastDist;
			break;
		default: 
			throw Exception("Event type not implemented");
	}

	return static_cast<int>(d*dt/std::abs(d*dt));
}//====================================================

/**
 *	@brief Print out a discription of the event
 */
void Event::printStatus() const{
	printf("Event: Type = %s, Trigger Dir = %d, KillSim = %s\n", getTypeStr(), triggerDir, 
		bStop ? "YES" : "NO");
	printf("  Dist: %e Last Dist: %e\n", dist, lastDist);
}//====================================================

/**
 *  @brief Reset the event to avoid any confusion when a simulation is rerun with the same event
 */
void Event::reset(){
	triggerCount = 0;
	dist = 100000;
	lastDist = 100000;
	theTime = 0;
	state.clear();
}//====================================================






}// END of Astrohelion namespace