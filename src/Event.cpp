/**
 *	\file Event.cpp
 *	\brief Data object that stores information about a simulation event
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

#include "Event.hpp"

#include "AsciiOutput.hpp"
#include "SysData_bc4bp.hpp"
#include "BodyData.hpp"
#include "Calculations.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"
 
#include <cmath>

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  \brief Basic constructor
 *  \details Use one of the two createEvent() functions to initialize the rest
 *  of the event object
 */
Event::Event(){}

/**
 *	\brief Create an event
 *
 *	Note that creating a CRASH event using this constructor will default to a crash
 *	with Primary #0 and a minimum acceptable distance of zero; to specify a different 
 *	primary and miss distance, use the customizable constructor.
 *
 *	\param t the event type
 *	\param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	\param willStop whether or not this event should stop the integration
 */
Event::Event(Event_tp t, int dir, bool willStop){
	createEvent(t, dir, willStop);
}//================================================


/**
 *	\brief Create an event with custom specifications
 *	
 *	Rather than using the default parameters, this constructor allows you to
 *	create more specialized events.
 *
 *	\param t the event type
 *	\param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	\param willStop whether or not this event should stop the integration
 *	\param params a vector of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect.
 *
 *	@see Event_tp::Event_tp
 *	\throws Exception if the dynamic model does not support this event type
 *	\throws Exception if the event type is not recognized
 *	\throws Exception if data values refer to invalid indices
 */
Event::Event(Event_tp t, int dir , bool willStop, std::vector<double> params){
	initEvent(t, dir, willStop, params);
}//==========================================

/**
 *	\brief copy constructor
 */
Event::Event(const Event &ev){
	copyMe(ev);
}//====================================================

/**
 *	\brief Destructor
 */
Event::~Event(){}

/**
 *	\brief Create an event
 *
 *	Note that creating a CRASH event using this constructor will default to a crash
 *	with Primary #0 and a minimum acceptable distance of zero; to specify a different 
 *	primary and miss distance, use the customizable constructor.
 *
 *	\param t the event type
 *	\param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	\param willStop whether or not this event should stop the integration
 *	\throws Exception if this constructor is called for an event type that requires data
 *	\throws Exception if the event type is not recognized
 */
void Event::createEvent(Event_tp t, int dir, bool willStop){
	switch(t){
		case Event_tp::SIM_TOF:
		case Event_tp::SIM_COMPTIME:
		case Event_tp::SIM_ERR:
		case Event_tp::YZ_PLANE:
		case Event_tp::XZ_PLANE:
		case Event_tp::XY_PLANE:
		case Event_tp::CRASH:
		case Event_tp::MASS:
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
 *	\brief Create an event with custom specifications
 *	
 *	Rather than using the default parameters, this constructor allows you to
 *	create more specialized events.
 *
 *	\param t the event type
 *	\param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	\param willStop whether or not this event should stop the integration
 *	\param params a vector of doubles that give the constructor extra information. No
 *	specific size is required, but params must have at least as many elements as the 
 *	event type will expect (otherwise it will read uninitialized memory).
 *
 *	@see Event_tp::Event_tp
 */
void Event::createEvent(Event_tp t, int dir, bool willStop, std::vector<double> params){
	initEvent(t, dir, willStop, params);
}//====================================================
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
		case Event_tp::SIM_TOF:
		case Event_tp::SIM_COMPTIME:
		case Event_tp::SIM_ERR:
			break;	// leave the default: conType = Constraint_tp::NONE;
		case Event_tp::YZ_PLANE:
		case Event_tp::XZ_PLANE:
		case Event_tp::XY_PLANE:
		case Event_tp::MASS:
			conType = Constraint_tp::ENDSEG_STATE;
			break;
		case Event_tp::CRASH:
			conType = Constraint_tp::ENDSEG_MAX_DIST;
			break;
		case Event_tp::JC:
			conType = Constraint_tp::ENDSEG_JC;
			break;
		case Event_tp::APSE:
			conType = Constraint_tp::ENDSEG_APSE;
			break;
		case Event_tp::DIST:
			conType = Constraint_tp::ENDSEG_DIST;
			break;
		default: 
			throw Exception("Event::initEvent: Creating event with no type");
	}
}//====================================================

/**
 *  \brief Initialize the event
 *  \details Store parameters from <code>paramsIn</code> in a newly 
 *  constructed Constraint object and compute any other required values
 * 
 *  \param pSys pointer to a system data object
 */
void Event::initialize(const SysData* pSys){
	pSysData = const_cast<SysData*>(pSys);

	if(!pSysData)
		throw Exception("Event::initialize: System data pointer has not been set!");

	if(! pSys->getDynamicsModel()->supportsEvent(type)){
		throw Exception("Event::initialize: The current dynamic model does not support this event type");
	}

	// double data[] = {NAN, NAN, NAN, NAN, NAN, NAN};	// six empty elements
	std::vector<double> data(pSysData->getDynamicsModel()->getCoreStateSize(), NAN);
	switch(type){
		case Event_tp::YZ_PLANE: data[0] = paramsIn[0]; break;	// x = specified value
		case Event_tp::XZ_PLANE: data[1] = paramsIn[0];	break;	// y = specified value
		case Event_tp::XY_PLANE: data[2] = paramsIn[0]; break;	// z = specified value
		case Event_tp::CRASH:
		{
			data[0] = paramsIn[0];	// Index of primary
			if(data[0] < pSysData->getNumPrimaries()){
				// Get body data, compute crash distance
			    BodyData primData(pSysData->getPrimID(static_cast<int>(data[0])));
			    data[1] = (primData.getBodyRad() + primData.getMinFlyBy())/pSysData->getCharL();
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
		case Event_tp::MASS:
			data.assign(7,NAN);
			data[6] = paramsIn[0];	// mass = specified value
			break;
		default: break;	// Do nothing
	}

	conData.clear();
	conData.insert(conData.begin(), data.begin(), data.end());
}//====================================================

//---------------------------------------------------------------------------------------
//      Operator Functions
//---------------------------------------------------------------------------------------

/**
 *	\brief Assignment operator
 */
Event& Event::operator =(const Event &ev){
	copyMe(ev);
	return *this;
}//====================================================

/**
 *	\brief Comparison operator
 *	\param lhs left-hand-side object
 *	\param rhs right-hand-side object
 *	\return true if the two events are identical
 */
bool operator ==(const Event &lhs, const Event &rhs){
	bool same = lhs.type == rhs.type &&
		lhs.triggerDir == rhs.triggerDir &&
		lhs.bStop == rhs.bStop &&
		lhs.paramsIn == rhs.paramsIn;

	return same;
}//====================================================

/**
 *	\brief Comparison operator
 *	\param lhs left-hand-side object
 *	\param rhs right-hand-side object
 *	\return true if the two events are not identical
 */
bool operator !=(const Event &lhs, const Event &rhs){
	return !(lhs == rhs);
}//====================================================

//---------------------------------------------------------------------------------------
//      Set and Get Functions
//---------------------------------------------------------------------------------------

/**
 *	\brief Retrieve the trigger direction for this Event
 *	\return the trigger direction for this event; -1 for negative, +1
 *	for positive, 0 for both/either
 */
int Event::getTriggerDir() const { return triggerDir; }

/**
 *	\brief Retrieve the type associated with this Event object
 *	\return the event type
 */
Event_tp Event::getType() const { return type; }

/**
 *  \brief Retrieve a human-readable string representing the event type
 *	\return a human-readable string representing the event type
 */
const char* Event::getTypeStr() const{ return getEventTpStr(type); }

/**
 *  \brief Determine whether or not this event will stop the integration
 *	\return whether or not this event will stop the integration
 */
bool Event::stopOnEvent() const { return bStop; }

/**
 *  \brief Retrieve the type of constraint this event will use to target the exact event occurence
 *	\return the type of constraint this event will use to target the exact event occurence
 */
Constraint_tp Event::getConType() const { return conType; }

/**
 *  \brief Retrieve the constraint data used to target this exact event
 *	\return the constraint data used to target this exact event
 */
std::vector<double> Event::getConData() const { return conData; }

/**
 *  \brief Return the system data pointer
 *  \details This pointer can only be set by runnining the initialize() function
 *  on the event or by copying another event
 *  \return the system data pointer
 */
const SysData* Event::getSysData() { return pSysData; }

/**
 *	\brief Retrieve the current trigger count, or the number of times
 *	this event has been triggered during the current simulation
 *	\return the trigger count
 */
int Event::getTriggerCount() const { return triggerCount; }

/**
 *	\brief Retrieve the number of triggers this event can have before 
 *	the simulation will be stopped (if applicable)
 *	\return the stopping trigger count
 */
int Event::getStopCount() const { return stopCount; }

/**
 *	\brief Increment the trigger counter by +1
 */
void Event::incrementCount(){ triggerCount++; }

/**
 *	\brief Set the trigger direction for this event
 *	\param d the direction: +1 for positive, -1 for negative, 0 for both/either
 */
void Event::setTriggerDir(int d){ triggerDir = d; }

/**
 *	\brief Set the number of triggers this event can endure before the simulation
 *	is forced to stop (if applicable, i.e. if stopOnEvent() = true)
 *	\param c the maximum number of triggers; simulation will be stopped when this
 *	number of triggers occurs (not after)
 */
void Event::setStopCount(int c){ stopCount = c; }

/**
 *  \brief Set the flag that determines whether a simulation ends when the event occurs
 * 
 *  \param s Whether or not the simulation should stop when this event is triggered
 */
void Event::setStopOnEvent(bool s){ bStop = s; }

//---------------------------------------------------------------------------------------
//      Event Computations
//---------------------------------------------------------------------------------------

/**
 *	\brief Determine (roughly) whether or not this event has occured between the
 *	previous trajectory state and the current one.
 *
 *	\param y the current integrated state (6 elements)
 *	\param len number of elements in <tt>y</tt>
 *	\param t non-dimensional time associated with state <tt>y</tt>
 *	\param tDir direction that time is being propagated: +1 for forward, -1 for negative
 *	\return whether or not the trajectory has passed through this event
 */
bool Event::crossedEvent(const double* y, unsigned int len, double t, int tDir) const{
	double newDist = getDist(y, len, t);

	// See if we have crossed (in pos. or neg. direction)
	if(newDist*dist < 0){ // have different signs
		if(triggerDir == 0){
			return true;
		}else{
			return triggerDir == getDir(tDir);
		}
	}
	return false;
}//============================================

/**
 *	\brief Update the distance variable, which will later be compared in 
 *	the <tt>crossedEvent()</tt> function to determine whether or not 
 *	the integration has crossed the event
 *
 *	\param y the state vector; must have at least the core states
 *	\param len number of elements stored in y
 *	\param t non-dimensional time associated with state <tt>y</tt>
 */
void Event::updateDist(const double* y, unsigned int len, double t){	
	// update the dist variable using information from y
	lastDist = dist;
	dist = getDist(y, len, t);
}//======================================

/**
 *	\brief Compute the distance from the input state to the event
 *	\param y a state vector representing the current integration state
 *	\param len number of elements in y
 *	\param t non-dimensional time associated with state <tt>y</tt>
 *	
 *	\return the distance
 *	
 *	\throws Exception if the event type associated with this event is not implemented
 *	\throws Exception if the system data pointer has not been initialized via the initialize() function
 */
double Event::getDist(const double *y, unsigned int len, double t) const{
	if(!pSysData)
		throw Exception("Event::getDist: SysData pointer has not been initialized; please call initialize() function.");

	if(len < pSysData->getDynamicsModel()->getCoreStateSize())
		throw Exception("Event::getDist: Input state must contain at least as many elements as the core state size.");

	double d = 0;
	switch(type){
		case Event_tp::SIM_TOF:
		case Event_tp::SIM_COMPTIME:
		case Event_tp::SIM_ERR:
			d = 1;	// Will never use event location to find these events; they are monitored by the simEngine individually
			break;
		case Event_tp::YZ_PLANE:
			d = y[0] - conData[0];
			break;
		case Event_tp::XZ_PLANE:
			d = y[1] - conData[1];
			break;
		case Event_tp::XY_PLANE:
			d = y[2] - conData[2];
			break;
		case Event_tp::CRASH:
		{
			if(len > 2){
				double pos[3] = {0};
				pSysData->getDynamicsModel()->getPrimPos(t, pSysData, static_cast<int>(conData[0]), pos);
				
				double dx = y[0] - pos[0];
				double dy = y[1] - pos[1];
				double dz = y[2] - pos[2];
				d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			}else{
				throw Exception("Event::getDist: input state is too short!");
			}
			break;
		}
		case Event_tp::JC:
		{
			if(len > 5){
				const SysData_cr3bp *crSys = static_cast<const SysData_cr3bp *>(pSysData);
				d = DynamicsModel_cr3bp::getJacobi(y, crSys->getMu()) - conData[0];
			}else{
				throw Exception("Event::getDist: input state is too short!");
			}
			break;
		}
		case Event_tp::APSE:
		{
			if(len > 5){
				int Pix = static_cast<int>(conData[0]);
				d = pSysData->getDynamicsModel()->getRDot(Pix, t, y, pSysData);
			}else{
				throw Exception("Event::getDist: input state is too short!");
			}
			break;
		}
		case Event_tp::DIST:
		{
			if(len > 2){
				double pos[3] = {0};
				pSysData->getDynamicsModel()->getPrimPos(t, pSysData, static_cast<int>(conData[0]), pos);
				
				double dx = y[0] - pos[0];
				double dy = y[1] - pos[1];
				double dz = y[2] - pos[2];
				d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			}else{
				throw Exception("Event::getDist: input state is too short!");
			}
			break;
		}
		case Event_tp::MASS:
			if(len > 6){
				d = y[6] - conData[6];
			}else{
				throw Exception("Event::getDist: input state is too short!");
			}

			break;
		default:
			throw Exception("Event::getDist: Event type not implemented");
	}

	return d;
}//====================================================

/**
 *	\brief Get the direction of propagation.
 *	\details This function should be called after updateDist() to retrieve the most recent
 *	information.
 *
 *	\param tDir direction that time is being propagated (+1 for forward, -1 for negative)
 *	\return positive or negative one to correspond with the sign
 *	\throws Exception if the event type associated with this event is not implemented
 */
int Event::getDir(int tDir) const{
	double d = dist-lastDist;
	return static_cast<int>(d*tDir/std::abs(d*tDir));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	\brief copy the event
 *	\param ev an event
 */
void Event::copyMe(const Event &ev){
	type = ev.type;
	triggerDir = ev.triggerDir;
	triggerCount = ev.triggerCount;
	stopCount = ev.stopCount;
	bStop = ev.bStop;
	dist = ev.dist;
	lastDist = ev.lastDist;
	conType = ev.conType;
	conData = ev.conData;
	paramsIn = ev.paramsIn;
	pSysData = ev.pSysData;	// COPY ADDRESS (ptr)
}//=============================================

/**
 *  \brief Retrieve a human-readable string representing an event type
 *	\return a human-readable string representing an event type
 */
const char* Event::getEventTpStr(Event_tp t){
	switch(t){
		case Event_tp::NONE: return "NONE"; break;
		case Event_tp::SIM_TOF: return "SimEngine Time-of-Flight"; break;
		case Event_tp::SIM_COMPTIME: return "SimEngine Computation Timeout"; break;
		case Event_tp::SIM_ERR: return "SimEngine Error"; break;
		case Event_tp::YZ_PLANE: return "yz-plane"; break;
		case Event_tp::XZ_PLANE: return "xz-plane"; break;
		case Event_tp::XY_PLANE: return "xy-plane"; break;
		case Event_tp::CRASH: return "crash"; break;
		case Event_tp::JC: return "jacobi constant"; break;
		case Event_tp::APSE: return "apse"; break;
		case Event_tp::DIST: return "distance"; break;
		case Event_tp::MASS: return "mass"; break;
		default: return "UNDEFINED!"; break;
	}
}//====================================================

/**
 *	\brief Print out a discription of the event
 */
void Event::printStatus() const{
	printf("Event: Type = %s, Trigger Dir = %d, KillSim = %s\n", getTypeStr(), triggerDir, 
		bStop ? "YES" : "NO");
	printf("  Dist: %e Last Dist: %e\n", dist, lastDist);
}//====================================================

/**
 *  \brief Reset the event to avoid any confusion when a simulation is rerun with the same event
 */
void Event::reset(){
	triggerCount = 0;
	dist = 100000;
	lastDist = 100000;
}//====================================================



}// END of Astrohelion namespace