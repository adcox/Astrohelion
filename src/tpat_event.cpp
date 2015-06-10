/**
 *	@file tpat_event.cpp
 */

#include "tpat_event.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_body_data.hpp"
#include "tpat_calculations.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>
#include <iostream>

using namespace std;

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
		default: 
			printWarn("tpat_event() Creating event with no type!\n");
			throw;
	}
}//================================================


/**
 *	@brief Create an event with custom specifications
 *	
 *	Rather than using the default parameters, this constructor allows you to
 *	create more specialized events. Details below
 *
 *	For the _PLANE crossing events, pass a single double into the <tt>params</tt>
 *	field to specify the location of the plane. For example, choosing an event
 *	with type YZ_PLANE and passing a <tt>param</tt> of 0.5 will create an event
 *	that fires when the trajectory passes through the <tt>x=0.5</tt> plane.
 *
 *	For the CRASH event, pass an array with the first element specifying the 
 *	primary index (0 for P1, 1 for P2, etc.) and the second element specifying
 *	the minimum acceptable distance; distances below this value (non-dim) will
 *	trigger the crash event
 *
 *	@param data a system data object that describes the system this event will occur in
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 *	@param params an array of doubles that give the constructor extra information
 *	
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
		default: 
			printWarn("tpat_event::initEvent() Creating event with no type!\n");
			throw;
	}

	double data[] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN};	// seven empty elements
	switch(type){
		case YZ_PLANE: data[0] = params[0]; break;	// x = 0
		case XZ_PLANE: data[1] = params[0];	break;	// y = 0
		case XY_PLANE: data[2] = params[0]; break;	// z = 0
		case CRASH:
		{
			data[0] = params[0];	// Index of primary
			if(data[0] < sysData->getNumPrimaries()){
				// Get body data, compute crash distance
			    tpat_body_data primData(sysData->getPrimary((int)(data[0])));
			    data[1] = (primData.getRadius() + primData.getMinFlyBy())/sysData->getCharL();
			}else{
				printErr("Cannot access primary #%d for crash event\n", data[0]);
				throw;
			}
			break;
		}
		default: break;	// Do nothing
	}
	conData.insert(conData.begin(), data, data+7);
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
}//=============================================

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
		case NONE: { return "NONE"; break; }
		case YZ_PLANE: { return "yz-plane"; break; }
		case XZ_PLANE: { return "xz-plane"; break; }
		case XY_PLANE: { return "xy-plane"; break; }
		case CRASH: { return "crash"; break; }
		default: { return "UNDEFINED!"; break; }
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
 *	@return the node index for the constrained node. This is hard-coded to be 1, the second
 *	node in a two-node nodeset
 */
int tpat_event::getConNode() const { return 1; }

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
//-----------------------------------------------------
//      Computations, etc.
//-----------------------------------------------------

/**
 *	@brief Determine (roughly) whether or not this event has occured between the
 *	previous trajectory state and the current one.
 *
 *	@param y the current integrated state (6 elements)
 *	@return whether or not the trajectory has passed through this event
 */
bool tpat_event::crossedEvent(double y[6], double t) const{
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
void tpat_event::updateDist(double y[6], double t){	
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
double tpat_event::getDist(double y[6], double t) const{
	double d = 0;
	switch(type){
		case YZ_PLANE: d = conData[0] - y[0]; break;
		case XZ_PLANE: d = conData[1] - y[1]; break;
		case XY_PLANE: d = conData[2] - y[2]; break;
		case CRASH:
		{
			vector<double> primPos;
			switch(sysData->getType()){
				case tpat_sys_data::CR3BP_SYS:
				{
					tpat_cr3bp_sys_data *crSysData = static_cast<tpat_cr3bp_sys_data*>(sysData);
					primPos.assign(6,0);
					primPos[0] = -1*crSysData->getMu();
					primPos[3] = 1 - crSysData->getMu();
					break;
				}
				case tpat_sys_data::BCR4BPR_SYS:
				{
					tpat_bcr4bpr_sys_data *bcSysData = static_cast<tpat_bcr4bpr_sys_data*>(sysData);
					primPos.assign(9, 0);
					bcr4bpr_getPrimaryPos(t, *bcSysData, &(primPos[0]));
					break;
				}
				default: 
					printErr("tpat_event::getDist() Unsupported system type for crash\n");
					throw;
			}

			int Pix = (int)(conData[0]);
			double dx = y[0] - primPos[Pix*3 + 0];
			double dy = y[1] - primPos[Pix*3 + 1];
			double dz = y[2] - primPos[Pix*3 + 2];
			d = sqrt(dx*dx + dy*dy + dz*dz) - conData[1];
			break;
		}
		default:
			printErr("tpat_event::getDist() Event type not implemented: %s\n", getTypeStr());
			throw;
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
int tpat_event::getDir(double y[6], double t) const{
	double d = 0;
	double dt = t - theTime;
	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case YZ_PLANE: d = y[0] - state[0]; break;
		case XZ_PLANE: d = y[1] - state[1]; break;
		case XY_PLANE: d = y[2] - state[2]; break;
		case CRASH: d = dist - lastDist; break;
		default: 
			printErr("tpat_event::getDir() Event type not implemented: %s\n", getTypeStr());
			throw;
	}

	return (int)(d*dt/abs(d*dt));
}//==============================================

/**
 *	@brief Print out a discription of the event
 */
void tpat_event::printStatus() const{
	printf("Event: Type = %s, Trigger Dir = %d, KillSim = %s\n", getTypeStr(), triggerDir, 
		stop ? "YES" : "NO");
	printf("  Dist: %e Last Dist: %e\n", dist, lastDist);
}//======================================


