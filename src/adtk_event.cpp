/**
 *	@file adtk_event.cpp
 */

#include "adtk_event.hpp"
#include "adtk_utilities.hpp"
 
#include <cmath>
#include <iostream>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/** @brief Default, do-nothing constructor */
adtk_event::adtk_event(){}

/**
 *	@brief Create an event
 *	@param t the event type
 *	@param dir direction (+/-/both) the event will trigger on. +1 indices (+)
 *	direction, -1 (-) direction, and 0 both directions.
 *	@param willStop whether or not this event should stop the integration
 */
adtk_event::adtk_event(event_t t, int dir, bool willStop){
	type = t;
	direction = dir;
	stop = willStop;

	// Create constraint data based on the type
	conType = adtk_constraint::STATE;
	double data[] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN};	// seven empty elements
	switch(type){
		case YZ_PLANE: data[0] = 0; break;	// x = 0
		case XZ_PLANE: data[1] = 0;	break;	// y = 0
		case XY_PLANE: data[2] = 0; break;	// z = 0
		default: break;	// Do nothing
	}
	conData.insert(conData.begin(), data, data+7);
}//================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the direction for this event; -1 for negative, +1
 *	for positive, 0 for both/either
 */
int adtk_event::getDir() const { return direction; }

/**
 *	@return the event type
 */
adtk_event::event_t adtk_event::getType() const { return type; }

/**
 *	@return a human-readable string representing the event type
 */
const char* adtk_event::getTypeStr(){
	switch(type){
		case NONE: { return "NONE"; break; }
		case YZ_PLANE: { return "yz-plane"; break; }
		case XZ_PLANE: { return "xz-plane"; break; }
		case XY_PLANE: { return "xy-plane"; break; }
		default: { return "UNDEFINED!"; break; }
	}
}//========================================

/**
 *	@return the time associated with this event
 */
double adtk_event::getTime() const { return theTime; }

/**
 *	@return a pointer to the state vector object; useful for in-place reading or writing
 */
std::vector<double>* adtk_event::getState() { return &state; }

/**
 *	@return whether or not this event will stop the integration
 */
bool adtk_event::stopOnEvent() const { return stop; }

/**
 *	@return the type of constraint this event will use to target the exact event occurence
 */
adtk_constraint::constraint_t adtk_event::getConType() const { return conType; }

/**
 *	@return the constraint data used to target this exact event
 */
std::vector<double> adtk_event::getConData() const { return conData; }

/**
 *	@return the node index for the constrained node. This is hard-coded to be 1, the second
 *	node in a two-node nodeset
 */
int adtk_event::getConNode() const { return 1; }

/**
 *	@brief Set the direction for this event
 *	@param d the direction: +1 for positive, -1 for negative, 0 for both/either
 */
void adtk_event::setDir(int d){ direction = d; }

//-----------------------------------------------------
//      Computations, etc.
//-----------------------------------------------------

/**
 *	@brief Determine (roughly) whether or not this event has occured between the
 *	previous trajectory state and the current one.
 *
 *	@param y the current integrated state (6 elements)
 *	@return whehter or not the trajectory has passed through this event
 */
bool adtk_event::crossedEvent(double y[6]) const{
	double newDist = getDist(y);

	// See if we have crossed (in pos. or neg. direction)
	if(newDist*dist < 0){ // have different signs
		if(direction == 0){
			return true;
		}else{
			return direction == getDir(y);
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
 */
void adtk_event::updateDist(double y[6]){
	// update the dist variable using information from y
	dist = getDist(y);

	// Save the state from y for later comparison
	state.clear();
	state.insert(state.begin(), y, y+6);
}//======================================

/**
 *	@brief Compute the distance from the input state to the event
 *	@param y a 6-element state vector representing the current integration state
 *	@return the distance
 */
double adtk_event::getDist(double y[6]) const{
	double d = 0;
	switch(type){
		case YZ_PLANE: d = 0 - y[0]; break;
		case XZ_PLANE: d = 0 - y[1]; break;
		case XY_PLANE: d = 0 - y[2]; break;
		default:
			printErr("Event type not implemented!\n");
			throw;
	}

	return d;
}//=====================================

/**
 *	@brief Get the direction of propagation for the event by comparing an input state <tt>y</tt>
 *	to the state stored in the <tt>state</tt> variable (which was updated last iteration)
 *
 *	@param y a 6-element state vector
 *	@return positive or negative one to correspond with the sign
 */
int adtk_event::getDir(double y[6]) const{
	double d = 0;
	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case YZ_PLANE: d = y[0] - state[0]; break;
		case XZ_PLANE: d = y[1] - state[1]; break;
		case XY_PLANE: d = y[2] - state[2]; break;
		default: 
			printErr("Event type not implemented!\n");
			throw;
	}

	return (int)(d/abs(d));
}//==============================================

