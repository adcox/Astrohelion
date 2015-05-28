/**
 *	@file adtk_event.cpp
 */

#include "adtk_event.hpp"

#include <cmath>
#include <iostream>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_event::adtk_event(){}

adtk_event::adtk_event(event_t t, int dir, bool willStop, int ix){
	type = t;
	direction = dir;
	stop = willStop;
	index = ix;

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

int adtk_event::getDir() const { return direction; }

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

double adtk_event::getTime() const { return theTime; }

std::vector<double>* adtk_event::getState() { return &state; }

int adtk_event::getIndex() const { return index; }

bool adtk_event::stopOnEvent() const { return stop; }

adtk_constraint::constraint_t adtk_event::getConType() const { return conType; }

std::vector<double> adtk_event::getConData() const { return conData; }

int adtk_event::getConNode() const { return 1; }

void adtk_event::setDir(int d){ direction = d; }

void adtk_event::setType(event_t t){ type = t; }

void adtk_event::setIndex(int ix){ index = ix; }

//-----------------------------------------------------
//      Computations, etc.
//-----------------------------------------------------

bool adtk_event::crossedEvent(double y[6]){
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
 *	Update the distance variable, which will later be compared in 
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
 *	Compute the distance from the input state to the event
 *	@param y a 6-element state vector representing the current integration state
 *	@return the distance
 */
double adtk_event::getDist(double y[6]){
	double d = 0;
	switch(type){
		case YZ_PLANE: d = 0 - y[0]; break;
		case XZ_PLANE: d = 0 - y[1]; break;
		case XY_PLANE: d = 0 - y[2]; break;
		default:
			fprintf(stderr, "Event type not implemented!\n");
			throw;
	}

	return d;
}//=====================================

/**
 *	Get the direction of propagation for the event by comparing an input state <tt>y</tt>
 *	to the state stored in the <tt>state</tt> variable (which was updated last iteration)
 *
 *	@param y a 6-element state vector
 *	@return positive or negative one to correspond with the sign
 */
int adtk_event::getDir(double y[6]){
	double d = 0;
	// Compute distance from old point (in state) to new point (in y)
	switch(type){
		case YZ_PLANE: d = y[0] - state[0]; break;
		case XZ_PLANE: d = y[1] - state[1]; break;
		case XY_PLANE: d = y[2] - state[2]; break;
		default: 
			fprintf(stderr, "Event type not implemented!\n");
			throw;
	}

	return (int)(d/abs(d));
}//==============================================

