/**
 *  @file tpat_model.cpp
 *	@brief Defines behavior for a dynamic model
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_model.hpp"

#include "tpat_correction_engine.hpp"
#include "tpat_event.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset.hpp"
#include "tpat_traj.hpp"
#include "tpat_sys_data.hpp"

#include <algorithm>
#include <cmath>

/**
 *	@brief Default constructor
 *	@param type the model type
 */
tpat_model::tpat_model(dynamicModel_t type){
	modelType = type;
}//===========================================

/**
 *	@brief Copy constructor
 */
tpat_model::tpat_model(const tpat_model &m){
	copyMe(m);
}//===========================================

/**
 *	@brief Deconstructor
 */
tpat_model::~tpat_model(){
	allowedCons.clear();
	allowedEvents.clear();
}//==============================

/**
 *	@brief Copy Operator
 */	
tpat_model& tpat_model::operator =(const tpat_model &m){
	copyMe(m);
	return *this;
}//============================================

/**
 *	@brief Copies all data from a dynamic model into this one
 *	@param m another dynamic model
 */
void tpat_model::copyMe(const tpat_model &m){
	modelType = m.modelType;
	coreStates = m.coreStates;
	stmStates = m.stmStates;
	extraStates = m.extraStates;
}//============================================

/**
 *	@brief Retrieve the number of core states
 *	@return the number of core states
 */
int tpat_model::getCoreStateSize() const { return coreStates; }

/**
 *	@brief Retrieve the number of STM elements stored in the state vector
 *	@return the number of STM elements stored in the state vector
 */
int tpat_model::getSTMStateSize() const { return stmStates; }

/**
 *	@brief Retrieve the number of extra states stored after the core states and STM elements
 *	@return the number of extra states stored after the core states and STM elements
 */
int tpat_model::getExtraStateSize() const { return extraStates; }

/**
 *	@brief Determine whether the specified constraint type is supported in this model
 *	@return whether or not the specified constraint type is supported in this model
 */
bool tpat_model::supportsCon(tpat_constraint::constraint_t type) const{
	return std::find(allowedCons.begin(), allowedCons.end(), type) != allowedCons.end();
}//===================================================

/**
 *	@brief Determine whether the specified event type is supported in this model
 *	@return whether or not the specified event type is supported in this model
 */
bool tpat_model::supportsEvent(tpat_event::event_t type) const{
	return std::find(allowedEvents.begin(), allowedEvents.end(), type) != allowedEvents.end();
}//===================================================

/**
 *	@brief Initialize the corrector's design vector with position and velocity states,
 *	and times-of-flight.
 *
 *	Derived models may replace this function or call it and then append more design 
 * 	variables.
 *
 *	@param it a pointer to the corrector's iteration data structure
 *	@param set a pointer to the nodeset being corrected
 */
void tpat_model::corrector_initDesignVec(iterationData *it, tpat_nodeset *set){
	// Create the initial state vector
	it->X.clear();

	// Copy in the position and velocity states for each node
	for(int n = 0; n < set->getNumNodes(); n++){
		std::vector<double> state = set->getNode(n).getPosVelState();
		it->X.insert(it->X.end(), state.begin(), state.end());
	}

	// Append the TOF for each node (except the last one, which isn't propagated)
	if(it->varTime){
		for(int n = 0; n < set->getNumNodes()-1; n++){
			it->X.insert(it->X.end(), set->getNode(n).getTOF());
		}
	}
}//============================================================

/**
 *	@brief Create continuity constraints for the correction algorithm; this function
 *	creates position and velocity constraints.
 *
 *	Defived models may replace this function or call it and then append more constraints
 *	for other variables that may be continuous, such as time (non-autonomous systems)
 *	or mass.
 *
 *	@param it a pointer to the corrector's iteration data structure
 *	@param set a pointer to the nodeset being corrected
 */	
void tpat_model::corrector_createContCons(iterationData *it, tpat_nodeset *set){
	// Create position and velocity constraints
    for(int n = 1; n < set->getNumNodes(); n++){
        // Get a vector specifying which velocity states are continuous from the nodeset
        std::vector<bool> velCon = set->getNode(n).getVelCon();
        // Force all positions to be continuous, use specified velocity continuities
        std::vector<double> contStates {1,1,1,(double)velCon[0],(double)velCon[1],(double)velCon[2]};
        // Create a constraint
        tpat_constraint continuity(tpat_constraint::CONT_PV, n, contStates);
        // Save constraint to total constraint vector
        it->allCons.push_back(continuity);
    }

	// Next, create other continuity constraints via the CONT_EX constraint type
	// These may include time continuity for non-autonomous systems, or mass continuity
}//============================================================

/**
 *	@brief Retrieve the initial conditions for a segment that the correction
 *	engine will integrate.
 *
 *	Derived models may replace this function to change how the initial conditions
 *	are chosen from the design vector.
 *
 *	@param it a pointer to the corrector's iteration data structure
 *	@param set a pointer to the nodeset being corrected
 *	@param n the index of the node that serves as the initial state
 *	@param ic a pointer to a 6-element initial state array
 *	@param t0 a pointer to a double representing the initial time (epoch)
 *	@param tof a pointer to a double the time-of-flight on the segment.
 */
void tpat_model::corrector_getSimICs(iterationData *it, tpat_nodeset *set, int n,
	double *ic, double *t0, double *tof){
	
	std::copy(it->X.begin()+6*n, it->X.begin()+6*(n+1), ic);
	*tof = it->varTime ? it->X[6*it->numNodes+n] : set->getTOF(n);
	*t0 = 0;
}//============================================================

/**
 *	@brief Compute constraint function and partial derivative values for a constraint
 *	
 *	This function provides a framework for most of the constraints available to all models.
 *	Given an iteration data object, a constraint, and the index of that constraint, we can
 *	compute the value of the constraint function(s) associated with the constraint as well
 *	as the partial derivatives that relate the constraint functions to design variables.
 *	Although this model defines behavior for many of the constraints, it is up to the model
 *	developer to ensure any new constraints are implemented in a derived model, and any
 *	existing constraints are treated fully and appropriately by the new, derived model.
 *
 *	@param it a pointer to the corrector's iteration data structure
 *	@param con the constraint being applied
 *	@param c the index of the constraint within the total constraint vector (which is, in
 *	turn, stored in the iteration data)
 */	
void tpat_model::corrector_applyConstraint(iterationData *it, tpat_constraint con, int c){

	int row0 = it->conRows[c];

	switch(con.getType()){
		case tpat_constraint::CONT_PV:
			// Apply position-velocity continuity constraints
			corrector_targetPosVelCons(it, con, row0);
			break;
		case tpat_constraint::CONT_EX:
			// Apply extra continuity constraints
			corrector_targetExContCons(it, con, row0);
			break;
		case tpat_constraint::STATE:
			corrector_targetState(it, con, row0);
			break;
		case tpat_constraint::MATCH_ALL:
			corrector_targetMatchAll(it, con, row0);
			break;
		case tpat_constraint::MATCH_CUST:
			corrector_targetMatchCust(it, con, row0);
			break;
		case tpat_constraint::MAX_DIST:
		case tpat_constraint::MIN_DIST:
		case tpat_constraint::DIST:
			corrector_targetDist(it, con, c);
			break;
		case tpat_constraint::DELTA_V:
		case tpat_constraint::MAX_DELTA_V:
			corrector_targetDeltaV(it, con, c);
			break;
		case tpat_constraint::TOF:
			corrector_targetTOF(it, con, row0);
			break;
		case tpat_constraint::APSE:
			corrector_targetApse(it, con, row0);
			break;	
		default: break;
	}
}//=========================================================

/**
 *	@brief Compute position and velocity constraint values and partial derivatives
 *
 *	This function computes and stores the default position continuity constraints as well
 *	as velocity constraints for all nodes marked continuous in velocity. The delta-Vs
 *	between arc segments and node states are recorded, and the partial derivatives of each
 *	node with respect to other nodes and integration time are all computed
 *	and placed in the appropriate spots in the Jacobian matrix.
 *
 *	Derived models may replace this function.
 *
 *	@param it a pointer to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void tpat_model::corrector_targetPosVelCons(iterationData* it, tpat_constraint con, int row0){
	int n = con.getNode();
	if(n == 0)
		throw tpat_exception("tpat_model::corrector_targetPosVelCons: Cannot constraint node 0 to be continuous with node -1");

	// Get info about the arc that was integrated to reach node n
	std::vector<double> lastState = it->allSegs.at(n-1).getState(-1);
	std::vector<double> lastAccel = it->allSegs.at(n-1).getAccel(-1);
	tpat_matrix stm = it->allSegs.at(n-1).getSTM(-1);
	std::vector<double> conData = con.getData();
	
	// Loop through conData
	for(size_t s = 0; s < conData.size(); s++){
		if(!isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			it->FX[row0+s] = lastState[s] - it->X[6*n+s];

			// Loop through all design variables for this node (6) and compute partials of F w.r.t. x
			
			for(size_t x = 0; x < 6; x++){
				// put STM elements into DF matrix
				it->DF[it->totalFree*(row0+s) + 6*(n-1)+x] = stm.at(s,x);
				// Negative identity matrix
				if(s == x)
					it->DF[it->totalFree*(row0+s) + 6*n+x] = -1;
			}

			// Compute partials of F w.r.t. times-of-flight
			// Columns of DF based on time constraints
			if(it->varTime){
				// Column of state derivatives: [vel; accel]
				if(s < 3)
					it->DF[it->totalFree*(row0+s) + 6*it->numNodes+n-1] = lastState[s+3];
				else{					
					it->DF[it->totalFree*(row0+s) + 6*it->numNodes+n-1] = lastAccel[s-3];
				}
			}
		}
	}
}//=========================================================

/**
 *	@brief Computes continuity constraints for constraints with the <tt>CONT_EX</tt> type.
 *
 *	In this base model, no behavior is defined for extra constraints. It is intended to enforce
 *	continuity constraints like epoch (time) continuity, mass continuity, etc.
 *
 *	@param it a pointer to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void tpat_model::corrector_targetExContCons(iterationData *it, tpat_constraint con, int row0){
	// Do absoluately nothing
	(void)it;
	(void)con;
	(void)row0;
}

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>STATE</tt>.
 *
 *	This method *should* provide full state constraining for any model; the STM and identity 
 *	matrices are used to relate node states and integrated states.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_model::corrector_targetState(iterationData* it, tpat_constraint con, int row0){
	std::vector<double> conData = con.getData();
	int n = con.getNode();
	// Allow user to constrain all 6 states
	
	int count = 0; 	// Count # rows since some may be skipped (NAN)
	for(int s = 0; s < ((int)con.getData().size()); s++){
		if(!isnan(conData[s])){
			if(s < 6){
				it->FX[row0+count] = it->X[6*n+s] - conData[s];
				it->DF[it->totalFree*(row0 + count) + 6*n + s] = 1;
				count++;
			}else{
				throw tpat_exception("State constraints must have <= 6 elements");
			}
		}
	}
}//=================================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>MATCH_ALL</tt>
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are applied
 *	to the Jacobian matrix.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_model::corrector_targetMatchAll(iterationData* it, tpat_constraint con, int row0){
	// Only allow matching 6 states, not TOF (state 7)
	int n = con.getNode();
	int cn = con.getData()[0];
	for(int row = 0; row < 6; row++){
		// Constrain the states of THIS node to be equal to the node 
		// with index stored in conData[0]
		it->FX[row0+row] = it->X[6*n+row] - it->X[6*cn+row];

		// Partial of this constraint wrt THIS node = I
		it->DF[it->totalFree*(row0 + row) + 6*n + row] = 1;

		// Partial of this constraint wrt other node = -I
		it->DF[it->totalFree*(row0 + row) + 6*cn+row] = -1;
	}
}//=============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>MATCH_CUST</tt>
 *
 *	This method *should* provide full functionality for any model; Only 1's and 0's are applied
 *	to the Jacobian matrix.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_model::corrector_targetMatchCust(iterationData* it, tpat_constraint con, int row0){
	std::vector<double> conData = con.getData();
	int n = con.getNode();
	int count = 0;
	// Only allow matching 6 states, not epoch time (state 7)
	for(int s = 0; s < 6; s++){
		if(!isnan(conData[s])){
			int cn = conData[0];
			it->FX[row0 + count] = it->X[6*n+s] - it->X[6*cn+s];

			// partial of this constraint wrt THIS node = 1
			it->DF[it->totalFree*(row0 + count) + 6*n+s] = 1;

			// partial of this constraint wrt other node = -1
			it->DF[it->totalFree*(row0 + count) + 6*cn+s] = -1;

			count++;
		}
	}
}//===============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>DIST</tt>, 
 *	<tt>MIN_DIST</tt>, or <tt>MAX_DIST</tt>
 *
 *	This method *should* provide full functionality for any model; It calls the getPrimPos() 
 *	functions, which all models define and uses dynamic-independent computations to populate
 *	the constraint vector and Jacobian matrix.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param c the index of this constraint in the constraint vector object
 */
void tpat_model::corrector_targetDist(iterationData* it, tpat_constraint con, int c){

	std::vector<double> conData = con.getData();
	int n = con.getNode();
	int Pix = (int)(conData[0]);	// index of primary
	int row0 = it->conRows[c];
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time
	tpat_sys_data *sysData = it->sysData;

	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, sysData);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it->X[6*n+0] - primPos[Pix*3+0];
	double dy = it->X[6*n+1] - primPos[Pix*3+1];
	double dz = it->X[6*n+2] - primPos[Pix*3+2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

	// Compute difference between desired distance and true distance
	it->FX[row0] = h - conData[1];

	// Partials with respect to node position states
	it->DF[it->totalFree*row0 + 6*n + 0] = dx/h;
	it->DF[it->totalFree*row0 + 6*n + 1] = dy/h;
	it->DF[it->totalFree*row0 + 6*n + 2] = dz/h;

	// Extra stuff for inequality constraints
	if(con.getType() == tpat_constraint::MIN_DIST || 
		con.getType() == tpat_constraint::MAX_DIST ){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
			it->slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
		int sign = con.getType() == tpat_constraint::MAX_DIST ? 1 : -1;

		// Subtract squared slack variable from constraint
		it->FX[row0] += sign*it->X[slackCol]*it->X[slackCol];

		// Partial with respect to slack variable
		it->DF[it->totalFree*row0 + slackCol] = sign*2*it->X[slackCol];
	}
}// End of targetDist() =========================================

/**
 *	@brief Compute partials and constraints for all nodes constrained with <tt>DELTA_V</tt> or
 *	<tt>MIN_DELTA_V</tt>
 *
 *	Because the delta-V constraint applies to the entire trajectory, the constraint function values
 *	and partial derivatives must be computed for each node along the trajectory. This function
 *	takes care of all of them at once.
 *
 *	This function computes the partials and constraints for an autonomous system; Non-autonomous
 *	dynamic models should implement a similar copy of this function that computes partials
 *	with respect to their additional design variables.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 *	@param c the index of the first row for this constraint
 */
void tpat_model::corrector_targetDeltaV(iterationData* it, tpat_constraint con, int c){

	int row0 = it->conRows[c];

	// Compute total deltaV magnitude
	double totalDV = 0;
	for(int n = 0; n < it->numNodes-1; n++){
		// compute magnitude of DV between node n and n+1
		double dvMag = sqrt(it->deltaVs[n*3]*it->deltaVs[n*3] +
		it->deltaVs[n*3+1]*it->deltaVs[n*3+1] + 
		it->deltaVs[n*3+2]*it->deltaVs[n*3+2]);

		totalDV += dvMag;

		// Compute parial w.r.t. node n+1 (where velocity is discontinuous)
		double dFdq_ndf_data[] = {0, 0, 0, -1*it->deltaVs[n*3]/dvMag, 
			-1*it->deltaVs[n*3+1]/dvMag, -1*it->deltaVs[n*3+2]/dvMag};
		tpat_matrix dFdq_n2(1, 6, dFdq_ndf_data);

		// Get info about the final state/accel of the integrated segment
		tpat_matrix stm = it->allSegs[n].getSTM(-1);
		std::vector<double> state_dot_data;
		std::vector<double> lastState = it->allSegs[n].getState(-1);
		std::vector<double> lastAccel = it->allSegs[n].getAccel(-1);
		state_dot_data.insert(state_dot_data.end(), lastState.begin()+3, lastState.begin()+6);
		state_dot_data.insert(state_dot_data.end(), lastAccel.begin(), lastAccel.end());

		// Partial w.r.t. integrated path (newSeg) from node n
		tpat_matrix dFdq_nf = -1*dFdq_n2*stm;

		// Compute partial w.r.t. integration time n
		tpat_matrix state_dot(6, 1, &(state_dot_data[0]));
		tpat_matrix dFdt_n = -1*dFdq_n2 * state_dot;

		for(int i = 0; i < 6; i++){
			it->DF[it->totalFree*row0 + 6*(n+1) + i] = dFdq_n2.at(0, i);
			it->DF[it->totalFree*row0 + 6*n + i] = dFdq_nf.at(0, i);
		}
	}
	
	// Copute the difference between the actual deltaV and the desired deltaV
	it->FX[row0] = totalDV - con.getData()[0];

	if(con.getType() == tpat_constraint::DELTA_V){
		it->FX[row0] -= con.getData()[0];
	}else if(con.getType() == tpat_constraint::MAX_DELTA_V){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(),
			it->slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());

		/* FIRST ITERATION ONLY (before the targeter computes a value for beta)
		Set the slack variable such that the constraint will evaluate to zero if 
		the actual deltaV is less than the required deltaV */
		if(it->count == 0){
			it->X[slackCol] = sqrt(std::abs(con.getData()[0] - it->FX[row0]));
		}
		it->FX[row0] -= con.getData()[0] - it->X[slackCol]*it->X[slackCol];
		it->DF[it->totalFree*row0 + slackCol] = 2*it->X[slackCol];
	}
}//==============================================

/**
 *	@brief Compute partials and constraint function values for time-of-flight constraints
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are
 *	used to relate TOFs.
 */
void tpat_model::corrector_targetTOF(iterationData *it, tpat_constraint con, int row0){
	// Sum all TOF for total, set partials w.r.t. integration times equal to one
	for(int i = 0; i < it->numNodes-1; i++){
		it->FX[row0] += it->X[6*it->numNodes+i];
		it->DF[it->totalFree*row0 + 6*it->numNodes+i] = 1;
	}
	
	// subtract the desired TOF from the constraint to finish its computation
	it->FX[row0] -= con.getData()[0];
}//===============================================

/**
 *	@brief Compute partials and constraint function values for apse constraints
 *
 *	This method *should* provide full functionality for any autonomous model. Non-
 *	autonomous models will need to modify the function to account for epoch time
 */
void tpat_model::corrector_targetApse(iterationData *it, tpat_constraint con, int row0){
	std::vector<double> conData = con.getData();
	int n = con.getNode();
	int Pix = (int)(conData[0]);	// index of primary
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time
	tpat_sys_data *sysData = it->sysData;

	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, sysData);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it->X[6*n+0] - primPos[Pix*3+0];
	double dy = it->X[6*n+1] - primPos[Pix*3+1];
	double dz = it->X[6*n+2] - primPos[Pix*3+2];

	// Constraint function: r_dot = 0 
	it->FX[row0] = dx*(it->X[6*n+3]) + dy*(it->X[6*n+4]) + dz*(it->X[6*n+5]);

	// Partials of F w.r.t. node state
	it->DF[it->totalFree*row0 + 6*n+0] = it->X[6*n+3];
	it->DF[it->totalFree*row0 + 6*n+1] = it->X[6*n+4];
	it->DF[it->totalFree*row0 + 6*n+2] = it->X[6*n+5];
	it->DF[it->totalFree*row0 + 6*n+3] = dx;
	it->DF[it->totalFree*row0 + 6*n+4] = dy;
	it->DF[it->totalFree*row0 + 6*n+5] = dz;
}//===============================================




