/**
 *  \file DynamicsModel.cpp
 *	\brief Defines behavior for a dynamic model
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "DynamicsModel.hpp"

#include "MultShootEngine.hpp"
#include "Event.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Nodeset.hpp"
#include "Traj.hpp"
#include "SysData.hpp"
#include "Utilities.hpp"
 
#include <algorithm>
#include <cmath>

namespace astrohelion{
/**
 *	\brief Default constructor
 *	\param type the model type
 */
DynamicsModel::DynamicsModel(DynamicsModel_tp type){
	modelType = type;
}//===========================================

/**
 *	\brief Copy constructor
 */
DynamicsModel::DynamicsModel(const DynamicsModel &m){
	copyMe(m);
}//===========================================

/**
 *	\brief Deconstructor
 */
DynamicsModel::~DynamicsModel(){}

/**
 *	\brief Copy Operator
 */	
DynamicsModel& DynamicsModel::operator =(const DynamicsModel &m){
	copyMe(m);
	return *this;
}//============================================

/**
 *	\brief Copies all data from a dynamic model into this one
 *	\param m another dynamic model
 */
void DynamicsModel::copyMe(const DynamicsModel &m){
	modelType = m.modelType;
	coreStates = m.coreStates;
	extraStates = m.extraStates;
}//============================================

/**
 *	\brief Retrieve the number of core states
 *	\return the number of core states
 */
unsigned int DynamicsModel::getCoreStateSize() const { return coreStates; }

/**
 *	\brief Retrieve the number of extra states stored after the core states and STM elements
 *	\return the number of extra states stored after the core states and STM elements
 */
unsigned int DynamicsModel::getExtraStateSize() const { return extraStates; }

/**
 *	\brief Determine whether the specified constraint type is supported in this model
 *	\return whether or not the specified constraint type is supported in this model
 */
bool DynamicsModel::supportsCon(Constraint_tp type) const{
	return std::find(allowedCons.begin(), allowedCons.end(), type) != allowedCons.end();
}//===================================================

/**
 *	\brief Determine whether the specified event type is supported in this model
 *	\return whether or not the specified event type is supported in this model
 */
bool DynamicsModel::supportsEvent(Event_tp type) const{
	return std::find(allowedEvents.begin(), allowedEvents.end(), type) != allowedEvents.end();
}//===================================================

/**
 *  \brief Determine the time derivative of the magnitude of a vector from a primary to
 *  the body of interest, non-dimensional units
 *  \details This derivation assumes the primary of interest is fixed in the frame the
 *  spacecraft coordinates are expressed in. If this is not the case (i.e., the primary
 *  moves in the working frame), this function will need to be overridden
 * 
 *  \param Pix Index of the primary
 *  \param t Non-dimensional epoch to compute r-dot at
 *  \param state six-element non-dimensional spacecraft state
 *  \param sys system data object
 *  \return Time derivative of the magnitude of a vector from a primary to the body
 *  of interest, non-dimensional velocity units
 */
double DynamicsModel::getRDot(int Pix, double t, const double *state, const SysData *sys) const{
	
	std::vector<double> primPos = getPrimPos(t, sys);
    double dx = state[0] - primPos[3*Pix+0];
    double dy = state[1] - primPos[3*Pix+1];
    double dz = state[2] - primPos[3*Pix+2];

    std::vector<double> primVel = getPrimVel(t, sys);
    double num = dx*(state[3] - primVel[3*Pix+0]) + dy*(state[4] - primVel[3*Pix+1])+ dz*(state[5] - primVel[3*Pix+2]);

    return num/sqrt(dx*dx + dy*dy + dz*dz);
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Create default events for a simulation run
 *  \details These events are intended to prevent numerical issues, e.g., to avoid singularities.
 * 
 *  \param pSys pointer to system data object
 *  \return A vector of events to use in the simulation
 */
std::vector<Event> DynamicsModel::sim_makeDefaultEvents(const SysData *pSys) const{
	std::vector<Event> events;
	// Create events that prevent the s/c from flying through planets
	for(int p = 0; p < pSys->getNumPrimaries(); p++){
        // Put primary index # into an array, create event
        std::vector<double> Pix {static_cast<double>(p)};
        // Event e(Event_tp::CRASH, 0, true, Pix);
        events.push_back(Event(Event_tp::CRASH, 0, true, Pix));
        // events.push_back(e);
    }
    return events;
}//==================================================

/**
 *  \brief Save data from the simulation to a trajectory object
 *  \details [long description]
 * 
 *  \param y pointer to state data array
 *  \param t time
 *  \param traj pointer to trajectory in which the data is stored
 *  \param params Structure that contains parameters used in the integration
 */
void DynamicsModel::sim_saveIntegratedData(const double *y, double t, Traj* traj, EOM_ParamStruct *params) const{
    int id = traj->addNode(Node(y, coreStates, t));
	if(id > 0){
		double tof = t - traj->getNodeRef_const(id-1).getEpoch();
		traj->addSeg(Segment(id-1, id, tof, y+coreStates, coreStates*coreStates, params->ctrlLawID));
	}
}//====================================================

/**
 *  \brief Create a node on the Trajectory
 *  \details [long description]
 * 
 *  \param y pointer to state data array
 *  \param t time
 *  \param traj pointer to trajectory in which the data is stored
 *  \param params Structure that contains parameters used in the integration
 *  \param evt Event type that occured at the Node
 *  \return ID of the node once it is added to the Trajectory
 *  
 *  \todo Incorporate the Event-saving part of this function
 */
int DynamicsModel::sim_addNode(Node &node, const double *y, double t, Traj* traj, EOM_ParamStruct *params, Event_tp tp) const{
	(void) y;
	(void) t;
	(void) params;
	(void) tp;
	return traj->addNode(node);
}//====================================================

int DynamicsModel::sim_addSeg(Segment &seg, const double *y, double t, Traj* traj, EOM_ParamStruct *params) const{
	(void) y;
	(void) t;
	(void) params;
	return traj->addSeg(seg);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *	\brief Initialize the corrector's design vector with position and velocity states,
 *	and times-of-flight.
 *
 *	Derived models may replace this function or call it and then append more design 
 * 	variables.
 *
 *	\param it a pointer to the corrector's iteration data structure
 *	\param set a pointer to the nodeset being corrected
 *	
 *	\throws Exception if equalArcTime is set to true (within the correction engine) and
 *	the nodeset to be corrected includes segments that are propagated in both forward and
 *	reverse time.
 */
void DynamicsModel::multShoot_initDesignVec(MultShootData *it, const Nodeset *set) const{
	// Create the initial state vector
	it->X.clear();
	it->freeVarMap.clear();

	// Copy in the state vector for each node
	int rowNum = 0;
	for(int n = 0; n < set->getNumNodes(); n++){
		std::vector<double> state = set->getStateByIx(n);
		rowNum = it->X.size();
		it->X.insert(it->X.end(), state.begin(), state.end());
		MSVarMap_Key key(MSVar_tp::STATE, set->getNodeRefByIx_const(n).getID());
		it->freeVarMap[key] = MSVarMap_Obj(key, rowNum, state.size());
	}

	if(it->bVarTime){		
		if(it->bEqualArcTime){
			// Make sure all times-of-flight have the same sign
			for(int s = 1; s < set->getNumSegs(); s++){
				if(set->getTOFByIx(s) * set->getTOFByIx(s-1) < 0)
					throw Exception("DynamicsModel::multShoot_initDesignVec: EqualArcTime is ON and times-of-flight have different signs... cannot proceed");
			}

			// Append the total TOF for the arc
			MSVarMap_Key key(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
			it->freeVarMap[key] = MSVarMap_Obj(key, static_cast<int>(it->X.size()));
			it->X.insert(it->X.end(), set->getTotalTOF());
		}else{
			// Append the TOF for each segment
			for(int s = 0; s < set->getNumSegs(); s++){
				MSVarMap_Key key(MSVar_tp::TOF, set->getSegRefByIx_const(s).getID());
				it->freeVarMap[key] = MSVarMap_Obj(key, static_cast<int>(it->X.size()));
				it->X.insert(it->X.end(), set->getTOFByIx(s));
			}
		}
	}
}//============================================================

/**
 *	\brief Create continuity constraints for the correction algorithm; this function
 *	creates position and velocity constraints.
 *
 *	Derived models may replace this function or call it and then append more constraints
 *	for other variables that may be continuous, such as time (non-autonomous systems)
 *	or mass. This function assumes the state elements at indices 3, 4, and 5 represent
 *	velocity and their continuity is governed by the boolean flags stored in a segment.
 *
 *	\param it a pointer to the corrector's iteration data structure
 *	\param set a pointer to the nodeset being corrected
 */	
void DynamicsModel::multShoot_createContCons(MultShootData *it, const Nodeset *set) const{
	// Create position and velocity constraints
	for(int s = 0; s < set->getNumSegs(); s++){
		// Force all positions to be continuous
		std::vector<double> contStates(coreStates, 1);
		if(set->getSegRefByIx_const(s).getTerminus() != Linkable::INVALID_ID){	
			// Get a vector specifying which velocity states are continuous
			std::vector<bool> velCon = set->getSegRefByIx_const(s).getVelCon();
			// If not continuous, put NAN into the constraint data; else unity
			contStates[3] = velCon[0] ? 1 : NAN;
			contStates[4] = velCon[1] ? 1 : NAN;
			contStates[5] = velCon[2] ? 1 : NAN;
			
			// Create a constraint
			Constraint con(Constraint_tp::CONT_PV, set->getSegRefByIx_const(s).getID(), contStates);

			// Save constraint to constraint vector
			it->allCons.push_back(con);
		}
	}

	// Next, create other continuity constraints via the CONT_EX constraint type
	// These may include time continuity for non-autonomous systems, or mass continuity
}//============================================================

/**
 *	\brief Retrieve the initial conditions for a segment that the correction
 *	engine will integrate.
 *
 *	Derived models may replace this function to change how the initial conditions
 *	are chosen from the design vector.
 *
 *	\param it a pointer to the corrector's iteration data structure
 *	\param set a pointer to the nodeset being corrected
 *	\param s the ID of the segment being propagated
 *	\param ic a pointer to a 6-element initial state array
 *	\param t0 a pointer to a double representing the initial time (epoch)
 *	\param tof a pointer to a double the time-of-flight on the segment.
 *	@todo No need to input Nodeset as well as the MultShootData pointer
 */
void DynamicsModel::multShoot_getSimICs(const MultShootData *it, const Nodeset *set, int s,
	double *ic, double *t0, double *tof) const{
	(void) set;

	// Retrieve  representative object and get data from free var vec
	MSVarMap_Obj state = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRef_const(s).getOrigin());
	std::copy(it->X.begin()+state.row0, it->X.begin()+state.row0 + state.nRows, ic);

	if(it->bVarTime){
		// Retrieve  representative object and get data from free var vec
		MSVarMap_Obj tof_obj = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
			it->bEqualArcTime ? -1 : s);
		*tof = it->bEqualArcTime ? it->X[tof_obj.row0]/(it->nodeset->getNumSegs()) : it->X[tof_obj.row0];
	}else{
		*tof = it->nodeset->getTOF(s);
	}
	*t0 = 0;
}//============================================================

/**
 *  \brief Compute the value of a slack variable for an inequality constraint.
 *  \details Computing the value of the slack variable can avoid unneccessary 
 *  shooting iterations when the inequality constraint is already met. If the 
 *  inequality constraint is met, the value returned by this function will make
 *  the constraint function evaluate to zero.
 *  
 *  Note: This function should be called after the state variable vector has 
 *  been initialized by the multiple shooting algorithm
 * 
 *  \param it the MultShootData object associated with the multiple shooting process
 *  \param con the inequality constraint for which the slack variable is being computed
 * 
 *  \return The value of the slack variable that minimizes the constraint function
 *  without setting the slack variable to zero
 *  \throws Exception if the constraint type does not include a slack variable
 */
double DynamicsModel::multShoot_getSlackVarVal(const MultShootData *it, Constraint con) const{
	switch(con.getType()){
		case Constraint_tp::MAX_DIST:
		case Constraint_tp::MIN_DIST:
			return multShoot_targetDist_compSlackVar(it, con);
		case Constraint_tp::MAX_DELTA_V:
			return multShoot_targetDeltaV_compSlackVar(it, con);
		default:
			throw Exception("DynamicsModel::multShoot_getSlackVarVal: Cannot compute slack variable values for equality constraints");
	}
}//===========================================================

/**
 *	\brief Compute constraint function and partial derivative values for a constraint
 *	
 *	This function provides a framework for most of the constraints available to all models.
 *	Given an iteration data object, a constraint, and the index of that constraint, we can
 *	compute the value of the constraint function(s) associated with the constraint as well
 *	as the partial derivatives that relate the constraint functions to design variables.
 *	Although this model defines behavior for many of the constraints, it is up to the model
 *	developer to ensure any new constraints are implemented in a derived model, and any
 *	existing constraints are treated fully and appropriately by the new, derived model.
 *
 *	\param it a pointer to the corrector's iteration data structure
 *	\param con the constraint being applied
 *	\param c the index of the constraint within the total constraint vector (which is, in
 *	turn, stored in the iteration data)
 */	
void DynamicsModel::multShoot_applyConstraint(MultShootData *it, Constraint con, int c) const{

	int row0 = it->conRows[c];

	switch(con.getType()){
		case Constraint_tp::CONT_PV:
			// Apply position-velocity continuity constraints
			multShoot_targetCont_State(it, con, row0);
			break;
		case Constraint_tp::SEG_CONT_PV:
			multShoot_targetCont_State_Seg(it, con, row0);
			break;
		case Constraint_tp::SEG_CONT_EX:
			multShoot_targetCont_Ex_Seg(it, con, row0);
			break;
		case Constraint_tp::CONT_EX:
			// Apply extra continuity constraints
			multShoot_targetCont_Ex(it, con, row0);
			break;
		case Constraint_tp::STATE:
			multShoot_targetState(it, con, row0);
			break;
		case Constraint_tp::MATCH_ALL:
			multShoot_targetMatchAll(it, con, row0);
			break;
		case Constraint_tp::MATCH_CUST:
			multShoot_targetMatchCust(it, con, row0);
			break;
		case Constraint_tp::MAX_DIST:
		case Constraint_tp::MIN_DIST:
		case Constraint_tp::DIST:
			multShoot_targetDist(it, con, c);
			break;
		case Constraint_tp::DELTA_V:
		case Constraint_tp::MAX_DELTA_V:
			multShoot_targetDeltaV(it, con, c);
			break;
		case Constraint_tp::TOF:
			multShoot_targetTOF(it, con, row0);
			break;
		case Constraint_tp::APSE:
			multShoot_targetApse(it, con, row0);
			break;	
		default: break;
	}
}//=========================================================

/**
 *	\brief Compute state continuity constraint values and partial derivatives
 *
 *	This function computes and stores the default state continuity constraints. The partial 
 *	derivatives of each node with respect to other nodes and integration time are all 
 *	computed and placed in the appropriate spots in the Jacobian matrix.
 *
 *	Derived models may replace this function.
 *
 *	\param it a pointer to the correctors iteration data structure
 *	\param con the constraint being applied
 *	\param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_State(MultShootData* it, Constraint con, int row0) const{
	int segID = con.getID();	// get segment ID
	std::vector<double> conData = con.getData();

	// Get index of segment
	int segIx = it->nodeset->getSegIx(segID);
	std::vector<double> lastState = it->propSegs[segIx].getStateByIx(-1);
	MatrixXRd stm = it->propSegs[segIx].getSTMByIx(-1);

	// Get index of origin node
	MSVarMap_Obj state0_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRef_const(segID).getOrigin());
	MSVarMap_Obj statef_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRef_const(segID).getTerminus());

	// Loop through conData
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			it->FX[row0+s] = lastState[s] - it->X[statef_var.row0+s];

			// Loop through all design variables for this node (6) and compute partials of F w.r.t. x
			for(unsigned int x = 0; x < coreStates; x++){
				// put STM elements into DF matrix
				it->DF[it->totalFree*(row0+s) + state0_var.row0 + x] = stm(s,x);
				// Negative identity matrix
				if(s == x)
					it->DF[it->totalFree*(row0+s) + statef_var.row0+x] = -1;
			}

			// Compute partials of F w.r.t. times-of-flight
			// Columns of DF based on time constraints
			if(it->bVarTime){
				std::vector<double> lastDeriv = it->propSegs[segIx].getStateDerivByIx(-1);

				// If equal arc time is enabled, place a 1/(n-1) in front of all time derivatives
				double timeCoeff = it->bEqualArcTime ? 1.0/(it->numNodes - 1) : 1.0;

				MSVarMap_Obj tofVar = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
					it->bEqualArcTime ? Linkable::INVALID_ID : segID);
				// If equal arc time is enabled, all time derivatives are in one column
				// int timeCol = it->bEqualArcTime ? 6*it->numNodes : 6*it->numNodes+segIx;
				
				// Column of state time derivatives: [vel; accel; other time derivatives]
				it->DF[it->totalFree*(row0+s) + tofVar.row0] = timeCoeff*lastDeriv[s];
			}
		}
	}
}//====================================================

/**
 *  \brief Compute segment state continuity constraint values and partial derivatives.
 *  
 *  This function computes and stores the position and velocity continuity constraints for
 *  segment-to-segment links. The partial derivatives of this constraint with respect to
 *  the origin node states and times-of-flight on each segment are computed and stored
 *  in the DF matrix. 
 *
 *	Derived models may extend or replace this function.
 *
 *	\param it a pointer to the correctors iteration data structure
 *	\param con the constraint being applied
 *	\param row0 the first row this constraint applies to
 *	
 *	\throws Exception if the constraint data is not properly formatted
 */
void DynamicsModel::multShoot_targetCont_State_Seg(MultShootData *it, Constraint con, int row0) const{
	int segID1 = con.getID();
	int ix = 0;
	int segID2 = static_cast<int>(con.getFirstDataValue(&ix));
	if(ix < 0)
		throw Exception("DynamicsModel::multShoot_targetCont_State_Seg: No segment ID was located in the cosntraint data vector");

	std::vector<double> conData = con.getData();

	// Get segment index
	int segIx1 = it->nodeset->getSegIx(segID1);
	int segIx2 = it->nodeset->getSegIx(segID2);

	MSVarMap_Obj state01_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRef_const(segID1).getOrigin());
	MSVarMap_Obj state02_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRef_const(segID2).getOrigin());

	std::vector<double> state1 = it->propSegs[segIx1].getStateByIx(-1);
	std::vector<double> state2 = it->propSegs[segIx2].getStateByIx(-1);
	MatrixXRd stm1 = it->propSegs[segIx1].getSTMByIx(-1);
	MatrixXRd stm2 = it->propSegs[segIx2].getSTMByIx(-1);

	std::vector<double> lastDeriv1, lastDeriv2;
	MSVarMap_Obj tof1_var(MSVar_tp::TOF), tof2_var(MSVar_tp::TOF);
	double timeCoeff = 1;
	if(it->bVarTime){
		lastDeriv1 = it->propSegs[segIx1].getStateDerivByIx(-1);
		lastDeriv2 = it->propSegs[segIx2].getStateDerivByIx(-1);

		tof1_var = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
			it->bEqualArcTime ? Linkable::INVALID_ID : segID1);
		tof2_var = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
			it->bEqualArcTime ? Linkable::INVALID_ID : segID2);

		// If equal arc time is enabled, place a 1/(n-1) in front of all time derivatives
		timeCoeff = it->bEqualArcTime ? 1.0/(it->numNodes - 1) : 1.0;
	}

	// Loop through conData
	int count = 0;
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			it->FX[row0+count] = (state1[s]- state2[s]);

			// Loop through all six states of the two origin nodes and compute partials w.r.t. state variables
			for(unsigned int x = 0; x < coreStates; x++){
				// put STM elements into DF matrix
				it->DF[it->totalFree*(row0+count) + state01_var.row0+x] = stm1(s,x);
				it->DF[it->totalFree*(row0+count) + state02_var.row0+x] = -stm2(s,x);
			}

			// Compute partials of F w.r.t. times-of-flight
			if(it->bVarTime){
				// Column of state derivatives: [vel; accel; other time derivatives]
				it->DF[it->totalFree*(row0+count) + tof1_var.row0] = timeCoeff*lastDeriv1[s];
				it->DF[it->totalFree*(row0+count) + tof2_var.row0] = -timeCoeff*lastDeriv2[s];
			}

			count++;
		}
	}
}//====================================================

/**
 *	\brief Computes continuity constraints for constraints with the <tt>CONT_EX</tt> type.
 *
 *	In this base model, no behavior is defined for extra constraints. It is intended to enforce
 *	continuity constraints like epoch (time) continuity, mass continuity, etc.
 *
 *	\param it a pointer to the correctors iteration data structure
 *	\param con the constraint being applied
 *	\param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_Ex(MultShootData *it, Constraint con, int row0) const{
	// Do absoluately nothing
	(void)it;
	(void)con;
	(void)row0;
}//======================================================

/**
 *	\brief Computes continuity constraints for constraints with the <tt>Constraint_tp::SEG_CONT_EX</tt> type.
 *
 *	In this base model, no behavior is defined for extra constraints. It is intended to enforce
 *	continuity constraints like epoch (time) continuity, mass continuity, etc.
 *
 *	\param it a pointer to the correctors iteration data structure
 *	\param con the constraint being applied
 *	\param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_Ex_Seg(MultShootData *it, Constraint con, int row0) const{
	// Do absoluately nothing
	(void)it;
	(void)con;
	(void)row0;
}//======================================================
/**
 *	\brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::STATE</tt>.
 *
 *	This method *should* provide full state constraining for any model; the STM and identity 
 *	matrices are used to relate node states and integrated states.
 *
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con the constraint being applied
 *	\param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetState(MultShootData* it, Constraint con, int row0) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());

	int count = 0; 	// Count # rows since some may be skipped (NAN)
	for(unsigned int s = 0; s < con.getData().size(); s++){
		if(!std::isnan(conData[s])){
			it->FX[row0+count] = it->X[state_var.row0+s] - conData[s];
			it->DF[it->totalFree*(row0 + count) + state_var.row0 + s] = 1;
			count++;
		}
	}	
}//=================================================

/**
 *	\brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::MATCH_ALL</tt>
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are applied
 *	to the Jacobian matrix.
 *
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con a copy of the constraint object
 *	\param row0 the index of the row this constraint begins at
 *	
 *	\throws Exception if the constraint data vector is not properly formatted
 */
void DynamicsModel::multShoot_targetMatchAll(MultShootData* it, Constraint con, int row0) const{
	if(con.getData().size() < 1)
		throw Exception("DynamicsModel::multShoot_targetMatchAll: No segment ID was located in the cosntraint data vector");

	// Only allow matching 6 states, not TOF (state 7)
	MSVarMap_Obj state1_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
	MSVarMap_Obj state2_var = it->getVarMap_obj(MSVar_tp::STATE, con.getData()[0]);
	
	for(unsigned int row = 0; row < coreStates; row++){
		// Constrain the states of THIS node to be equal to the node 
		// with index stored in conData[0]
		it->FX[row0+row] = it->X[state1_var.row0+row] - it->X[state2_var.row0+row];

		// Partial of this constraint wrt THIS node = I
		it->DF[it->totalFree*(row0 + row) + state1_var.row0 + row] = 1;

		// Partial of this constraint wrt other node = -I
		it->DF[it->totalFree*(row0 + row) + state2_var.row0 + row] = -1;
	}
}//=============================================

/**
 *	\brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::MATCH_CUST</tt>
 *
 *	This method *should* provide full functionality for any model; Only 1's and 0's are applied
 *	to the Jacobian matrix.
 *
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con a copy of the constraint object
 *	\param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetMatchCust(MultShootData* it, Constraint con, int row0) const{
	int ix = 0;
	int ID2 = static_cast<int>(con.getFirstDataValue(&ix));
	if(ix < 0)
		throw Exception("DynamicsModel::multShoot_targetMatchCust: No segment ID was located in the cosntraint data vector");

	std::vector<double> conData = con.getData();
	MSVarMap_Obj state1_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
	MSVarMap_Obj state2_var = it->getVarMap_obj(MSVar_tp::STATE, ID2);
	int count = 0;
	
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			it->FX[row0 + count] = it->X[state1_var.row0+s] - it->X[state2_var.row0+s];

			// partial of this constraint wrt THIS node = 1
			it->DF[it->totalFree*(row0 + count) + state1_var.row0+s] = 1;

			// partial of this constraint wrt other node = -1
			it->DF[it->totalFree*(row0 + count) + state2_var.row0+s] = -1;

			count++;
		}
	}
}//===============================================

/**
 *	\brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::DIST</tt>, 
 *	<tt>Constraint_tp::MIN_DIST</tt>, or <tt>Constraint_tp::MAX_DIST</tt>
 *
 *	This method *should* provide full functionality for any autonomous model; It calls the getPrimPos() 
 *	functions, which all models define and uses dynamic-independent computations to populate
 *	the constraint vector and Jacobian matrix. Nonautonomous models will need to include time
 *	dependencies.
 *
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con a copy of the constraint object
 *	\param c the index of this constraint in the constraint vector object
 */
void DynamicsModel::multShoot_targetDist(MultShootData* it, Constraint con, int c) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
	int Pix = static_cast<int>(conData[0]);	// index of primary
	int row0 = it->conRows[c];
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, it->sysData);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it->X[state_var.row0+0] - primPos[Pix*3+0];
	double dy = it->X[state_var.row0+1] - primPos[Pix*3+1];
	double dz = it->X[state_var.row0+2] - primPos[Pix*3+2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

	// Compute difference between desired distance and true distance
	it->FX[row0] = h - conData[1];

	// Partials with respect to node position states
	it->DF[it->totalFree*row0 + state_var.row0 + 0] = dx/h;
	it->DF[it->totalFree*row0 + state_var.row0 + 1] = dy/h;
	it->DF[it->totalFree*row0 + state_var.row0 + 2] = dz/h;

	// Extra stuff for inequality constraints
	if(con.getType() == Constraint_tp::MIN_DIST || 
		con.getType() == Constraint_tp::MAX_DIST ){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
			it->slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
		int sign = con.getType() == Constraint_tp::MAX_DIST ? 1 : -1;

		// Subtract squared slack variable from constraint
		it->FX[row0] += sign*it->X[slackCol]*it->X[slackCol];

		// Partial with respect to slack variable
		it->DF[it->totalFree*row0 + slackCol] = sign*2*it->X[slackCol];
	}
}// End of targetDist() =========================================

/**
 *  \brief Compute the value of the slack variable for inequality distance constraints
 *  \details This function computes a value for the slack variable in an
 *  inequality distance constraint. If the constraint is already met by the initial
 *  design, using this value will prevent the multiple shooting algorithm from
 *  searching all over for the propper value.
 * 
 *  \param it the iteration data object for the multiple shooting process
 *  \param con the constraint the slack variable applies to
 *  \return the value of the slack variable
 */
double DynamicsModel::multShoot_targetDist_compSlackVar(const MultShootData* it, Constraint con) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
	int Pix = static_cast<int>(conData[0]);	// index of primary	
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, it->sysData);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it->X[state_var.row0+0] - primPos[Pix*3+0];
	double dy = it->X[state_var.row0+1] - primPos[Pix*3+1];
	double dz = it->X[state_var.row0+2] - primPos[Pix*3+2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance
	int sign = con.getType() == Constraint_tp::MAX_DIST ? 1 : -1;
    double diff = conData[1] - h;

    /*  If diff and sign have the same sign (+/-), then the constraint
     *  is satisfied, so compute the value of the slack variable that 
     *  sets the constraint function equal to zero. Otherwise, choose 
     *  a small value of the slack variable but don't set it to zero as 
     *  that will make the partials zero and will prevent the mulitple
     *  shooting algorithm from updating the slack variable
     */
    return diff*sign > 0 ? sqrt(std::abs(diff)) : 1e-4;
}//==========================================================

/**
 *	\brief Compute partials and constraints for all nodes constrained with <tt>Constraint_tp::DELTA_V</tt> or
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
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con the constraint being applied
 *	\param c the index of the first row for this constraint
 */
void DynamicsModel::multShoot_targetDeltaV(MultShootData* it, Constraint con, int c) const{

	int row0 = it->conRows[c];

	// Don't allow dividing by zero, but otherwise scale by value to keep order ~1
	double dvMax = con.getData()[0] == 0 ? 1 : con.getData()[0];

	// Compute total deltaV magnitude
	double totalDV = 0;
	for(int s = 0; s < it->nodeset->getNumSegs(); s++){
		// compute magnitude of DV between segment s and its terminal point
		// This takes the form v_n,f - v_n+1,0
		double dvx = it->deltaVs[s*3];
		double dvy = it->deltaVs[s*3+1];
		double dvz = it->deltaVs[s*3+2];
		double dvMag = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

		// If dvMag is zero, don't bother computing partials; they're all equal to zero but the
		// computation will try to divide by zero which is problematic...
		if(dvMag > 0){
			totalDV += dvMag;

			// Compute parial w.r.t. terminus node (where velocity is discontinuous)
			double dFdq_n2_data[] = {0, 0, 0, -dvx/dvMag, -dvy/dvMag, -dvz/dvMag};
			Eigen::RowVectorXd dFdq_n2 = Eigen::Map<Eigen::RowVectorXd>(dFdq_n2_data, 1, 6);

			// Get info about the final state/accel of the integrated segment
			MatrixXRd stm = it->propSegs[s].getSTMByIx(-1);

			// Partial w.r.t. integrated path (newSeg) from origin node
			Eigen::RowVectorXd dFdq_nf = -1*dFdq_n2*stm;

			MSVarMap_Obj state0_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRefByIx_const(s).getOrigin());
			MSVarMap_Obj statef_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getSegRefByIx_const(s).getTerminus());

			for(unsigned int i = 0; i < 6; i++){
				it->DF[it->totalFree*row0 + statef_var.row0 + i] += dFdq_n2(0, i)/dvMax;
				it->DF[it->totalFree*row0 + state0_var.row0 + i] += dFdq_nf(0, i)/dvMax;
			}

			// Compute partial w.r.t. segment time-of-flight
			if(it->bVarTime){
				// Derivative of the final state of segment s
				std::vector<double> state_dot_data = it->propSegs[s].getStateDerivByIx(-1);
				Eigen::VectorXd state_dot = Eigen::Map<Eigen::VectorXd>(&(state_dot_data[0]), 6, 1);

				double timeCoeff = it->bEqualArcTime ? 1.0/(it->numNodes - 1) : 1.0;
				MSVarMap_Obj tof_var = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
					it->bEqualArcTime ? Linkable::INVALID_ID : it->nodeset->getSegRefByIx_const(s).getID());

				Eigen::RowVectorXd dFdt_n = -1*dFdq_n2 * state_dot;
				it->DF[it->totalFree*row0 + tof_var.row0] = timeCoeff*dFdt_n(0)/dvMax;
			}
		}
	}
	
	// Copute the difference between the actual deltaV and the desired deltaV
	it->FX[row0] = con.getData()[0] == 0 ? totalDV : totalDV/dvMax - 1;

	if(con.getType() == Constraint_tp::MAX_DELTA_V){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(),
			it->slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
		
		// printf("Constraint_tp::MAX_DELTA_V Constraint Details:\n");
		// printf("Total DV")
		it->FX[row0] += it->X[slackCol]*it->X[slackCol];
		it->DF[it->totalFree*row0 + slackCol] = 2*it->X[slackCol];
	}
}//==============================================

/**
 *  \brief Compute the slack variable value for a delta-V constraint
 *  \details This function currently returns a hard-coded value of 1e-2
 * 
 *  \param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con the constraint being applied
 * 
 *  \return Ideal value of the slack variable
 */
double DynamicsModel::multShoot_targetDeltaV_compSlackVar(const MultShootData *it, Constraint con) const{
	// double totalDV = 0;
	// for(int n = 0; n < it->numNodes-1; n++){
	// 	// compute squared magnitude of DV between node n and n+1
	// 	// This takes the form v_n,f - v_n+1,0
	// 	totalDV += it->deltaVs[n*3]*it->deltaVs[n*3] +
	// 		it->deltaVs[n*3+1]*it->deltaVs[n*3+1] + 
	// 		it->deltaVs[n*3+2]*it->deltaVs[n*3+2];
	// }
	// // Value of the constraint function
	// double F = totalDV - con.getData()[0];
	
	// No info about delta-Vs between segments exists yet because nothing
	// has been integrated... return a constant for now
	(void) it;
	(void) con;
	return 1e-2;

	// If F < 0, the constraint is satisfied, so choose a slack variable
	// that sets F = 0; else choose a small slack variable value
	// return F < 0 ? sqrt(std::abs(F)) : 1e-4;
}//====================================================

/**
 *	\brief Compute partials and constraint function values for time-of-flight constraints
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are
 *	used to relate TOFs.
 *	
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con a copy of the constraint object
 *	\param row0 the index of the row this constraint begins at
 *	
 *	\throws Exception if variable time is set to OFF
 */
void DynamicsModel::multShoot_targetTOF(MultShootData *it, Constraint con, int row0) const{
	if(! it->bVarTime)
		throw Exception("DynamicsModel::multShoot_targetTOF: Cannot target TOF when variable time is off!");

	if(it->bEqualArcTime){
		MSVarMap_Obj tof_var = it->getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
		it->FX[row0] = it->X[tof_var.row0];
		it->DF[it->totalFree*row0 + tof_var.row0] = 1;
	}else{
		// Sum all TOF for total, set partials w.r.t. integration times equal to one
		for(int s = 0; s < it->nodeset->getNumSegs(); s++){
			MSVarMap_Obj tof_var = it->getVarMap_obj(MSVar_tp::TOF, it->nodeset->getSegRefByIx_const(s).getID());
			it->FX[row0] += it->X[tof_var.row0];
			it->DF[it->totalFree*row0 + tof_var.row0] = 1;
		}
	}
	
	// subtract the desired TOF from the constraint to finish its computation
	it->FX[row0] -= con.getData()[0];
}//====================================================

/**
 *	\brief Compute partials and constraint function values for apse constraints
 *
 *	This method *should* provide full functionality for any autonomous model. Non-
 *	autonomous models will need to modify the function to account for epoch time
 *	
 *	\param it a pointer to the class containing all the data relevant to the corrections process
 *	\param con a copy of the constraint object
 *	\param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetApse(MultShootData *it, Constraint con, int row0) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
	int Pix = static_cast<int>(conData[0]);	// index of primary
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time
	
	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, it->sysData);

	// Get distance between node and primary in x, y, and z-coordinates, use non-scaled coordinates
	double dx = it->X[state_var.row0+0] - primPos[Pix*3+0];
	double dy = it->X[state_var.row0+1] - primPos[Pix*3+1];
	double dz = it->X[state_var.row0+2] - primPos[Pix*3+2];
	double vx = it->X[state_var.row0+3];
	double vy = it->X[state_var.row0+4];
	double vz = it->X[state_var.row0+5];

	// Constraint function: r_dot = 0 (using non-scaled coordinates)
	it->FX[row0] = dx*vx + dy*vy + dz*vz;

	// Partials of F w.r.t. node state
	it->DF[it->totalFree*row0 + state_var.row0+0] = vx;
	it->DF[it->totalFree*row0 + state_var.row0+1] = vy;
	it->DF[it->totalFree*row0 + state_var.row0+2] = vz;
	it->DF[it->totalFree*row0 + state_var.row0+3] = dx;
	it->DF[it->totalFree*row0 + state_var.row0+4] = dy;
	it->DF[it->totalFree*row0 + state_var.row0+5] = dz;
}//====================================================




}// END of Astrohelion namespace