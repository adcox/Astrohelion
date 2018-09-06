/**
 *  @file DynamicsModel.cpp
 *	@brief Defines behavior for a dynamic model
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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
#include "Arcset.hpp"
#include "SysData.hpp"
#include "Utilities.hpp"
 
#include <algorithm>
#include <cmath>

#include "AsciiOutput.hpp"

namespace astrohelion{
/**
 *	@brief Default constructor
 *	@param type the model type
 */
DynamicsModel::DynamicsModel(){}

/**
 *	@brief Copy constructor
 */
DynamicsModel::DynamicsModel(const DynamicsModel &m){
	copyMe(m);
}//===========================================

/**
 *	@brief Deconstructor
 */
DynamicsModel::~DynamicsModel(){}

/**
 *	@brief Copy Operator
 */	
DynamicsModel& DynamicsModel::operator =(const DynamicsModel &m){
	copyMe(m);
	return *this;
}//============================================

/**
 *	@brief Copies all data from a dynamic model into this one
 *	@param m another dynamic model
 */
void DynamicsModel::copyMe(const DynamicsModel &m){
	coreDim = m.coreDim;
	extraDim = m.extraDim;
}//============================================

/**
 *	@brief Retrieve the number of core states
 *	@return the number of core states
 */
unsigned int DynamicsModel::getCoreStateSize() const { return coreDim; }

/**
 *	@brief Retrieve the number of extra states stored after the core states and STM elements
 *	@return the number of extra states stored after the core states and STM elements
 */
unsigned int DynamicsModel::getExtraStateSize() const { return extraDim; }

/**
 *	@brief Determine whether the specified constraint type is supported in this model
 *	@return whether or not the specified constraint type is supported in this model
 */
bool DynamicsModel::supportsCon(Constraint_tp type) const{
	return std::find(allowedCons.begin(), allowedCons.end(), type) != allowedCons.end();
}//===================================================

/**
 *	@brief Determine whether the specified event type is supported in this model
 *	@return whether or not the specified event type is supported in this model
 */
bool DynamicsModel::supportsEvent(Event_tp type) const{
	return std::find(allowedEvents.begin(), allowedEvents.end(), type) != allowedEvents.end();
}//===================================================

/**
 *  @brief Determine whether or not the model supports a specific control law
 *  @details By default, a DynamicsModel does not support any control laws.
 * 
 *  @param pLaw Pointer to a control law
 *  @return Whether or not the Dynamics model supports the specified control law
 */
bool DynamicsModel::supportsControl(const ControlLaw *pLaw) const{
	if(pLaw){
		// By default, a non-nullptr control law is not supported. Special code
		return false;
	}else{
		// If pLaw = nullptr, no control is applied; this is always allowed
		return true;
	}
}//===================================================

/**
 *  @brief Construct a new control law and allocated it on the stack.
 *  @details Each dynamic model will return a pointer to the specific control
 *  law applicable to the system / model
 *  
 *  @param id the control law ID
 *  @param params the parameters that create the control law
 *  
 *  @return A pointer to a control law object. The object has been allocated
 *  on the stack so the delete() function must be employed to free the memory
 */
ControlLaw* DynamicsModel::createControlLaw(unsigned int id,
	const std::vector<double> &params) const{

	return new ControlLaw(id, params);
}//====================================================

/**
 *  @brief Determine the time derivative of the magnitude of a vector from a primary to
 *  the body of interest, non-dimensional units
 *  @details This derivation assumes the primary of interest is fixed in the frame the
 *  spacecraft coordinates are expressed in. If this is not the case (i.e., the primary
 *  moves in the working frame), this function will need to be overridden
 * 
 *  @param Pix Index of the primary
 *  @param t Non-dimensional epoch to compute r-dot at
 *  @param state six-element non-dimensional spacecraft state
 *  @param sys system data object
 *  @return Time derivative of the magnitude of a vector from a primary to the body
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
 *  @brief Create default events for a simulation run
 *  @details These events are intended to prevent numerical issues, e.g., to avoid singularities.
 * 
 *  @param pSys pointer to system data object
 *  @return A vector of events to use in the simulation
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
 *  @brief Create a node on the Arcset
 * 
 *  @param node a reference to the node that is being added
 *  @param y pointer to full state data array
 *  @param t current integration time
 *  @param arcset pointer to trajectory in which the data is stored
 *  @param params pointer to a structure that contains parameters used in the 
 *  integration
 *  @param tp Event type that occured at the Node
 *  
 *  @return ID of the node once it is added to the Trajectory
 */
int DynamicsModel::sim_addNode(Node &node, const double *y, double t, 
	Arcset* arcset, EOM_ParamStruct *params, Event_tp tp) const{

	(void) y;
	(void) t;
	
	node.setTriggerEvent(tp);

	// Save the control law information for each node
	if(params->pCtrlLaw){
		unsigned int ctrl_dim = params->pCtrlLaw->getNumStates();
		if(ctrl_dim > 0){
			node.setExtraParamVec(PARAMKEY_CTRL, std::vector<double>(y + coreDim, 
				y + coreDim + ctrl_dim));
		}
	}

	return arcset->addNode(node);
}//====================================================

/**
 *  @brief Add a segment to the Arcset
 *  @details [long description]
 * 
 *  @param seg Reference to the segment to add to the Arcset
 *  @param y pointer to full state data array
 *  @param t current integration time
 *  @param arcset pointer to the Arcset in which the data is stored
 *  @param params pointer to structure that contains parameters used in the integration
 *  @return the ID of the segment once it is added to the arcset
 */
int DynamicsModel::sim_addSeg(Segment &seg, const double *y, double t, Arcset* arcset, EOM_ParamStruct *params) const{
	(void) y;
	(void) t;

	seg.setCtrlLaw(params->pCtrlLaw);
	return arcset->addSeg(seg);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *	@brief Initialize the corrector's design vector with position and velocity states,
 *	and times-of-flight.
 *
 *	Derived models may replace this function or call it and then append more design 
 * 	variables.
 *
 *	@param it a reference to the corrector's iteration data structure
 *	
 *	@throws Exception if equalArcTime is set to true (within the correction engine) and
 *	the nodeset to be corrected includes segments that are propagated in both forward and
 *	reverse time.
 */
void DynamicsModel::multShoot_initDesignVec(MultShootData &it) const{
	// Create the initial state vector
	it.X.clear();
	it.freeVarMap.clear();

	// Copy in the state vector for each node
	int rowNum = 0;
	for(unsigned int n = 0; n < it.pArcIn->getNumNodes(); n++){
		// Determine if this state should be withheld from the free variable vector
		bool stateIsFree = true, ctrlIsFree = true;
		std::vector<Constraint> nodeCons = it.pArcIn->getNodeRefByIx_const(n).getConstraints();
		for(const Constraint &con : nodeCons){
			if(con.getType() == Constraint_tp::RM_STATE)
				stateIsFree = false;

			if(con.getType() == Constraint_tp::RM_CTRL)
				ctrlIsFree = false;
		}

		if(stateIsFree){
			// Insert the state into the free variable vector
			std::vector<double> state = it.pArcIn->getStateByIx(n);
			rowNum = it.X.size();
			it.X.insert(it.X.end(), state.begin(), state.end());
			// Save a key and object that represents the state and remembers its location in 
			// the free variable vector
			MSVarMap_Key key(MSVar_tp::STATE, it.pArcIn->getNodeRefByIx_const(n).getID());
			it.freeVarMap[key] = MSVarMap_Obj(key, rowNum, state.size());
		}else{
			// Create a key and object for the free variable map that will tell multiple
			// shooting functions that the state is not part of the free variable vector
			MSVarMap_Key key(MSVar_tp::STATE, it.pArcIn->getNodeRefByIx_const(n).getID());
			it.freeVarMap[key] = MSVarMap_Obj(key, -1, it.pArcIn->getSysData()->getDynamicsModel()->getCoreStateSize());
		}

		try{
			std::vector<double> ctrlStates = it.pArcIn->getNodeRefByIx_const(n).getExtraParamVec(PARAMKEY_CTRL);

			if(ctrlIsFree){
				rowNum = it.X.size();
				it.X.insert(it.X.end(), ctrlStates.begin(), ctrlStates.end());

				MSVarMap_Key key(MSVar_tp::CTRL, it.pArcIn->getNodeRefByIx_const(n).getID());
				it.freeVarMap[key] = MSVarMap_Obj(key, rowNum, ctrlStates.size());
			}else{
				MSVarMap_Key key(MSVar_tp::CTRL, it.pArcIn->getNodeRefByIx_const(n).getID());
				it.freeVarMap[key] = MSVarMap_Obj(key, -1, ctrlStates.size());
			}
		}catch(const Exception &e){
			// Extra parameter vector doesn't exist; no action necessary
		}
	}

	if(to_underlying(it.tofTp) > 0){	// TOF is not fixed
		switch(it.tofTp){
			case MSTOF_tp::VAR_FREE:
			{
				// Append the TOF for each segment
				for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
					MSVarMap_Key key(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
					it.freeVarMap[key] = MSVarMap_Obj(key, static_cast<int>(it.X.size()));
					it.X.insert(it.X.end(), it.pArcIn->getTOFByIx(s));
				}

				break;
			}
			case MSTOF_tp::VAR_EQUALARC:
			{
				// Make sure all times-of-flight have the same sign
				for(unsigned int s = 1; s < it.pArcIn->getNumSegs(); s++){
					if(it.pArcIn->getTOFByIx(s) * it.pArcIn->getTOFByIx(s-1) < 0)
						throw Exception("DynamicsModel::multShoot_initDesignVec: EqualArcTime is ON and times-of-flight have different signs... cannot proceed");
				}

				// Append the total TOF for the arc
				MSVarMap_Key key(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
				it.freeVarMap[key] = MSVarMap_Obj(key, static_cast<int>(it.X.size()));
				it.X.insert(it.X.end(), it.pArcIn->getTotalTOF());
				break;
			}
			case MSTOF_tp::VAR_FIXSIGN:
			{
				// Append the sqrt(TOF) for each segment
				for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
					MSVarMap_Key key(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
					it.freeVarMap[key] = MSVarMap_Obj(key, static_cast<int>(it.X.size()));
					it.X.insert(it.X.end(), sqrt(std::abs(it.pArcIn->getTOFByIx(s))));
				}

				break;
			}
			default:
				throw Exception("DynamicsModel::multShoot_initDesignVec: Unhandled time type");
		}
	}
}//============================================================

/**
 *	@brief Create continuity constraints for the correction algorithm; this function
 *	creates position and velocity constraints.
 *
 *	Derived models may replace this function or call it and then append more constraints
 *	for other variables that may be continuous, such as time (non-autonomous systems)
 *	or mass. This function assumes the state elements at indices 3, 4, and 5 represent
 *	velocity and their continuity is governed by the boolean flags stored in a segment.
 *
 *	@param it a reference to the corrector's iteration data structure
 */	
void DynamicsModel::multShoot_createContCons(MultShootData &it) const{
	// Create position and velocity constraints
	for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
		// Force all positions to be continuous
		std::vector<double> contStates(coreDim, 1);
		if(it.pArcIn->getSegRefByIx_const(s).getTerminus() != Linkable::INVALID_ID){	
			// Get a vector specifying which velocity states are continuous
			std::vector<bool> velCon = it.pArcIn->getSegRefByIx_const(s).getVelCon();
			// If not continuous, put NAN into the constraint data; else unity
			contStates[3] = velCon[0] ? 1 : NAN;
			contStates[4] = velCon[1] ? 1 : NAN;
			contStates[5] = velCon[2] ? 1 : NAN;
			
			// Create a constraint
			Constraint con(Constraint_tp::CONT_PV, it.pArcIn->getSegRefByIx_const(s).getID(), contStates);

			// Save constraint to constraint vector
			it.allCons.push_back(con);
		}
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
 *	@param it a reference to the corrector's iteration data structure
 *	@param s the ID of the segment being propagated
 *	@param ic a pointer to the initial state array
 *	@param ctrl0 a pointer to the initial control state array
 *	@param t0 a pointer to a double representing the initial time (epoch)
 *	@param tof a pointer to a double the time-of-flight on the segment.
 */
void DynamicsModel::multShoot_getSimICs(const MultShootData &it, int s,
	double *ic, double *ctrl0, double *t0, double *tof) const{

	int segOriginID = it.pArcIn->getSegRef_const(s).getOrigin();
	ControlLaw *pLaw = it.pArcIn->getSegRef_const(s).getCtrlLaw();

	// Retrieve  representative object for state and get data from free var vec
	MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, segOriginID);
	if(state_var.row0 == -1){
		std::vector<double> stateVec = it.pArcIn->getState(state_var.key.id);
		std::copy(stateVec.begin(), stateVec.end(), ic);
	}else{
		std::copy(it.X.begin()+state_var.row0, it.X.begin()+state_var.row0 + state_var.nRows, ic);
	}

	if(pLaw){
		// Retrieve representative object for control and get data from free var vec
		try{
			MSVarMap_Obj ctrl_var = it.getVarMap_obj(MSVar_tp::CTRL, segOriginID);

			// If the control variables are not part of the free variable vector, retrieve from input nodeset:
			if(ctrl_var.row0 == -1){
				std::vector<double> ctrl = it.pArcIn->getNodeRef_const(segOriginID).getExtraParamVec(PARAMKEY_CTRL);
				std::copy(ctrl.begin(), ctrl.end(), ctrl0);
			}else{
				std::copy(it.X.begin()+ctrl_var.row0, it.X.begin()+ctrl_var.row0 + ctrl_var.nRows, ctrl0);
			}
		}catch(const Exception &e){ /* No need to handle exception */ }
	}

	switch(it.tofTp){
		case MSTOF_tp::FIXED:
			*tof = it.pArcIn->getTOF(s);
			break;
		case MSTOF_tp::VAR_FREE:
		{
			MSVarMap_Obj tofObj = it.getVarMap_obj(MSVar_tp::TOF, s);
			*tof = it.X[tofObj.row0];	
			break;
		}
		case MSTOF_tp::VAR_EQUALARC:
		{
			MSVarMap_Obj tofObj = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
			*tof = it.X[tofObj.row0]/(it.pArcIn->getNumSegs());
			break;
		}
		case MSTOF_tp::VAR_FIXSIGN:
		{
			MSVarMap_Obj tofObj = it.getVarMap_obj(MSVar_tp::TOF, s);
			*tof = astrohelion::sign(it.pArcIn->getTOF(s))*(it.X[tofObj.row0])*(it.X[tofObj.row0]);
			break;
		}
	}
	
	// Most currently implemented systems are autonomous, thus, epoch is not a variable
	*t0 = 0;//it.pArcIn->getEpoch(state_var.key.id);
}//============================================================

/**
 *  @brief Compute the value of a slack variable for an inequality constraint.
 *  @details Computing the value of the slack variable can avoid unneccessary 
 *  shooting iterations when the inequality constraint is already met. If the 
 *  inequality constraint is met, the value returned by this function will make
 *  the constraint function evaluate to zero.
 *  
 *  Note: This function should be called after the state variable vector has 
 *  been initialized by the multiple shooting algorithm
 * 
 *  @param it the MultShootData object associated with the multiple shooting process
 *  @param con the inequality constraint for which the slack variable is being computed
 * 
 *  @return The value of the slack variable that minimizes the constraint function
 *  without setting the slack variable to zero
 *  @throws Exception if the constraint type does not include a slack variable
 */
double DynamicsModel::multShoot_getSlackVarVal(const MultShootData &it, const Constraint& con) const{
	switch(con.getType()){
		case Constraint_tp::MAX_DIST:
		case Constraint_tp::MIN_DIST:
			return multShoot_targetDist_compSlackVar(it, con);
		case Constraint_tp::ENDSEG_MAX_DIST:
		case Constraint_tp::ENDSEG_MIN_DIST:
			return multShoot_targetDist_endSeg_compSlackVar(it, con);
		case Constraint_tp::MAX_DELTA_V:
			return multShoot_targetDeltaV_compSlackVar(it, con);
		default:
			throw Exception("DynamicsModel::multShoot_getSlackVarVal: Cannot compute slack variable values for equality constraints");
	}
}//===========================================================

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
 *	@param it a reference to the corrector's iteration data structure
 *	@param con the constraint being applied
 *	@param c the index of the constraint within the total constraint vector (which is, in
 *	turn, stored in the iteration data)
 */	
void DynamicsModel::multShoot_applyConstraint(MultShootData& it, const Constraint& con, int c) const{

	int row0 = it.conRows[c];

	switch(con.getType()){
		case Constraint_tp::CONT_PV:
			multShoot_targetCont_State(it, con, row0);		// Apply position-velocity continuity constraints
			break;
		case Constraint_tp::CONT_CTRL:
			multShoot_targetCont_Ctrl(it, con, row0);		// Control state continuity
			break;
		case Constraint_tp::SEG_CONT_PV:
			multShoot_targetCont_State_Seg(it, con, row0);
			break;
		case Constraint_tp::SEG_CONT_EX:
			multShoot_targetCont_Ex_Seg(it, con, row0);
			break;
		case Constraint_tp::CONT_EX:
			multShoot_targetCont_Ex(it, con, row0);			// Apply extra continuity constraints
			break;
		case Constraint_tp::STATE:
			multShoot_targetState(it, con, row0);
			break;
		case Constraint_tp::CTRL:
			multShoot_targetCtrl(it, con, row0);
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
		case Constraint_tp::TOF_TOTAL:
			multShoot_targetTOF(it, con, row0);
			break;
		case Constraint_tp::APSE:
			multShoot_targetApse(it, con, row0);
			break;	
		case Constraint_tp::ENDSEG_STATE:
			multShoot_targetState_endSeg(it, con, row0);
			break;
		case Constraint_tp::ENDSEG_APSE:
			multShoot_targetApse_endSeg(it, con, row0);
			break;
		case Constraint_tp::ENDSEG_MAX_DIST:
		case Constraint_tp::ENDSEG_MIN_DIST:
		case Constraint_tp::ENDSEG_DIST:
			multShoot_targetDist_endSeg(it, con, c);
			break;
		default: break;
	}
}//=========================================================

/**
 *	@brief Compute core state continuity constraint values and partial derivatives
 *
 *	This function computes and stores the default core state continuity constraints. The partial 
 *	derivatives of each node with respect to other nodes and integration time are all 
 *	computed and placed in the appropriate spots in the Jacobian matrix.
 *
 *	Derived models may replace this function.
 *
 *	@param it a reference to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_State(MultShootData& it, const Constraint& con, int row0) const{
	int segID = con.getID();	// get segment ID
	std::vector<double> conData = con.getData();

	// Get index of segment
	int segIx = it.pArcIn->getSegIx(segID);
	std::vector<double> lastState = it.propSegs[segIx].getStateByIx(-1);
	MatrixXRd stm = it.propSegs[segIx].getSTMByIx(-1);

	// Get index of origin node
	MSVarMap_Obj state0_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(segID).getOrigin());
	MSVarMap_Obj statef_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(segID).getTerminus());

	// Loop through conData
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			double statef_value = statef_var.row0 == -1 ? it.pArcIn->getNodeRef_const(statef_var.key.id).getStateRef_const()[s] : it.X[statef_var.row0 + s];
			it.FX[row0+s] = lastState[s] - statef_value;

			// Loop through all design variables for this node and compute partials of F w.r.t. state[x]
			for(unsigned int x = 0; x < coreDim; x++){
				// put STM elements into DF matrix
				if(state0_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0+s, state0_var.row0 + x, stm(s,x)));

				// Negative identity matrix
				if(statef_var.row0 != -1 && s == x){
					it.DF_elements.push_back(Tripletd(row0+s, statef_var.row0+x, -1.0));
				}
			}

			// Compute partials of F w.r.t. times-of-flight
			// Columns of DF based on time constraints
			if(to_underlying(it.tofTp) > 0){	// Time is variable
				std::vector<double> lastDeriv = it.propSegs[segIx].getStateDerivByIx(-1);
				double timeCoeff = 1;
				MSVarMap_Obj tofVar;

				switch(it.tofTp){
					case MSTOF_tp::VAR_FIXSIGN:
						tofVar = it.getVarMap_obj(MSVar_tp::TOF, segID);
						timeCoeff = astrohelion::sign(it.pArcIn->getTOF(segID))*2*it.X[tofVar.row0];
						break;
					case MSTOF_tp::VAR_FREE:
						tofVar = it.getVarMap_obj(MSVar_tp::TOF, segID);
						break;
					case MSTOF_tp::VAR_EQUALARC:
						tofVar = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
						timeCoeff = 1.0/(it.pArcIn->getNumSegs());
						break;
					default:
						throw Exception("DynamicsModel::multShoot_targetCont_State: Unhandled time type");
				}
				
				// Column of state time derivatives: [vel; accel; other time derivatives]
				it.DF_elements.push_back(Tripletd(row0+s, tofVar.row0, timeCoeff*lastDeriv[s]));
			}

		}
	}

	// If the control law is nontrivial, find the control variable and compute partials 
	// of the continuity constraint w.r.t. control states
	if(ControlLaw *pLaw = it.pArcIn->getSegRef_const(segID).getCtrlLaw()){
		// Check to make sure the number of states is nontrivial; if numStates = 0, there won't be any MSVarMap_obj stored
		if(pLaw->getNumStates() > 0){
			MSVarMap_Obj ctrl0_var = it.getVarMap_obj(MSVar_tp::CTRL, it.pArcIn->getSegRef_const(segID).getOrigin());

			if(ctrl0_var.row0 != -1){
				for(unsigned int s = 0; s < conData.size(); s++){
					if(!std::isnan(conData[s])){
						for(unsigned int c = 0; c < pLaw->getNumStates(); c++){
							// put STM elements into DF matrix
							it.DF_elements.push_back(Tripletd(row0+s, ctrl0_var.row0+c, stm(s, coreDim + c)));
						}
					}
				}
			}
		}
	}
}//====================================================

/**
 *	@brief Compute control state continuity constraint values and partial 
 *	derivatives
 *
 *	@details This function computes and stores control state continuity 
 *	constraints. The partial derivatives of each node with respect to other 
 *	nodes and integration times are all computed and placed in the appropriate 
 *	spots in the Jacobian matrix.
 *
 *	Derived models may replace this function.
 *
 *	@param it a reference to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_Ctrl(MultShootData& it, 
	const Constraint& con, int row0) const{
	
	int segID = con.getID();	// get segment ID
	std::vector<double> conData = con.getData();

	int segIx = it.pArcIn->getSegIx(segID);
	const Segment &propSeg = it.propSegs[segIx].getSegRefByIx_const(0);
	const Segment &inputSeg = it.pArcIn->getSegRef_const(segID);

	if(propSeg.getCtrlLaw() == nullptr)
		throw Exception("DynamicsModel::multShoot_targetCont_Ctrl: "
			"Constrained segment has a nullptr control law");
	

	const unsigned int ctrlDim = propSeg.getCtrlLaw()->getNumStates();

	std::vector<double> lastFullState = propSeg.getStateByRow(-1);
	
	MatrixXRd stm = propSeg.getSTM();

	MSVarMap_Obj ctrl0_var = it.getVarMap_obj(MSVar_tp::CTRL, 
		inputSeg.getOrigin());
	MSVarMap_Obj nodeState0_var = it.getVarMap_obj(MSVar_tp::STATE, 
		inputSeg.getOrigin());

	MSVarMap_Obj ctrlf_var = it.getVarMap_obj(MSVar_tp::CTRL, 
		inputSeg.getTerminus());

	std::vector<double> ctrlf_static_vector = 
		it.pArcIn->getNodeRef_const(ctrlf_var.key.id).\
		getExtraParamVec(PARAMKEY_CTRL);

	// Loop through conData
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			double ctrlf_value = ctrlf_var.row0 == -1 ?\
				ctrlf_static_vector[s] : it.X[ctrlf_var.row0 + s];

			it.FX[row0+s] = lastFullState[coreDim + s] - ctrlf_value;

			// Compute partials of F w.r.t. origin node state
			for(unsigned int x = 0; x < coreDim; x++){
				// put STM elements into DF matrix
				if(nodeState0_var.row0 != -1){
					it.DF_elements.push_back(Tripletd(row0+s, 
						nodeState0_var.row0 + x, stm(coreDim + s, x)));
				}
			}

			// Compute partials of F w.r.t. origin and terminal node ctrl states
			for(unsigned int x = 0; x < ctrlDim; x++){
				// For origin node: Put STM elements into DF matrix
				if(ctrl0_var.row0 != -1){
					it.DF_elements.push_back(Tripletd(row0+s, 
						ctrl0_var.row0 + x, stm(coreDim + s, coreDim + x)));
				}

				// For terminal node: Negative identity
				if(ctrlf_var.row0 != -1 && s == x){
					it.DF_elements.push_back(Tripletd(row0+s, 
						ctrlf_var.row0 + x, -1.0));
				}
			}

			// Compute partials of F w.r.t. times-of-flight
			// Columns of DF based on time constraints
			if(to_underlying(it.tofTp) > 0){
				/*
				 * For now, all time-derivatives of control laws are zero
				 *
				 * TODO - Will need to retrieve constrol state time derivatives 
				 * from propagated segment
				 * * Instantaneous time derivatives are available via 
				 *   ControlLaw.getTimeDeriv()
				 */

				// std::vector<double> lastDeriv = propSeg.getStateDerivByIx(-1);

				// // If equal arc time is enabled, place a 1/(n-1) in front of  
				// // all time derivatives
				// double timeCoeff = it.bEqualArcTime ?\
				// 	1.0/(it.pArcIn->getNumSegs()) : 1.0;

				// MSVarMap_Obj tofVar = it.getVarMap_obj(it.bEqualArcTime ?\
				// MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
				// 	it.bEqualArcTime ? Linkable::INVALID_ID : segID);
				
				// // Column of state time derivatives: [vel; accel; other 
				// // time derivatives]
				// it.DF_elements.push_back(Tripletd(row0+s, tofVar.row0, 
				// timeCoeff*lastDeriv[s]));
			}

		}
	}
}//====================================================

/**
 *  @brief Compute segment state continuity constraint values and partial derivatives.
 *  
 *  This function computes and stores the position and velocity continuity constraints for
 *  segment-to-segment links. The partial derivatives of this constraint with respect to
 *  the origin node states and times-of-flight on each segment are computed and stored
 *  in the DF matrix. 
 *
 *	Derived models may extend or replace this function.
 *
 *	@param it a reference to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 *	
 *	@throws Exception if the constraint data is not properly formatted
 */
void DynamicsModel::multShoot_targetCont_State_Seg(MultShootData& it, const Constraint& con, int row0) const{
	int segID1 = con.getID();
	int ix = 0;
	int segID2 = static_cast<int>(con.getFirstDataValue(&ix));
	if(ix < 0)
		throw Exception("DynamicsModel::multShoot_targetCont_State_Seg: No segment ID was located in the cosntraint data vector");

	std::vector<double> conData = con.getData();

	// Get segment index
	int segIx1 = it.pArcIn->getSegIx(segID1);
	int segIx2 = it.pArcIn->getSegIx(segID2);

	MSVarMap_Obj state01_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(segID1).getOrigin());
	MSVarMap_Obj state02_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(segID2).getOrigin());

	std::vector<double> state1 = it.propSegs[segIx1].getStateByIx(-1);
	std::vector<double> state2 = it.propSegs[segIx2].getStateByIx(-1);
	MatrixXRd stm1 = it.propSegs[segIx1].getSTMByIx(-1);
	MatrixXRd stm2 = it.propSegs[segIx2].getSTMByIx(-1);

	std::vector<double> lastDeriv1, lastDeriv2;
	MSVarMap_Obj tof1_var(MSVar_tp::TOF), tof2_var(MSVar_tp::TOF);
	double timeCoeff1 = 1, timeCoeff2 = 1;
	if(to_underlying(it.tofTp) > 0){	// Time is variable
		lastDeriv1 = it.propSegs[segIx1].getStateDerivByIx(-1);
		lastDeriv2 = it.propSegs[segIx2].getStateDerivByIx(-1);

		tof1_var = it.getVarMap_obj(it.tofTp == MSTOF_tp::VAR_EQUALARC ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
			it.tofTp == MSTOF_tp::VAR_EQUALARC ? Linkable::INVALID_ID : segID1);
		tof2_var = it.getVarMap_obj(it.tofTp == MSTOF_tp::VAR_EQUALARC ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
			it.tofTp == MSTOF_tp::VAR_EQUALARC ? Linkable::INVALID_ID : segID2);

		switch(it.tofTp){
			case MSTOF_tp::VAR_FREE: break; 	// Leave both coefficients as unity
			case MSTOF_tp::VAR_FIXSIGN:
				timeCoeff1 = astrohelion::sign(it.pArcIn->getTOF(segID1))*2*it.X[tof1_var.row0];
				timeCoeff2 = astrohelion::sign(it.pArcIn->getTOF(segID2))*2*it.X[tof2_var.row0];
				break;
			case MSTOF_tp::VAR_EQUALARC:
				timeCoeff1 = 1/(it.pArcIn->getNumSegs());
				timeCoeff2 = timeCoeff1;
				break;
			default:
				throw Exception("DynamicsModel::multShoot_targetCont_State_Seg: Unhandled time type");
		}
	}

	// Loop through conData
	int count = 0;
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			// This state is constrained to be continuous; compute error
			it.FX[row0+count] = (state1[s]- state2[s]);

			// Loop through all six states of the two origin nodes and compute partials w.r.t. state variables
			for(unsigned int x = 0; x < coreDim; x++){
				// put STM elements into DF matrix
				if(state01_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0+count, state01_var.row0+x, stm1(s,x)));

				if(state02_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0+count, state02_var.row0+x, -stm2(s,x)));
			}

			// Compute partials of F w.r.t. times-of-flight
			if(to_underlying(it.tofTp) > 0){	// Time is variable
				// Column of state derivatives: [vel; accel; other time derivatives]
				it.DF_elements.push_back(Tripletd(row0+count, tof1_var.row0, timeCoeff1*lastDeriv1[s]));
				it.DF_elements.push_back(Tripletd(row0+count, tof2_var.row0, -timeCoeff2*lastDeriv2[s]));
			}

			count++;
		}
	}

	// If the control law is nontrivial, find the control variable and compute partials 
	// of the continuity constraint w.r.t. control states
	if(ControlLaw *pLaw = it.pArcIn->getSegRef_const(segID1).getCtrlLaw()){
		// Check to make sure the number of states is nontrivial; if numStates = 0, there won't be any MSVarMap_obj stored
		if(pLaw->getNumStates() > 0){
			MSVarMap_Obj ctrl1_var = it.getVarMap_obj(MSVar_tp::CTRL, it.pArcIn->getSegRef_const(segID1).getOrigin());

			if(ctrl1_var.row0 != -1){
				count = 0;
				for(unsigned int s = 0; s < conData.size(); s++){
					if(!std::isnan(conData[s])){
						for(unsigned int c = 0; c < pLaw->getNumStates(); c++){
							// put STM elements into DF matrix
							it.DF_elements.push_back(Tripletd(row0+count, ctrl1_var.row0+c, stm1(s, coreDim + c)));
						}
					}
				}
			}
		}
	}

	// Same for the other segment
	if(ControlLaw *pLaw = it.pArcIn->getSegRef_const(segID2).getCtrlLaw()){
		// Check to make sure the number of states is nontrivial; if numStates = 0, there won't be any MSVarMap_obj stored
		if(pLaw->getNumStates() > 0){
			MSVarMap_Obj ctrl2_var = it.getVarMap_obj(MSVar_tp::CTRL, it.pArcIn->getSegRef_const(segID2).getOrigin());

			if(ctrl2_var.row0 != -1){
				count = 0;
				for(unsigned int s = 0; s < conData.size(); s++){
					if(!std::isnan(conData[s])){
						for(unsigned int c = 0; c < pLaw->getNumStates(); c++){
							// put STM elements into DF matrix
							it.DF_elements.push_back(Tripletd(row0+count, ctrl2_var.row0+c, stm1(s, coreDim + c)));
						}
					}
				}
			}
		}
	}
}//====================================================

/**
 *	@brief Computes continuity constraints for constraints with the `CONT_EX` type.
 *
 *	In this base model, no behavior is defined for extra constraints. It is intended to enforce
 *	continuity constraints like epoch (time) continuity, mass continuity, etc.
 *
 *	@param it a reference to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_Ex(MultShootData& it, const Constraint& con, int row0) const{
	// Do absoluately nothing
	(void)it;
	(void)con;
	(void)row0;
}//======================================================

/**
 *	@brief Computes continuity constraints for constraints with the `Constraint_tp::SEG_CONT_EX` type.
 *
 *	In this base model, no behavior is defined for extra constraints. It is intended to enforce
 *	continuity constraints like epoch (time) continuity, mass continuity, etc.
 *
 *	@param it a reference to the correctors iteration data structure
 *	@param con the constraint being applied
 *	@param row0 the first row this constraint applies to
 */
void DynamicsModel::multShoot_targetCont_Ex_Seg(MultShootData& it, const Constraint& con, int row0) const{
	// Do absoluately nothing
	(void)it;
	(void)con;
	(void)row0;
}//======================================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with `Constraint_tp::STATE`.
 *
 *	This method *should* provide full state constraining for any model.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 *	@param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetState(MultShootData& it, const Constraint& con, int row0) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());

	if(conData.size() > coreDim)
		throw Exception("DynamicsModel::multShoot_targetState: ConData has too many states");
	
	if(state_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetState: Cannot constrain state that is not in the free variable vector.");

	int count = 0; 	// Count # rows since some may be skipped (NAN)
	for(unsigned int s = 0; s < con.getData().size(); s++){
		if(!std::isnan(conData[s])){
			it.FX[row0+count] = it.X[state_var.row0+s] - conData[s];
			it.DF_elements.push_back(Tripletd(row0+count, state_var.row0+s, 1.0));
			count++;
		}
	}
}//=================================================

/**
 *  @brief Compute constraint value and partial derivatives of the constraint
 * 
 *  @param it a reference to iteration data object
 *  @param con Constraint object
 *  @param row0 The row of the constraint within the constraint vector
 */
void DynamicsModel::multShoot_targetState_endSeg(MultShootData& it, const Constraint& con, int row0) const{
	std::vector<double> conData = con.getData();
	int segID = con.getID();
	int segIx = it.pArcIn->getSegIx(segID);
	
	// Get object representing origin of segment
	MSVarMap_Obj prevNode_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(con.getID()).getOrigin());

	MSVarMap_Obj tof_var;
	if(to_underlying(it.tofTp) > 0){
		switch(it.tofTp){
			case MSTOF_tp::VAR_FREE:
			case MSTOF_tp::VAR_FIXSIGN:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF, con.getID());
				break;
			case MSTOF_tp::VAR_EQUALARC:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
				break;
			default:
				throw Exception("DynamicsModel::multShoot_targetState_endSeg: Unhandled time type");
		}
	}
	
	// Data associated with the previous node and the propagated segment
	std::vector<double> lastState = it.propSegs[segIx].getStateByIx(-1);
	std::vector<double> lastDeriv = it.propSegs[segIx].getStateDerivByIx(-1);
	MatrixXRd stm = it.propSegs[segIx].getSTMByIx(-1);

	if(conData.size() > coreDim)
		throw Exception("DynamicsModel::multShoot_targetState_endSeg: ConData has too many states");

	int count = 0; 	// Count # rows since some may be skipped (NAN)
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			it.FX[row0+count] = lastState[s] - conData[s];

			// Partials of F w.r.t. previous node states
			for(unsigned int x = 0; x < coreDim; x++){
				// put STM elements into DF matrix
				if(prevNode_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0+count, prevNode_var.row0 + x, stm(s,x)));
			}

			// Partials of F w.r.t. time-of-flight
			if(to_underlying(it.tofTp) > 0){
				double timeCoeff = 1;
				switch(it.tofTp){
					case MSTOF_tp::VAR_FREE: break;	// Leave timeCoeff = 1
					case MSTOF_tp::VAR_FIXSIGN:
						timeCoeff = astrohelion::sign(it.pArcIn->getTOF(con.getID()))*2*it.X[tof_var.row0];
						break;
					case MSTOF_tp::VAR_EQUALARC:
						timeCoeff = 1.0/(it.pArcIn->getNumSegs());
						break;
					default:
						throw Exception("DynamicsModel::multShoot_targetState_endSeg: Unhandled time type");
				}

				// Column of state time derivatives: {vel; accel};
				it.DF_elements.push_back(Tripletd(row0+count, tof_var.row0, timeCoeff*lastDeriv[s]));
			}
			count++;
		}
	}

	// Partials w.r.t. control states
	if(ControlLaw *pLaw = it.pArcIn->getSegRef_const(segID).getCtrlLaw()){
		if(pLaw->getNumStates() > 0){
			MSVarMap_Obj ctrl_var = it.getVarMap_obj(MSVar_tp::CTRL, it.pArcIn->getSegRef_const(segID).getOrigin());

			if(ctrl_var.row0 != -1){
				for(unsigned int s = 0; s < conData.size(); s++){
					if(!std::isnan(conData[s])){
						for(unsigned int c = 0; c < pLaw->getNumStates(); c++){
							// put STM elements into DF matrix
							it.DF_elements.push_back(Tripletd(row0+s, ctrl_var.row0+c, stm(s, coreDim + c)));
						}
					}
				}
			}
		}
	}
}//=================================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with `Constraint_tp::CTRL`.
 *
 *	This method *should* provide full ctrl state constraining for any model.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 *	@param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetCtrl(MultShootData& it, const Constraint& con, int row0) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj ctrl_var = it.getVarMap_obj(MSVar_tp::CTRL, con.getID());

	if(ctrl_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetCtrl: Cannot consrain ctrl that is not in the free variable vector.");

	int count = 0;	// Count # rows b/c some may be skipped (NAN)
	for(unsigned int s = 0; s < con.getData().size(); s++){
		if(!std::isnan(conData[s])){
			it.FX[row0+count] = it.X[ctrl_var.row0+s] - conData[s];
			it.DF_elements.push_back(Tripletd(row0+count, ctrl_var.row0+s, 1.0));
			count++;
		}
	}
}//====================================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with `Constraint_tp::MATCH_ALL`
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are applied
 *	to the Jacobian matrix.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 *	
 *	@throws Exception if the constraint data vector is not properly formatted
 */
void DynamicsModel::multShoot_targetMatchAll(MultShootData& it, const Constraint& con, int row0) const{
	if(con.getData().size() < 1)
		throw Exception("DynamicsModel::multShoot_targetMatchAll: No segment ID was located in the cosntraint data vector");

	// Only allow matching 6 states, not TOF (state 7)
	MSVarMap_Obj state1_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());
	MSVarMap_Obj state2_var = it.getVarMap_obj(MSVar_tp::STATE, con.getData()[0]);
	
	if(state1_var.row0 == -1 && state2_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetMatchAll: Neither state vector is free; constraint is uncontrollabel.");

	const double *state1 = state1_var.row0 == -1 ? &(it.pArcIn->getNodeRef_const(state1_var.key.id).getStateRef_const().front()) : &(it.X[state1_var.row0]);
	const double *state2 = state2_var.row0 == -1 ? &(it.pArcIn->getNodeRef_const(state2_var.key.id).getStateRef_const().front()) : &(it.X[state2_var.row0]);

	for(unsigned int row = 0; row < coreDim; row++){
		// Constrain the states of THIS node to be equal to the node 
		// with index stored in conData[0]
		it.FX[row0+row] = state1[row] - state2[row];

		// Partial of this constraint wrt THIS node = I
		if(state1_var.row0 != -1)
			it.DF_elements.push_back(Tripletd(row0+row, state1_var.row0+row, 1.0));

		// Partial of this constraint wrt other node = -I
		if(state2_var.row0 != -1)
			it.DF_elements.push_back(Tripletd(row0+row, state2_var.row0+row, -1.0));
	}
}//=============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with 
 *	`Constraint_tp::MATCH_CUST`
 *
 *	This method *should* provide full functionality for any model; Only 1's and 
 *	0's are applied to the Jacobian matrix.
 *
 *	@param it a reference to the class containing all the data relevant to the 
 *	corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetMatchCust(MultShootData& it, 
	const Constraint& con, int row0) const{

	int ix = 0;
	int ID2 = static_cast<int>(con.getFirstDataValue(&ix));
	if(ix < 0)
		throw Exception("DynamicsModel::multShoot_targetMatchCust: "
			"No segment ID was located in the cosntraint data vector");

	std::vector<double> conData = con.getData();
	MSVarMap_Obj state1_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());
	MSVarMap_Obj state2_var = it.getVarMap_obj(MSVar_tp::STATE, ID2);
	
	if(state1_var.row0 == -1 && state2_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetMatchAll: Neither "
			"state vector is free; constraint is uncontrollabel.");

	const double *state1 = state1_var.row0 == -1 ? 
		&(it.pArcIn->getNodeRef_const(state1_var.key.id).getStateRef_const().front()) : 
		&(it.X[state1_var.row0]);
	const double *state2 = state2_var.row0 == -1 ? 
		&(it.pArcIn->getNodeRef_const(state2_var.key.id).getStateRef_const().front()) : 
		&(it.X[state2_var.row0]);

	int count = 0;
	for(unsigned int s = 0; s < conData.size(); s++){
		if(!std::isnan(conData[s])){
			it.FX[row0 + count] = state1[s] - state2[s];

			// partial of this constraint wrt THIS node = 1
			if(state1_var.row0 != -1){
				it.DF_elements.push_back(Tripletd(row0+count, 
					state1_var.row0+s, 1.0));
			}

			// partial of this constraint wrt other node = -1
			if(state2_var.row0 != -1){
				it.DF_elements.push_back(Tripletd(row0+count, 
					state2_var.row0+s, -1.0));
			}

			count++;
		}
	}
}//===============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with `Constraint_tp::DIST`, 
 *	`Constraint_tp::MIN_DIST`, or `Constraint_tp::MAX_DIST`
 *
 *	This method *should* provide full functionality for any autonomous model; It calls the getPrimPos() 
 *	functions, which all models define and uses dynamic-independent computations to populate
 *	the constraint vector and Jacobian matrix. Nonautonomous models will need to include time
 *	dependencies.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param c the index of this constraint in the constraint vector object
 */
void DynamicsModel::multShoot_targetDist(MultShootData& it, const Constraint& con, int c) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());

	if(state_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetDist: State vector is not part of free variable vector; cannot target distance.");

	int Pix = static_cast<int>(conData[0]);	// index of primary
	int row0 = it.conRows[c];
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	// Get the primary position
	double primPos[3] = {0};
	getPrimPos(t, it.pArcIn->getSysData(), Pix, primPos);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it.X[state_var.row0+0] - primPos[0];
	double dy = it.X[state_var.row0+1] - primPos[1];
	double dz = it.X[state_var.row0+2] - primPos[2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

	// Compute difference between desired distance and true distance
	it.FX[row0] = h - conData[1];

	// Partials with respect to node position states
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+0, dx/h));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+1, dy/h));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+2, dz/h));

	// Extra stuff for inequality constraints
	if(con.getType() == Constraint_tp::MIN_DIST || 
		con.getType() == Constraint_tp::MAX_DIST ){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it.slackAssignCon.begin(), 
			it.slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it.totalFree - it.numSlack + (slackIx - it.slackAssignCon.begin());
		int sign = con.getType() == Constraint_tp::MAX_DIST ? 1 : -1;

		// Subtract squared slack variable from constraint
		it.FX[row0] += sign*it.X[slackCol]*it.X[slackCol];

		// Partial with respect to slack variable
		it.DF_elements.push_back(Tripletd(row0, slackCol, sign*2*it.X[slackCol]));
	}
}// End of targetDist() =========================================

/**
 *	@brief Compute partials and constraint functions for segments constrained with `Constraint_tp::ENDSEG_DIST`, 
 *	`Constraint_tp::ENDSEG_MIN_DIST`, or `Constraint_tp::ENDSEG_MAX_DIST`
 *
 *	This method *should* provide full functionality for any autonomous model; It calls the getPrimPos() 
 *	functions, which all models define and uses dynamic-independent computations to populate
 *	the constraint vector and Jacobian matrix. Nonautonomous models will need to include time
 *	dependencies.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param c the index of this constraint in the constraint vector object
 */
void DynamicsModel::multShoot_targetDist_endSeg(MultShootData& it, const Constraint& con, int c) const{
	std::vector<double> conData = con.getData();
	int segIx = it.pArcIn->getSegIx(con.getID());

	// Get object representing origin node state of segment
	MSVarMap_Obj prevNode_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(con.getID()).getOrigin());
	
	MSVarMap_Obj tof_var;
	double timeCoeff = 1;
	if(to_underlying(it.tofTp) > 0){
		switch(it.tofTp){
			case MSTOF_tp::VAR_FREE:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF, con.getID());
				break;
			case MSTOF_tp::VAR_FIXSIGN:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF, con.getID());
				timeCoeff = astrohelion::sign(it.pArcIn->getTOF(con.getID()))*2*it.X[tof_var.row0];
				break;
			case MSTOF_tp::VAR_EQUALARC:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
				timeCoeff = 1.0/(it.pArcIn->getNumSegs());
				break;
			default:
				throw Exception("DynamicsModel::multShoot_targetDist_endSeg: Unhandled time type");
		}
	}

	int Pix = static_cast<int>(conData[0]);	// index of primary
	int row0 = it.conRows[c];
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	// Get the primary position
	double primPos[3] = {0};
	getPrimPos(t, it.pArcIn->getSysData(), Pix, primPos);

	std::vector<double> lastState = it.propSegs[segIx].getStateByIx(-1);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = lastState[0] - primPos[0];
	double dy = lastState[1] - primPos[1];
	double dz = lastState[2] - primPos[2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

	// Compute difference between desired distance and true distance
	it.FX[row0] = h - conData[1];

	double dFdr_nf[3] = {dx/h, dy/h, dz/h};

	// Partials with respect to node position states
	if(prevNode_var.row0 != -1){
		MatrixXRd stm = it.propSegs[segIx].getSTMByIx(-1);

		// Do matrix multiplication with loops to avoid expensive vector allocation
		double sum;
		for(unsigned int c = 0; c < 6; c++){
			sum = 0;
			for(unsigned int r = 0; r < 3; r++){
				sum += dFdr_nf[r]*stm(r,c);
			}

			it.DF_elements.push_back(Tripletd(row0, prevNode_var.row0+c, sum));
		}
	}

	// Partials with respect to time-of-flight
	if(to_underlying(it.tofTp) > 0){
		std::vector<double> lastDeriv = it.propSegs[segIx].getStateDerivByIx(-1);

		// Compute dot product between dFdr_nf and final state derivative vector
		double dp = dFdr_nf[0]*lastDeriv[0] + dFdr_nf[1]*lastDeriv[1] + dFdr_nf[2]*lastDeriv[2];

		it.DF_elements.push_back(Tripletd(row0, tof_var.row0, timeCoeff*dp));
	}

	// Extra stuff for inequality constraints
	if(con.getType() == Constraint_tp::ENDSEG_MIN_DIST || 
		con.getType() == Constraint_tp::ENDSEG_MAX_DIST ){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it.slackAssignCon.begin(), 
			it.slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it.totalFree - it.numSlack + (slackIx - it.slackAssignCon.begin());
		int sign = con.getType() == Constraint_tp::ENDSEG_MAX_DIST ? 1 : -1;

		// Subtract squared slack variable from constraint
		it.FX[row0] += sign*it.X[slackCol]*it.X[slackCol];

		// Partial with respect to slack variable
		it.DF_elements.push_back(Tripletd(row0, slackCol, sign*2*it.X[slackCol]));
	}
}// End of targetDist_endSeg() ====================================

/**
 *  @brief Compute the value of the slack variable for inequality distance constraints
 *  @details This function computes a value for the slack variable in an
 *  inequality distance constraint. If the constraint is already met by the initial
 *  design, using this value will prevent the multiple shooting algorithm from
 *  searching all over for the propper value.
 * 
 *  @param it the iteration data object for the multiple shooting process
 *  @param con the constraint the slack variable applies to
 *  @return the value of the slack variable
 */
double DynamicsModel::multShoot_targetDist_compSlackVar(const MultShootData& it, const Constraint& con) const{
	std::vector<double> conData = con.getData();

	MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());
	int Pix = static_cast<int>(conData[0]);	// index of primary	
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	if(state_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetDist: State vector is not part of free variable vector; cannot target distance.");

	// Get the primary position
	double primPos[3] = {0};
	getPrimPos(t, it.pArcIn->getSysData(), Pix, primPos);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it.X[state_var.row0+0] - primPos[0];
	double dy = it.X[state_var.row0+1] - primPos[1];
	double dz = it.X[state_var.row0+2] - primPos[2];

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
 *  @brief Compute the value of the slack variable for inequality distance constraints
 *  @details This function computes a value for the slack variable in an
 *  inequality distance constraint. If the constraint is already met by the initial
 *  design, using this value will prevent the multiple shooting algorithm from
 *  searching all over for the propper value.
 * 
 *  @param it the iteration data object for the multiple shooting process
 *  @param con the constraint the slack variable applies to
 *  @return the value of the slack variable
 */
double DynamicsModel::multShoot_targetDist_endSeg_compSlackVar(const MultShootData& it, const Constraint& con) const{
	std::vector<double> conData = con.getData();

	std::vector<double> lastState = it.pArcIn->getSegRef_const(con.getID()).getStateByRow(-1);
	
	int Pix = static_cast<int>(conData[0]);	// index of primary
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time

	// Get the primary position
	double primPos[3] = {0};
	getPrimPos(t, it.pArcIn->getSysData(), Pix, primPos);

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = lastState[0] - primPos[0];
	double dy = lastState[1] - primPos[1];
	double dz = lastState[2] - primPos[2];

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
 *	@brief Compute partials and constraints for all nodes constrained with `Constraint_tp::DELTA_V` or
 *	`MIN_DELTA_V`
 *
 *	Because the delta-V constraint applies to the entire trajectory, the constraint function values
 *	and partial derivatives must be computed for each node along the trajectory. This function
 *	takes care of all of them at once.
 *
 *	This function computes the partials and constraints for an autonomous system; Non-autonomous
 *	dynamic models should implement a similar copy of this function that computes partials
 *	with respect to their additional design variables.
 *
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 *	@param c the index of the first row for this constraint
 */
void DynamicsModel::multShoot_targetDeltaV(MultShootData& it, const Constraint& con, int c) const{

	int row0 = it.conRows[c];

	// Don't allow dividing by zero, but otherwise scale by value to keep order ~1
	double dvMax = con.getData()[0] == 0 ? 1 : con.getData()[0];

	// Compute total deltaV magnitude
	double totalDV = 0;
	for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
		// compute magnitude of DV between segment s and its terminal point
		// This takes the form v_n,f - v_n+1,0
		double dvx = it.deltaVs[s*3];
		double dvy = it.deltaVs[s*3+1];
		double dvz = it.deltaVs[s*3+2];
		double dvMag = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

		// If dvMag is zero, don't bother computing partials; they're all equal to zero but the
		// computation will try to divide by zero which is problematic...
		if(dvMag > 0){
			totalDV += dvMag;

			// Compute parial w.r.t. terminus node (where velocity is discontinuous)
			double dFdq_n2_data[] = {0, 0, 0, -dvx/dvMag, -dvy/dvMag, -dvz/dvMag};
			Eigen::RowVectorXd dFdq_n2 = Eigen::Map<Eigen::RowVectorXd>(dFdq_n2_data, 1, 6);

			// Get info about the final state/accel of the integrated segment
			MatrixXRd stm = it.propSegs[s].getSTMByIx(-1);

			// Partial w.r.t. integrated path (newSeg) from origin node
			Eigen::RowVectorXd dFdq_nf = -1*dFdq_n2*stm;

			MSVarMap_Obj state0_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRefByIx_const(s).getOrigin());
			MSVarMap_Obj statef_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRefByIx_const(s).getTerminus());

			for(unsigned int i = 0; i < 6; i++){
				if(statef_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0, statef_var.row0+i, dFdq_n2(0,i)/dvMax));

				if(state0_var.row0 != -1)
					it.DF_elements.push_back(Tripletd(row0, state0_var.row0+i, dFdq_nf(0,i)/dvMax));
			}

			// Compute partial w.r.t. segment time-of-flight
			if(to_underlying(it.tofTp) > 0){	// Time is variable
				// Derivative of the final state of segment s
				std::vector<double> state_dot_data = it.propSegs[s].getStateDerivByIx(-1);
				Eigen::VectorXd state_dot = Eigen::Map<Eigen::VectorXd>(&(state_dot_data[0]), 6, 1);

				double timeCoeff = 1;
				MSVarMap_Obj tof_var;
				switch(it.tofTp){
					case MSTOF_tp::VAR_FIXSIGN:
						tof_var = it.getVarMap_obj(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
						timeCoeff = astrohelion::sign(it.pArcIn->getTOFByIx(s))*2*it.X[tof_var.row0];
						break;
					case MSTOF_tp::VAR_FREE:
						tof_var = it.getVarMap_obj(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
						break;
					case MSTOF_tp::VAR_EQUALARC:
						tof_var = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
						timeCoeff = 1.0/(it.pArcIn->getNumSegs());
						break;
					default:
						throw Exception("DynamicsModel::multShoot_targetDeltaV: Unhandled time type");
				}

				Eigen::RowVectorXd dFdt_n = -1*dFdq_n2 * state_dot;
				it.DF_elements.push_back(Tripletd(row0, tof_var.row0, timeCoeff*dFdt_n(0)/dvMax));
			}
		}
	}
	
	// Copute the difference between the actual deltaV and the desired deltaV
	it.FX[row0] = con.getData()[0] == 0 ? totalDV : totalDV/dvMax - 1;

	if(con.getType() == Constraint_tp::MAX_DELTA_V){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it.slackAssignCon.begin(),
			it.slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it.totalFree - it.numSlack + (slackIx - it.slackAssignCon.begin());
		
		// printf("Constraint_tp::MAX_DELTA_V Constraint Details:\n");
		// printf("Total DV")
		it.FX[row0] += it.X[slackCol]*it.X[slackCol];
		it.DF_elements.push_back(Tripletd(row0, slackCol, 2*it.X[slackCol]));
	}
}//==============================================

/**
 *  @brief Compute the slack variable value for a delta-V constraint
 *  @details This function currently returns a hard-coded value of 1e-2
 * 
 *  @param it a reference to the class containing all the data relevant to the corrections process
 *	@param con the constraint being applied
 * 
 *  @return Ideal value of the slack variable
 */
double DynamicsModel::multShoot_targetDeltaV_compSlackVar(const MultShootData& it, const Constraint& con) const{
	// double totalDV = 0;
	// for(int n = 0; n < it.numNodes-1; n++){
	// 	// compute squared magnitude of DV between node n and n+1
	// 	// This takes the form v_n,f - v_n+1,0
	// 	totalDV += it.deltaVs[n*3]*it.deltaVs[n*3] +
	// 		it.deltaVs[n*3+1]*it.deltaVs[n*3+1] + 
	// 		it.deltaVs[n*3+2]*it.deltaVs[n*3+2];
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
 *	@brief Compute partials and constraint function values for time-of-flight constraints
 *
 *	This method *should* provide full functionality for any model; only 1's and 0's are
 *	used to relate TOFs.
 *	
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 *	
 *	@throws Exception if variable time is set to OFF
 */
void DynamicsModel::multShoot_targetTOF(MultShootData& it, const Constraint& con, int row0) const{
	if(to_underlying(it.tofTp) <= 0)
		throw Exception("DynamicsModel::multShoot_targetTOF: Cannot target TOF when variable time is off!");

	switch(it.tofTp){
		case MSTOF_tp::VAR_FREE:
		{
			// Sum all TOF for total, set partials w.r.t. integration times equal to one
			for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
				MSVarMap_Obj tof_var = it.getVarMap_obj(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
				it.FX[row0] += it.X[tof_var.row0];
				it.DF_elements.push_back(Tripletd(row0, tof_var.row0, astrohelion::sign(it.X[tof_var.row0])));
			}

			// subtract the desired TOF from the constraint to finish its computation
			it.FX[row0] -= con.getData()[0];
			break;
		}
		case MSTOF_tp::VAR_FIXSIGN:
		{
			// Sum all TOF for total, set partials w.r.t. integration times equal to one
			for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
				MSVarMap_Obj tof_var = it.getVarMap_obj(MSVar_tp::TOF, it.pArcIn->getSegRefByIx_const(s).getID());
				it.FX[row0] += astrohelion::sign(it.pArcIn->getTOFByIx(s)) * it.X[tof_var.row0] * it.X[tof_var.row0];
				it.DF_elements.push_back(Tripletd(row0, tof_var.row0, 2*it.X[tof_var.row0] * astrohelion::sign(it.pArcIn->getTOFByIx(s))));
			}

			// subtract the desired TOF from the constraint to finish its computation
			it.FX[row0] -= con.getData()[0];
			break;
		}
		case MSTOF_tp::VAR_EQUALARC:
		{
			MSVarMap_Obj tof_var = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
			it.FX[row0] = it.X[tof_var.row0];
			it.DF_elements.push_back(Tripletd(row0, tof_var.row0, 1.0));

			// subtract the desired TOF from the constraint to finish its computation
			it.FX[row0] -= con.getData()[0];
			break;
		}
		default:
			throw Exception("DynamicsModel::multShoot_targetTOF: Unhandled time type");
	}
}//====================================================

/**
 *	@brief Compute partials and constraint function values for apse constraints
 *
 *	This method *should* provide full functionality for any autonomous model. Non-
 *	autonomous models will need to modify the function to account for epoch time
 *	
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetApse(MultShootData& it, const Constraint& con, int row0) const{
	std::vector<double> conData = con.getData();
	MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, con.getID());
	int Pix = static_cast<int>(conData[0]);	// index of primary
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time
	
	if(state_var.row0 == -1)
		throw Exception("DynamicsModel::multShoot_targetApse: State vector is not part of free variable vector; cannot target apse.");

	// Get the primary position
	std::vector<double> primPos = getPrimPos(t, it.pArcIn->getSysData());

	// Get distance between node and primary in x, y, and z-coordinates, use non-scaled coordinates
	double dx = it.X[state_var.row0+0] - primPos[Pix*3+0];
	double dy = it.X[state_var.row0+1] - primPos[Pix*3+1];
	double dz = it.X[state_var.row0+2] - primPos[Pix*3+2];
	double vx = it.X[state_var.row0+3];
	double vy = it.X[state_var.row0+4];
	double vz = it.X[state_var.row0+5];

	// Constraint function: r_dot = 0 (using non-scaled coordinates)
	it.FX[row0] = dx*vx + dy*vy + dz*vz;

	// Partials of F w.r.t. node state
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+0, vx));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+1, vy));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+2, vz));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+3, dx));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+4, dy));
	it.DF_elements.push_back(Tripletd(row0, state_var.row0+5, dz));
}//====================================================

/**
 *	@brief Compute partials and constraint function values for apse constraints
 *	on segment ends
 *
 *	This method *should* provide full functionality for any autonomous model. Non-
 *	autonomous models will need to modify the function to account for epoch time
 *	
 *	@param it a reference to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param row0 the index of the row this constraint begins at
 */
void DynamicsModel::multShoot_targetApse_endSeg(MultShootData& it, const Constraint& con, int row0) const{
	std::vector<double> conData = con.getData();
	int segIx = it.pArcIn->getSegIx(con.getID());

	// Get object representing origin of segment
	MSVarMap_Obj prevNode_var = it.getVarMap_obj(MSVar_tp::STATE, it.pArcIn->getSegRef_const(con.getID()).getOrigin());

	MSVarMap_Obj tof_var;
	double timeCoeff = 1;
	if(to_underlying(it.tofTp) > 0){
		switch(it.tofTp){
			case MSTOF_tp::VAR_FREE:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF, con.getID());
				break;
			case MSTOF_tp::VAR_FIXSIGN:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF, con.getID());
				timeCoeff = astrohelion::sign(it.pArcIn->getTOF(con.getID()))*2*it.X[tof_var.row0];
				break;
			case MSTOF_tp::VAR_EQUALARC:
				tof_var = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
				timeCoeff = 1.0/(it.pArcIn->getNumSegs());
				break;
			default:
				throw Exception("DynamicsModel::multShoot_targetDist_endSeg: Unhandled time type");
		}
	}	

	// Data associated with the previous node and the propagated segment
	std::vector<double> lastState = it.propSegs[segIx].getStateByIx(-1);

	int Pix = static_cast<int>(conData[0]);	// index of primary
	double t = 0;	// If the system is non-autonomous, this will need to be replaced with an epoch time
	double primPos[3] = {0};
	getPrimPos(t, it.pArcIn->getSysData(), Pix, primPos);

	// Get distance between node and primary in x, y, and z-coordinates, use non-scaled coordinates
	double dx = lastState[0] - primPos[0];
	double dy = lastState[1] - primPos[1];
	double dz = lastState[2] - primPos[2];
	double vx = lastState[3];
	double vy = lastState[4];
	double vz = lastState[5];

	// Constraint function: r_dot = 0 (using non-scaled coordinates)
	it.FX[row0] = dx*vx + dy*vy + dz*vz;

	// Partials of F w.r.t. propagated state at segment end
	double dFdq_nf[6] = {vx, vy, vz, dx, dy, dz};
	
	if(prevNode_var.row0 != -1){
		MatrixXRd stm = it.propSegs[segIx].getSTMByIx(-1);

		// Do the matrix multiplication with loops to avoid expensive vector allocation
		double sum;
		for(unsigned int c = 0; c < 6; c++){
			sum = 0;
			for(unsigned int r = 0; r < 6; r++){
				sum += dFdq_nf[r]*stm(r,c);
			}
			it.DF_elements.push_back(Tripletd(row0, prevNode_var.row0+c, sum));
		}
	}

	// Partials of F w.r.t. time-of-flight
	if(to_underlying(it.tofTp) > 0){
		std::vector<double> lastDeriv = it.propSegs[segIx].getStateDerivByIx(-1);

		// Compute dot product between dFdq_nf and final state derivative vector
		double dp = 0;
		for(unsigned int r = 0; r < 6; r++){
			dp += dFdq_nf[r]*lastDeriv[r];
		}
		it.DF_elements.push_back(Tripletd(row0, tof_var.row0, timeCoeff*dp));
	}
}//====================================================

/**
 *  @brief Take the final, corrected free variable vector `X` and create an output 
 *  nodeset
 *
 *  @param it an iteration data object containing all info from the corrections process
 */
void DynamicsModel::multShoot_createOutput(const MultShootData& it) const{

    std::vector<int> newNodeIDs;
    newNodeIDs.reserve(it.numNodes);
    
    unsigned int n = 0, s = 0,i = 0;
    for(n = 0; n < it.numNodes; n++){
        Node node = it.pArcIn->getNodeByIx(n);

        MSVarMap_Obj state_var = it.getVarMap_obj(MSVar_tp::STATE, node.getID());
        if(state_var.row0 != -1){
            node.setState(&(it.X[state_var.row0]), state_var.nRows);
        }

        // If no control law is applied, getVarMap_obj will throw an exception
        try{
            MSVarMap_Obj ctrl_var = it.getVarMap_obj(MSVar_tp::CTRL, node.getID());

            if(ctrl_var.row0 != -1){
                std::vector<double> ctrl = std::vector<double>(it.X.begin()+ctrl_var.row0, it.X.begin() + ctrl_var.row0 + ctrl_var.nRows);
                node.setExtraParamVec(PARAMKEY_CTRL, ctrl);
            }
        }catch(const Exception &e){ /* No need to report exception */ }

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(it.pArcOut->addNode(node));
    }

    // Update any constraints that refer nodeIDs in their data vectors
    for(n = 0; n < it.numNodes; n++){
    	Node &newNode = it.pArcOut->getNodeRefByIx(n);
    	std::vector<Constraint>& nodeCons = newNode.getConsRef();
    	for(Constraint &con : nodeCons){
    		if(con.dataStoresID()){
    			std::vector<double> data = con.getData();
    			for(i = 0; i < data.size(); i++){
    				if(!std::isnan(data[i])){
    					data[i] = newNodeIDs[it.pArcIn->getNodeIx(data[i])];
    				}
    			}
    			con.setData(data);
    		}
    	}
    }

    double tof;
    int newOrigID, newTermID;
    for(s = 0; s < it.pArcIn->getNumSegs(); s++){
        Segment seg = it.pArcIn->getSegByIx(s);

        if(to_underlying(it.tofTp) > 0){
        	switch(it.tofTp){
        		case MSTOF_tp::VAR_FREE:
        		{
        			MSVarMap_Obj tofVar = it.getVarMap_obj(MSVar_tp::TOF, seg.getID());
        			tof = it.X[tofVar.row0];
        			break;
        		}
        		case MSTOF_tp::VAR_FIXSIGN:
        		{
        			MSVarMap_Obj tofVar = it.getVarMap_obj(MSVar_tp::TOF, seg.getID());
        			tof = astrohelion::sign(it.pArcIn->getTOFByIx(s)) * (it.X[tofVar.row0])*(it.X[tofVar.row0]);
        			break;	
        		}
        		case MSTOF_tp::VAR_EQUALARC:
        		{
        			MSVarMap_Obj tofVar = it.getVarMap_obj(MSVar_tp::TOF_TOTAL, Linkable::INVALID_ID);
        			tof = it.X[tofVar.row0]/it.pArcIn->getNumSegs();
        			break;
        		}
        		default:
					throw Exception("DynamicsModel::multShoot_createOutput: Unhandled time type");
        	}
        }else{
            tof = seg.getTOF();
        }

        newOrigID = newNodeIDs[it.pArcIn->getNodeIx(seg.getOrigin())];
        int termID = seg.getTerminus();
        newTermID = termID == Linkable::INVALID_ID ? termID : newNodeIDs[it.pArcIn->getNodeIx(termID)];
        
        Segment newSeg(newOrigID, newTermID, tof);
        newSeg.setConstraints(seg.getConstraints());
        newSeg.setVelCon(seg.getVelCon());
        newSeg.setSTM(it.propSegs[s].getSTMByIx(-1));
        newSeg.setCtrlLaw(seg.getCtrlLaw());
        newSeg.setStateVector(it.propSegs[s].getSegRef_const(0).getStateVector());
        newSeg.setStateWidth(it.propSegs[s].getSegRef_const(0).getStateWidth());
        newSeg.setTimeVector(it.propSegs[s].getSegRef_const(0).getTimeVector());
        it.pArcOut->addSeg(newSeg);
    }

    // Determine the chronological order of the nodeset
    // it.pArcOut->print();
    // it.pArcOut->printInChrono();
    std::vector<ArcPiece> order = it.pArcOut->getChronoOrder();
    // Set the epoch of each node based on the time of flight from
    // the first node
    double epoch = NAN;
    for(i = 0; i < order.size(); i++){
        if(order[i].type == ArcPiece::Piece_tp::NODE){
            if(std::isnan(epoch)){
                // Copy the epoch value of the first node
                epoch = it.pArcOut->getNode(order[i].id).getEpoch();
            }else{
                // Set the epoch value of all other nodes
                it.pArcOut->getNodeRef(order[i].id).setEpoch(epoch);       
            }
        }
        if(order[i].type == ArcPiece::Piece_tp::SEG){
        	// Retrieve the first time in a chronological sense
        	double firstTime = it.pArcOut->getSegRef_const(order[i].id).getTOF() > 0 ? 
        		it.pArcOut->getSegRef_const(order[i].id).getTimeByIx(0) :
        		it.pArcOut->getSegRef_const(order[i].id).getTimeByIx(-1);

        	// If arcset begins with a segment, epoch will be NAN at first
        	if(std::isnan(epoch))
        		epoch = firstTime;

        	// Shift all segment times so that the first chronological time is consistent
        	// with the epoch
        	it.pArcOut->getSegRef(order[i].id).shiftAllTimes(epoch - firstTime);

            if(!std::isnan(epoch)){
                // When stepping through in chronological order, every step is
                // forward in time; negative TOFs are associated with segments that
                // flow opposite the chronological order; ignore sign here.
                epoch += std::abs(it.pArcOut->getSeg(order[i].id).getTOF());
            }
        }
    }

    // it.pArcOut->print();
    std::vector<Constraint> arcCons = it.pArcIn->getArcConstraints();
    for(i = 0; i < arcCons.size(); i++){
        it.pArcOut->addConstraint(arcCons[i]);
    }
}//======================================================

/**
 * @brief Check the analytically-computed partial derivatives in the A matrix 
 * by comparing with a central-differencing method
 * 
 * @param q Core and control state vector; must have enough elements to specify
 * all core and control states with the core size from the DynamicsModel and
 * the control size from the control law in the EOM_ParamStruct
 * @param t time at which to evaluate the equations of motion
 * @param pertSize size of the perturbation to use in central differencing
 * @param params structure containing information for the equations of motion,
 * including the system data and control law
 * @param verb Verbosity, defaults to Verbosity_tp::SOME_MSG
 * @param writeToFile Whether or not to save the computing A matrices to file;
 * they are saved to the current directory with filenames following the pattern
 * "FiniteDiff_A*.csv"
 * 
 * @return true if the numeric and analytic A matrices match to within a 
 * "reasonable" tolerance
 */
bool DynamicsModel::finiteDiff_checkAMat(double *q, double t, double pertSize, 
	EOM_ParamStruct *params, Verbosity_tp verb, bool writeToFile) const{

	if(params == nullptr){
		throw Exception("DynamicsModel::finiteDiff_checkAMat: params is null");
	}

	unsigned int ctrl_dim = 0;
	if(params->pCtrlLaw){
		ctrl_dim = params->pCtrlLaw->getNumStates();
	}

	unsigned int q_dim = coreDim + ctrl_dim;
	unsigned int vecLen = q_dim*(q_dim+1) + extraDim;
	std::vector<double> q0(vecLen, 0);
	
	// Assign arbitrary state values to state and control; extra variables left 
	// as zeros
	for(unsigned int i = 0; i < q_dim; i++)
		q0[i] = q[i];

	// Assign STM = Identity
	for(unsigned int i = q_dim; i < q_dim*(q_dim+1); i += q_dim+1)
		q0[i] = 1;

	MatrixXRd A_analytic(q_dim, q_dim), A_numeric(q_dim, q_dim);

	// Get the A matrix analytically
	std::vector<double> dq(vecLen, 0), qp(vecLen, 0), qdot(vecLen, 0);
	
	eom_fcn eoms = getFullEOM_fcn();

	eoms(t, &(q0[0]), &(dq[0]), params);
	A_analytic = Eigen::Map<MatrixXRd>(&(dq[q_dim]), q_dim, q_dim);

	for(unsigned int i = 0; i < q_dim; i++){
		qdot.assign(vecLen, 0);

		// Perturb the state vector in forward direction
		qp = q0;
		qp[i] += pertSize;

		// Compute the state derivative
		eoms(t, &(qp[0]), &(qdot[0]), params);
		Eigen::VectorXd Af = Eigen::Map<Eigen::VectorXd>(&(qdot[0]), q_dim, 1);

		// Perturb the state vector in the backward direction
		qp = q0;
		qdot.assign(vecLen, 0);
		qp[i] -= pertSize;
		
		// Compute state derivatives
		eoms(t, &(qp[0]), &(qdot[0]), params);
		Eigen::VectorXd Ab = Eigen::Map<Eigen::VectorXd>(&(qdot[0]), q_dim, 1);

		A_numeric.col(i) = (Af - Ab)/(2*pertSize);
	}

	MatrixXRd diff = A_analytic - A_numeric;
	diff = diff.cwiseAbs();	// Get coefficient-wise aboslute value

	// Relative differences (divide by magnitude of each value)
	MatrixXRd relDiff = MatrixXRd::Zero(diff.rows(), diff.cols());	

    // Divide each element by the magnitude of the DF element to get a relative 
    // difference magnitude
    for(int r = 0; r < diff.rows(); r++){
        for(int c = 0; c < diff.cols(); c++){
            // If one of the elements is zero, let the difference just be the 
            // difference; no division
            if(std::abs(A_analytic(r,c)) > 1e-13 && std::abs(A_numeric(r,c)) > 1e-13)
                relDiff(r,c) = diff(r,c)/std::abs(A_numeric(r,c));
        }
    }

	// Compute the largest coefficient in each row and column
    std::vector<unsigned int> rowMaxIndex(diff.rows(), 0), colMaxIndex(diff.cols(), 0);
    Eigen::VectorXd rowMax(diff.rows(), 1);
    Eigen::VectorXd colMax(diff.cols(), 1);
    for(unsigned int r = 0; r < diff.rows(); r++){
    	rowMax(r) = diff.row(r).maxCoeff(&rowMaxIndex[r]);
    }
	for(unsigned int c = 0; c < diff.cols(); c++){
		colMax(c) = diff.col(c).maxCoeff(&colMaxIndex[c]);
	}

	unsigned int rowMaxMaxIx = 0, colMaxMaxIx = 0;
    double rowMaxMax = rowMax.maxCoeff(&rowMaxMaxIx);
    double colMaxMax = colMax.maxCoeff(&colMaxMaxIx);
    int errScalar = 100000;

    bool goodA = true;

    if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
        if(verb >= Verbosity_tp::SOME_MSG)
        printVerbColor(verb >= Verbosity_tp::SOME_MSG, BOLDGREEN,
        	"No significant errors! (Abs diff is small)\n");

        goodA = goodA && true;
    }else{

    	if(rowMaxMax >= errScalar*pertSize){
    		// row error is too big
    		if(relDiff(rowMaxMaxIx, rowMaxIndex[rowMaxMaxIx]) < errScalar*pertSize){
    			// The relative difference is small enough
    			printVerbColor(verb >= Verbosity_tp::SOME_MSG, BOLDGREEN,
    				"No significant errors! (Abs diff is large, rel diff is small)\n");
    			goodA = goodA && true;
    		}else{
    			printColor(BOLDRED, "Significant errors (Abs diff and Rel diff)\n");
    			goodA = goodA && false;
    		}
    	}else{
    		// Column error is too big
    		if(relDiff(colMaxIndex[colMaxMaxIx], colMaxMaxIx) < errScalar*pertSize){
    			printVerbColor(verb >= Verbosity_tp::SOME_MSG, BOLDGREEN,
    				"No significant errors! (Abs diff is large, rel diff is small)\n");
    			goodA = goodA && true;
    		}else{
    			printColor(BOLDRED, "Significant errors (Abs diff and Rel diff)\n");
    			goodA = goodA && false;
    		}
    	}
    }

	if(verb >= Verbosity_tp::ALL_MSG){
		printColor(BLUE, "Analytic:\n");
		std::cout << A_analytic << std::endl;
		printColor(BLUE, "Numeric:\n");
		std::cout << A_numeric << std::endl;
		printColor(BLUE, "Difference:\n");
		std::cout << diff << std::endl;
	}

	if(writeToFile){
        toCSV(A_numeric, "FiniteDiff_Aest.csv");
        toCSV(A_analytic, "FiniteDiff_A.csv"); 
        toCSV(diff, "FiniteDiff_ADiff.csv");
    }

    return goodA;
}//======================================================

}// END of Astrohelion namespace