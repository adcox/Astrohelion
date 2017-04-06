/**
 *  \file DynamicsModel_bc4bp.cpp
 *  \brief Derivative of DynamicsModel, specific to BCR4BPR
 *  
 *  \author Andrew Cox
 *  \version May 25, 2016
 *  \copyright GNU GPL v3.0
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

#include "DynamicsModel_bc4bp.hpp"

#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "CorrectionEngine.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Nodeset_bc4bp.hpp"
#include "SysData_bc4bp.hpp"
#include "Traj_bc4bp.hpp"
#include "Event.hpp"
#include "Node.hpp"
#include "Utilities.hpp"

#include <cspice/SpiceUsr.h>
#include <gsl/gsl_errno.h>

namespace astrohelion{
/**
 *  \brief Construct a BCR4BP Dynamic DynamicsModel
 */
DynamicsModel_bc4bp::DynamicsModel_bc4bp() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP) {
    coreStates = 6;
    extraStates = 6;
    allowedCons.push_back(Constraint_tp::SP);
    allowedCons.push_back(Constraint_tp::SP_RANGE);
    allowedCons.push_back(Constraint_tp::SP_DIST);
    allowedCons.push_back(Constraint_tp::SP_MAX_DIST);
}//==============================================

/**
 *  \brief Copy constructor
 *  \param m a model reference
 */
DynamicsModel_bc4bp::DynamicsModel_bc4bp(const DynamicsModel_bc4bp &m) : DynamicsModel(m) {}

/**
 *  \brief Assignment operator
 *  \param m a model reference
 */
DynamicsModel_bc4bp& DynamicsModel_bc4bp::operator =(const DynamicsModel_bc4bp &m){
	DynamicsModel::operator =(m);
	return *this;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_bc4bp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_bc4bp::getFullEOM_fcn() const{
	return &fullEOMs;
}//==============================================

/**
 *  \brief Compute the positions of all primaries
 *
 *  \param t the epoch at which the computations occur
 *  \param sysData object describing the specific system
 *  \return an n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> DynamicsModel_bc4bp::getPrimPos(double t, const SysData *sysData) const{
    double primPos[9];
    const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(sysData);
    getPrimaryPos(t, bcSys, primPos);

    return std::vector<double>(primPos, primPos+9);
}//==============================================

/**
 *  \brief Compute the velocities of all primaries
 *
 *  \param t the epoch at which the computations occur
 *  \param sysData object describing the specific system
 *  \return an n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one velocity vector in non-dimensional units
 */
std::vector<double> DynamicsModel_bc4bp::getPrimVel(double t, const SysData *sysData) const{
    double primVel[9];
    const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(sysData);
    getPrimaryVel(t, bcSys, primVel);

    return std::vector<double>(primVel, primVel+9);
}//==============================================

/**
 *  \brief Retrieve the state derivative
 *  \details Evaluate the equations of motion to compute the state time-derivative at 
 *  the specified time and state
 * 
 *  \param t time parameter
 *  \param state state vector
 *  \param params structure containing parameters relevant to the integration
 *  \return the time-derivative of the state vector
 */
std::vector<double> DynamicsModel_bc4bp::getStateDeriv(double t, std::vector<double> state, EOM_ParamStruct *params) const{
    if(state.size() != coreStates)
        throw Exception("DynamicsModel_bc4bp::getStateDeriv: State size does not match the core state size specified by the dynamical model");

    // Compute the acceleration
    std::vector<double> dsdt(coreStates,0);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), params);
    
    return dsdt;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Takes an input state and time and saves the data to the trajectory
 *  \param y an array containing the core state and any extra states integrated
 *  by the EOM function, including STM elements.
 *  \param t the time at the current integration state
 *  \param traj a pointer to the trajectory we should store the data in
 *  \param params structure containing parameters relevant to the integration
 */
void DynamicsModel_bc4bp::sim_saveIntegratedData(const double* y, double t, Traj* traj, EOM_ParamStruct *params) const{

    DynamicsModel::sim_saveIntegratedData(y, t, traj, params);
    
    Traj_bc4bp *bcTraj = static_cast<Traj_bc4bp*>(traj);
    bcTraj->set_dqdTByIx(-1, y+42); // dqdT is stored in y(42:47)
}//=====================================================

/**
 *  \brief Use a correction algorithm to accurately locate an event crossing
 *
 *  The simulation engine calls this function if and when it determines that an event 
 *  has been crossed. To accurately locate the event, we employ differential corrections
 *  and find the exact event occurence in space and time.
 *
 *  \param event the event we're looking for
 *  \param traj a pointer to the trajectory the event should occur on
 *  \param ic the core state vector for this system
 *  \param t0 non-dimensional time at the beginning of the search arc
 *  \param tof the time-of-flight for the arc to search over
 *  \param params structure containing parameters relevant to the integration
 *  \param verbose whether or not we should be verbose with output messages
 *
 *  \return wether or not the event has been located. If it has, a new point
 *  has been appended to the trajectory's data vectors.
 */
bool DynamicsModel_bc4bp::sim_locateEvent(Event event, Traj *traj,
    const double *ic, double t0, double tof, EOM_ParamStruct *params, Verbosity_tp verbose) const{

    // **** Make sure you fix the epoch of the first node as well as the states
    double IC[7] = {0};
    std::copy(ic, ic+6, IC);
    IC[6] = t0;

    // Recast system data pointer
    const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp*>(params->sysData);

    // Create a nodeset for this particular type of system
    astrohelion::printVerb(verbose >= Verbosity_tp::ALL_MSG, "  Creating nodeset for event location\n");
    Nodeset_bc4bp eventNodeset(bcSys, IC, t0, tof, 2, Nodeset::TIME);   //bc4bp nodeset does not use control laws

    // Constraint to keep first node unchanged
    Constraint fixFirstCon(Constraint_tp::STATE, 0, IC, 7);

    // Constraint to enforce event
    Constraint eventCon(event.getConType(), 1, event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose == Verbosity_tp::ALL_MSG){ eventNodeset.print(); }

    astrohelion::printVerb(verbose >= Verbosity_tp::ALL_MSG, "  Applying corrections process to locate event\n");
    CorrectionEngine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbosity(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations
    
    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    Nodeset_bc4bp correctedNodes(bcSys);
    try{
        corrector.multShoot(&eventNodeset, &correctedNodes);
    }catch(DivergeException &e){
        if(verbose >= Verbosity_tp::SOME_MSG)
            astrohelion::printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(LinAlgException &e){
        if(verbose >= Verbosity_tp::SOME_MSG)
            astrohelion::printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    std::vector<double> state = correctedNodes.getStateByIx(-1);
    std::vector<double> stm = correctedNodes.getExtraParamVecByIx(-1, PARAMKEY_STM);
    std::vector<double> dqdT = correctedNodes.getExtraParamVecByIx(-1, PARAMKEY_STATE_EPOCH_DERIV);
    
    state.insert(state.end(), stm.begin(), stm.end());
    state.insert(state.end(), dqdT.begin(), dqdT.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOFByIx(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    sim_saveIntegratedData(&(state[0]), eventTime, traj, params);
    
    return true;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Initialize the corrector's design vector with position and velocity states,
 *  and times-of-flight.
 *
 *  Derived models may replace this function or call it and then append more design 
 *  variables.
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param set a pointer to the nodeset being corrected
 *  \throws Exception if equal arc times is turned ON; this has not been implemented for this system
 */
void DynamicsModel_bc4bp::multShoot_initDesignVec(MultShootData *it, const Nodeset *set) const{
    // Call base class to do most of the work
    DynamicsModel::multShoot_initDesignVec(it, set);

    // Append the Epoch for each node
    if(it->bVarTime){
        // epochs come after ALL the TOFs have been added
        const Nodeset_bc4bp *bcSet = static_cast<const Nodeset_bc4bp *>(set);
        for(int n = 0; n < bcSet->getNumNodes(); n++){
            MSVarMap_Key key(MSVar_tp::EPOCH, set->getNodeByIx(n).getID());
            it->freeVarMap[key] = MSVarMap_Obj(key, it->X.size());
            it->X.push_back(bcSet->getEpochByIx(n));
        }
    }
}//====================================================

/**
 *  \brief Perform model-specific initializations on the MultShootData object
 *  \param it pointer to the object to be initialized
 */
void DynamicsModel_bc4bp::multShoot_initIterData(MultShootData *it) const{
    Traj_bc4bp traj(static_cast<const SysData_bc4bp *>(it->sysData));
    it->propSegs.assign(it->nodeset->getNumSegs(), traj);
}//====================================================

/**
 *  \brief Create continuity constraints for the correction algorithm; this function
 *  creates position and velocity constraints.
 *
 *  This function overrides the base model's to add time continuity
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param set a pointer to the nodeset being corrected
 */ 
void DynamicsModel_bc4bp::multShoot_createContCons(MultShootData *it, const Nodeset *set) const{
    DynamicsModel::multShoot_createContCons(it, set);

    if(it->bVarTime){
        std::vector<double> zero {0};
        for(int s = 0; s < it->nodeset->getNumSegs(); s++){
            if(it->nodeset->getSegByIx(s).getTerminus() != Linkable::INVALID_ID){
                Constraint timeCont(Constraint_tp::CONT_EX, it->nodeset->getSegByIx(s).getID(), zero);   // 0 value is for Epoch Time
                it->allCons.push_back(timeCont);
            }
        }
    }
}//============================================================

/**
 *  \brief Retrieve the initial conditions for a segment that the correction
 *  engine will integrate.
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param set a pointer to the nodeset being corrected
 *  \param s the ID of the segment being propagated
 *  \param ic a pointer to a 6-element initial state array
 *  \param t0 a pointer to a double representing the initial time (epoch)
 *  \param tof a pointer to a double the time-of-flight on the segment.
 */
void DynamicsModel_bc4bp::multShoot_getSimICs(const MultShootData *it, const Nodeset *set, int s,
    double *ic, double *t0, double *tof) const{

    DynamicsModel::multShoot_getSimICs(it, set, s, ic, t0, tof);   // Perform default behavior

    // Compute and reverse-scale epoch
    const Nodeset_bc4bp *bcSet = static_cast<const Nodeset_bc4bp *>(set);

    if(it->bVarTime){
        MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSeg(s).getOrigin());
        *t0 = it->X[epochVar.row0];
    }else{
        *t0 = bcSet->getEpoch(s);
    }
}//============================================================

/**
 *  \brief Compute the value of a slack variable for an inequality constraint.
 *  \details Computing the value of the slack variable can avoid unneccessary 
 *  shooting iterations when the inequality constraint is already met. If the 
 *  inequality constraint is met, the value returned by this function will make
 *  the constraint function evaluate to zero.
 *  
 *  This method adds additional functionality for the BCR4BP and calls the base
 *  model version for most constraint types
 *  
 *  Note: This function should be called after the state variable vector has 
 *  been initialized by the multiple shooting algorithm
 * 
 *  \param it the MultShootData object associated with the multiple shooting process
 *  \param con the inequality constraint for which the slack variable is being computed
 * 
 *  \return The value of the slack variable that minimizes the constraint function
 *  without setting the slack variable to zero
 */
double DynamicsModel_bc4bp::multShoot_getSlackVarVal(const MultShootData *it, Constraint con) const{
    switch(con.getType()){
        case Constraint_tp::SP_RANGE:
            return multShoot_targetSPMag_compSlackVar(it, con);
        case Constraint_tp::SP_MAX_DIST:
            return multShoot_targetSP_maxDist_compSlackVar(it, con);
        default:
            return DynamicsModel::multShoot_getSlackVarVal(it, con);
    }
}//===========================================================

/**
 *  \brief Compute constraint function and partial derivative values for a constraint
 *  
 *  This function calls its relative in the DynamicsModel base class and appends additional
 *  instructions specific to the BCR4BPR
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param con the constraint being applied
 *  \param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void DynamicsModel_bc4bp::multShoot_applyConstraint(MultShootData *it, Constraint con, int c) const{

    // Let the base class do its thing first
    DynamicsModel::multShoot_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case Constraint_tp::APSE:
            multShoot_targetApse(it, con, row0);
            break;
        case Constraint_tp::SP:
            multShoot_targetSP(it, con, row0);
            break;
        case Constraint_tp::SP_RANGE:
            multShoot_targetSP_mag(it, con, c);
            break;
        case Constraint_tp::SP_DIST:
        case Constraint_tp::SP_MAX_DIST:
            multShoot_targetSP_dist(it, con, c);
            break;
        default: break;
    }
}//=========================================================

/**
 *  \brief Compute position and velocity constraint values and partial derivatives
 *
 *  This function computes and stores the default position continuity constraints as well
 *  as velocity constraints for all nodes marked continuous in velocity. The delta-Vs
 *  between arc segments and node states are recorded, and the partial derivatives of each
 *  node with respect to other nodes and integration time are all computed
 *  and placed in the appropriate spots in the Jacobian matrix.
 *
 *  This function replaces the one found in the base model
 *
 *  \param it a pointer to the correctors iteration data structure
 *  \param con the constraint being applied
 *  \param row0 the first row this constraint applies to
 */
void DynamicsModel_bc4bp::multShoot_targetCont_State(MultShootData* it, Constraint con, int row0) const{
    // Call base function first to do most of the work
    DynamicsModel::multShoot_targetCont_State(it, con, row0);

    // Add epoch dependencies for this model
    if(it->bVarTime){
        int segIx = it->nodeset->getSegIx(con.getID());
        std::vector<double> conData = con.getData();
        std::vector<double> last_dqdT = it->propSegs[segIx].getExtraParamVecByIx(-1, PARAMKEY_STATE_EPOCH_DERIV);

        MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(segIx).getOrigin());

        // Loop through conData
        int count = 0;
        for(unsigned int s = 0; s < conData.size(); s++){
            if(!std::isnan(conData[s])){
                // Epoch dependencies
                it->DF[it->totalFree*(row0+count) + epochVar.row0] = last_dqdT[s];
                count++;
            }
        }
    }
}//=========================================================

/**
 *  \brief Computes continuity constraints for constraints with the <tt>CONT_EX</tt> type.
 *
 *  This function overrides the base model function
 *
 *  \param it a pointer to the correctors iteration data structure
 *  \param con the constraint being applied
 *  \param row0 the first row this constraint applies to
 */
void DynamicsModel_bc4bp::multShoot_targetCont_Ex(MultShootData *it, Constraint con, int row0) const{
    /* Add time-continuity constraints if applicable; we need to match
    the epoch time of node n to the sum of node n-1's epoch and TOF */
    if(it->bVarTime){
        int segIx = it->nodeset->getSegIx(con.getID());
        MSVarMap_Obj T0_var = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(segIx).getOrigin()); 
        MSVarMap_Obj Tf_var = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(segIx).getTerminus());
        MSVarMap_Obj tof_var = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
            it->bEqualArcTime ? Linkable::INVALID_ID : con.getID());
        
        double T0 = it->X[T0_var.row0];
        double tof = it->bEqualArcTime ? it->X[tof_var.row0]/(it->nodeset->getNumSegs()) : it->X[tof_var.row0];
        double T1 = it->X[Tf_var.row0];

        it->FX[row0] = T1 - (T0 + tof);
        it->DF[it->totalFree*(row0) + tof_var.row0] = -1;
        it->DF[it->totalFree*(row0) + tof_var.row0] /= it->bEqualArcTime ? it->nodeset->getNumSegs() : 1.0;
        it->DF[it->totalFree*(row0) + T0_var.row0] = -1;
        it->DF[it->totalFree*(row0) + Tf_var.row0] = 1;
    }
}//=========================================================

/**
 *  \brief Computes continuity constraints for constraints with the <tt>Constraint_tp::SEG_CONT_EX</tt> type.
 *
 *  This function overrides the base model function
 *
 *  \param it a pointer to the correctors iteration data structure
 *  \param con the constraint being applied
 *  \param row0 the first row this constraint applies to
 */
void DynamicsModel_bc4bp::multShoot_targetCont_Ex_Seg(MultShootData *it, Constraint con, int row0) const{
    if(it->bVarTime){
        int segIx1 = it->nodeset->getSegIx(con.getID());
        int segIx2 = it->nodeset->getSegIx(con.getData()[0]);

        MSVarMap_Obj T0_var1 = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(segIx1).getOrigin());
        MSVarMap_Obj T0_var2 = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(segIx2).getOrigin());
        MSVarMap_Obj tof_var1 = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
            it->bEqualArcTime ? Linkable::INVALID_ID : con.getID());
        MSVarMap_Obj tof_var2 = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
            it->bEqualArcTime ? Linkable::INVALID_ID : con.getData()[0]);

        double T01 = it->X[T0_var1.row0];
        double T02 = it->X[T0_var2.row0];
        double tof1 = it->bEqualArcTime ? it->X[tof_var1.row0]/(it->nodeset->getNumSegs()) : it->X[tof_var1.row0];
        double tof2 = it->bEqualArcTime ? it->X[tof_var2.row0]/(it->nodeset->getNumSegs()) : it->X[tof_var2.row0];

        it->FX[row0] = T01 + tof1 - (T02 + tof2);
        it->DF[it->totalFree*row0 + T0_var1.row0] = 1;
        it->DF[it->totalFree*row0 + T0_var2.row0] = -1;
        it->DF[it->totalFree*row0 + tof_var1.row0] = 1;
        it->DF[it->totalFree*row0 + tof_var2.row0] = -1;
    }
}//=========================================================

/**
 *  \brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::STATE</tt>.
 *
 *  This method replaces the base class function and allows the user to constrain epoch as well
 *  as the configuration space states
 *
 *  \param it a pointer to the class containing all the data relevant to the corrections process
 *  \param con the constraint being applied
 *  \param row0 the index of the row this constraint begins at
 */
void DynamicsModel_bc4bp::multShoot_targetState(MultShootData* it, Constraint con, int row0) const{
    std::vector<double> conData = con.getData();
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());

    // Allow user to constrain 7 states
    int count = 0;  // Count # rows since some may be skipped (NAN)
    for(unsigned int s = 0; s < con.getData().size(); s++){
        if(!std::isnan(conData[s])){
            if(s < 6){
                it->FX[row0+count] = it->X[stateVar.row0+s] - conData[s];
                it->DF[it->totalFree*(row0 + count) + stateVar.row0 + s] = 1;
                count++;
            }else if(s == 6){
                // Allow constraining epoch
                if(it->bVarTime){
                    MSVarMap_Obj epoch_var = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
                    it->FX[row0+count] = it->X[epoch_var.row0] - conData[s];
                    it->DF[it->totalFree*(row0 + count) + epoch_var.row0] = 1;
                    count++;
                }
            }else{
                astrohelion::printErr("DynamicsModel_bc4bp::multShoot_targetState: constraint has more than six elements... ignoring the extras\n");
            }
        }
    }
}//=================================================

/**
 *  \brief Compute partials and constraint functions for nodes constrained with <tt>Constraint_tp::DIST</tt>, 
 *  <tt>Constraint_tp::MIN_DIST</tt>, or <tt>Constraint_tp::MAX_DIST</tt>
 *
 *  This method overrides the base class function to add functionality for epoch-time dependencies
 *
 *  \param it a pointer to the class containing all the data relevant to the corrections process
 *  \param con a copy of the constraint object
 *  \param c the index of this constraint in the constraint vector object
 */
void DynamicsModel_bc4bp::multShoot_targetDist(MultShootData* it, Constraint con, int c) const{

    std::vector<double> conData = con.getData();
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    int Pix = static_cast<int>(conData[0]);    // index of primary
    int row0 = it->conRows[c];

    // Get the node epoch either from the design vector or from the original set of nodes
    MSVarMap_Obj epochVar;
    double t0 = 0;
    if(it->bVarTime){
        epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        t0 = it->X[epochVar.row0];
    }else{
        t0 = it->nodeset->getEpoch(con.getID());
    }
    // MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
    // double t0 = it->bVarTime ? it->X[epochVar.row0]/sT : it->nodeset->getEpoch(con.getID());

    // Get the primary position
    std::vector<double> primPos = getPrimPos(t0, it->sysData);

    // Get distance between node and primary in x, y, and z-coordinates
    double dx = it->X[stateVar.row0+0] - primPos[Pix*3+0];
    double dy = it->X[stateVar.row0+1] - primPos[Pix*3+1];
    double dz = it->X[stateVar.row0+2] - primPos[Pix*3+2];

    double h = sqrt(dx*dx + dy*dy + dz*dz);     // true distance

    // Compute difference between desired distance and true distance
    it->FX[row0] = h - conData[1];

    // Partials with respect to node position states
    it->DF[it->totalFree*row0 + stateVar.row0 + 0] = dx/h;
    it->DF[it->totalFree*row0 + stateVar.row0 + 1] = dy/h;
    it->DF[it->totalFree*row0 + stateVar.row0 + 2] = dz/h;

    if(it->bVarTime){
        // Epoch dependencies from primary positions
        double dhdr_data[] = {-dx/h, -dy/h, -dz/h};
        std::vector<double> primVel = getPrimVel(t0, it->sysData);
    
        Eigen::RowVector3d dhdr = Eigen::Map<Eigen::RowVector3d>(dhdr_data, 1, 3);
        Eigen::Vector3d drdT = Eigen::Map<Eigen::Vector3d>(&(primVel[Pix*3]), 3, 1);
        
        double prod = dhdr*drdT;
        it->DF[it->totalFree*row0 + epochVar.row0] = prod;
    }

    // Extra stuff for inequality constraints
    if(con.getType() == Constraint_tp::MIN_DIST || 
        con.getType() == Constraint_tp::MAX_DIST ){
        // figure out which of the slack variables correspond to this constraint
        std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
            it->slackAssignCon.end(), c);

        // which column of the DF matrix the slack variable is in
        int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
        int sign = con.getType() == Constraint_tp::MAX_DIST ? 1 : -1;

        // printf("Dist from P%d is %f (%s %f)\n", Pix, h,
        //     con.getType() == Constraint_tp::MIN_DIST ? "Min" : "Max", conData[1]);
        // printf("  Slack Var^2 = %e\n", it->X[slackCol]*it->X[slackCol]);
        // Subtract squared slack variable from constraint
        it->FX[row0] += sign*it->X[slackCol]*it->X[slackCol];

        // Partial with respect to slack variable
        it->DF[it->totalFree*row0 + slackCol] = sign*2*it->X[slackCol];
    }
}// End of targetDist() =========================================

/**
 *  \brief  Compute the value of the slack variable for inequality distance constraints
 *  \details This function computes a value for the slack variable in an
 *  inequality distance constraint. If the constraint is already met by the initial
 *  design, using this value will prevent the multiple shooting algorithm from
 *  searching all over for the propper value.
 *  
 *  Note: This method overrides the base class function to add functionality for non-zero epochs
 * 
 *  \param it the iteration data object for the multiple shooting process
 *  \param con the constraint the slack variable applies to
 *  \return the value of the slack variable that minimizes the constraint function
 *  without setting the slack variable equal to zero
 */
double DynamicsModel_bc4bp::multShoot_targetDist_compSlackVar(const MultShootData* it, Constraint con) const{
    std::vector<double> conData = con.getData();
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    int Pix = static_cast<int>(conData[0]);    // index of primary

    // Get the node epoch either from the design vector or from the original set of nodes
    double t0 = 0;
    if(it->bVarTime){
        MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        t0 = it->X[epochVar.row0];
    }else{
        t0 = it->nodeset->getEpoch(con.getID());
    }

    // Get the primary position
    std::vector<double> primPos = getPrimPos(t0, it->sysData);

    // Get distance between node and primary in x, y, and z-coordinates
    double dx = it->X[stateVar.row0+0] - primPos[Pix*3+0];
    double dy = it->X[stateVar.row0+1] - primPos[Pix*3+1];
    double dz = it->X[stateVar.row0+2] - primPos[Pix*3+2];

    double h = sqrt(dx*dx + dy*dy + dz*dz);     // true distance
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
 *  \brief Compute partials and constraints for all nodes constrained with <tt>Constraint_tp::DELTA_V</tt> or
 *  <tt>MIN_DELTA_V</tt>
 *
 *  Because the delta-V constraint applies to the entire trajectory, the constraint function values
 *  and partial derivatives must be computed for each node along the trajectory. This function
 *  takes care of all of them at once.
 *
 *  This function overrides the base targeting function to add support for epoch dependencies
 *
 *  \param it a pointer to the class containing all the data relevant to the corrections process
 *  \param con the constraint being applied
 *  \param c the index of the first row for this constraint
 */
void DynamicsModel_bc4bp::multShoot_targetDeltaV(MultShootData* it, Constraint con, int c) const{
    // Call base function to take care of most of the constraint computations and partials
    DynamicsModel::multShoot_targetDeltaV(it, con, c);

    if(it->bVarTime){
        // Add partials w.r.t. epoch time
        int row0 = it->conRows[c];

        // Don't allow dividing by zero, but otherwise scale by value to keep order ~1
        double dvMax = con.getData()[0] == 0 ? 1 : con.getData()[0];

        for(int s = 0; s < it->nodeset->getNumSegs(); s++){
            // compute magnitude of DV between segment s and its terminal point
            // This takes the form v_n,f - v_n+1,0
            double dvx = it->deltaVs[s*3+0];
            double dvy = it->deltaVs[s*3+1];
            double dvz = it->deltaVs[s*3+2];
            double dvMag = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

            // Don't bother computing partials if this dv is zero; all partials will equal zero
            // but the computation would divide by zero too, which would be... unfortunate.
            if(dvMag > 0){
                // Compute parial w.r.t. node n+1 (where velocity is discontinuous)
                double dFdq_n2_data[] = {0, 0, 0, -dvx/dvMag, -dvy/dvMag, -dvz/dvMag};
                Eigen::RowVectorXd dFdq_n2 = Eigen::Map<Eigen::RowVectorXd>(dFdq_n2_data, 1, 6);

                // Compute partial w.r.t. epoch time n
                std::vector<double> last_dqdT = it->propSegs.at(s).getExtraParamVecByIx(-1, PARAMKEY_STATE_EPOCH_DERIV);
                Eigen::VectorXd dqdT = Eigen::Map<Eigen::VectorXd>(&(last_dqdT[0]), 6, 1);

                double dFdT_n = dFdq_n2*dqdT;
                MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, it->nodeset->getSegByIx(s).getOrigin());
                it->DF[it->totalFree*row0 + epochVar.row0] = -dFdT_n/dvMax;
            }
        }
    }
}//==============================================

void DynamicsModel_bc4bp::multShoot_targetApse(MultShootData *it, Constraint con, int row0) const{
    std::vector<double> conData = con.getData();
    int Pix = static_cast<int>(conData[0]);    // index of primary
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    // int nodeIx = it->nodeset->getNodeIx(con.getID());

    // Get the node epoch either from the design vector or from the original set of nodes
    MSVarMap_Obj epochVar;
    double t0 = 0;
    if(it->bVarTime){
        epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        t0 = it->X[epochVar.row0];
    }else{
        t0 = it->nodeset->getEpoch(con.getID());
    }

    const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(it->sysData);
    double primPos[9], primVel[9], primAccel[9];

    getPrimaryPos(t0, bcSys, primPos);
    getPrimaryVel(t0, bcSys, primVel);
    getPrimaryAccel(t0, bcSys, primAccel);

    double dx = it->X[stateVar.row0+0] - primPos[3*Pix+0];
    double dy = it->X[stateVar.row0+1] - primPos[3*Pix+1];
    double dz = it->X[stateVar.row0+2] - primPos[3*Pix+2];
    double dvx = it->X[stateVar.row0+3] - primVel[3*Pix+0];
    double dvy = it->X[stateVar.row0+4] - primVel[3*Pix+1];
    double dvz = it->X[stateVar.row0+5] - primVel[3*Pix+2];

    it->FX[row0] = dx*dvx + dy*dvy + dz*dvz;

    it->DF[it->totalFree*row0 + stateVar.row0+0] = dvx;
    it->DF[it->totalFree*row0 + stateVar.row0+1] = dvy;
    it->DF[it->totalFree*row0 + stateVar.row0+2] = dvz;
    it->DF[it->totalFree*row0 + stateVar.row0+3] = dx;
    it->DF[it->totalFree*row0 + stateVar.row0+4] = dy;
    it->DF[it->totalFree*row0 + stateVar.row0+5] = dz;

    if(it->bVarTime){
        it->DF[it->totalFree*row0 + epochVar.row0] = -1*(dvx*primVel[3*Pix+0] + dvy*primVel[3*Pix+1] + dvz*primVel[3*Pix+2]);
        it->DF[it->totalFree*row0 + epochVar.row0] -= (dx*primAccel[3*Pix+0] + dy*primAccel[3*Pix+1] + dz*primAccel[3*Pix+2]);
    }
}//===================================================

/**
 *  \brief Compute partials and constraint values for nodes constrained with <tt>SP</tt>
 *
 *  This function computes three constraint values and three rows of partials for the Jacobian.
 *  Each row/function corresponds to one position state. The FX and DF matrices are updated
 *  in place by editing their values stored in <tt>it</tt>
 *
 *  \param it the MultShootData object holding the current data for the corrections process
 *  \param con the constraint being applied
 *  \param row0 the index of the first row for this constraint
 */
void DynamicsModel_bc4bp::multShoot_targetSP(MultShootData* it, Constraint con, int row0) const{
    // int nodeIx = it->nodeset->getNodeIx(con.getID());
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());

    // Get the node epoch either from the design vector or from the original set of nodes
    MSVarMap_Obj epochVar;
    double t0 = 0;
    if(it->bVarTime){
        epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        t0 = it->X[epochVar.row0];
    }else{
        t0 = it->nodeset->getEpoch(con.getID());
    }

    const SysData_bc4bp *bcSysData = static_cast<const SysData_bc4bp *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(t0, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double *X = &(it->X[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+stateVar.row0, 3, 1);   // Position vector

    // Create relative position vectors between s/c and primaries
    Eigen::Vector3d r_p1 = r - primPos.row(0).transpose();
    Eigen::Vector3d r_p2 = r - primPos.row(1).transpose();
    Eigen::Vector3d r_p3 = r - primPos.row(2).transpose();

    double d1 = r_p1.norm();
    double d2 = r_p2.norm();
    double d3 = r_p3.norm();

    double k = bcSysData->getK();
    double mu = bcSysData->getMu();
    double nu = bcSysData->getNu();

    // Evaluate three constraint function values 
    Eigen::Vector3d conEval;
    conEval.noalias() = -(1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2,3) - nu*r_p3/pow(d3, 3);
    
    // printf("Node %d: Saddle Point Constraint Data:", n);
    // printf("Epoch = %.16f\n", epoch);
    // printf("s = [%.16f, %.16f %.16f], mag = %.16f\n", r_p1(0), r_p1(1), r_p1(2), d1);
    // printf("e = [%.16f, %.16f %.16f], mag = %.16f\n", r_p2(0), r_p2(1), r_p2(2), d2);
    // printf("m = [%.16f, %.16f %.16f], mag = %.16f\n", r_p3(0), r_p3(1), r_p3(2), d3);
    // printf("mu = %.16f\n", mu);
    // printf("nu = %.16f\n", nu);
    // printf("FX = [%.16f, %.16f %.16f]\n", conEval(0), conEval(1), conEval(2));

    // Parials w.r.t. node position r
    double dFdq_data[9] = {0};
    dFdq_data[0] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(0),2)/pow(d1,5)) -
        (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
        3*pow(r_p3(0),2)/pow(d3,5));     //dxdx
    dFdq_data[1] = (1/k - mu)*3*r_p1(0)*r_p1(1)/pow(d1,5) + 
        (mu - nu)*3*r_p2(0)*r_p2(1)/pow(d2,5) +
        nu*3*r_p3(0)*r_p3(1)/pow(d3,5);   //dxdy
    dFdq_data[2] = (1/k - mu)*3*r_p1(0)*r_p1(2)/pow(d1,5) +
        (mu - nu)*3*r_p2(0)*r_p2(2)/pow(d2,5) +
        nu*3*r_p3(0)*r_p3(2)/pow(d3,5);   //dxdz
    dFdq_data[3] = dFdq_data[1];    // dydx = dxdy
    dFdq_data[4] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(1),2)/pow(d1,5)) -
        (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(1),2)/pow(d2,5)) - 
        nu*(1/pow(d3,3) - 3*pow(r_p3(1),2)/pow(d3,5));   //dydy
    dFdq_data[5] = (1/k - mu)*3*r_p1(1)*r_p1(2)/pow(d1,5) +
        (mu - nu)*3*r_p2(1)*r_p2(2)/pow(d2,5) +
        nu*3*r_p3(1)*r_p3(2)/pow(d3,5);   //dydz
    dFdq_data[6] = dFdq_data[2];    //dzdx = dxdz
    dFdq_data[7] = dFdq_data[5];    //dzdy = dydz
    dFdq_data[8] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(2),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
            3*pow(r_p3(2),2)/pow(d3,5)); //dzdz

    Matrix3Rd dFdq = Eigen::Map<Matrix3Rd>(dFdq_data, 3, 3);

    // Copy data into the correct vectors/matrices
    double* conEvalPtr = conEval.data();
    double* dFdq_ptr = dFdq.data();

    double *FX = &(it->FX[0]);
    double *DF = &(it->DF[0]);

    std::copy(conEvalPtr, conEvalPtr+3, FX+row0);
    std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + stateVar.row0);
    std::copy(dFdq_ptr+3, dFdq_ptr+6, DF + it->totalFree*(row0+1) + stateVar.row0);
    std::copy(dFdq_ptr+6, dFdq_ptr+9, DF + it->totalFree*(row0+2) + stateVar.row0);

    if(it->bVarTime){
        // Get primary velocities at the specified epoch time
        std::vector<double> primVelData = getPrimVel(t0, it->sysData);
        Matrix3Rd primVel = Eigen::Map<Matrix3Rd>(&(primVelData[0]), 3, 3);

        // Compute partials of state w.r.t. primary positions; dont' compute partials
        // for P1 because its velocity is zero in the rotating frame
        double dfdr2_data[18] = {0};   double dfdr3_data[18] = {0};

        dfdr2_data[9] = -1/pow(d2,3) + 3*pow(r_p2(0),2)/pow(d2,5);        //dxdx2
        dfdr2_data[10] = 3*r_p2(0)*r_p2(1)/pow(d2,5);                  //dxdy2
        dfdr2_data[11] = 3*r_p2(0)*r_p2(2)/pow(d2,5);                  //dxdz2
        dfdr2_data[13] = -1/pow(d2,3) + 3*pow(r_p2(1),2)/pow(d2,5);       //dydy2
        dfdr2_data[14] = 3*r_p2(1)*r_p2(2)/pow(d2,5);                  //dydz2
        dfdr2_data[17] = -1/pow(d2,3) + 3*pow(r_p2(2),2)/pow(d2,5);       //dzdz2

        dfdr2_data[12] = dfdr2_data[10];      // Fill in symmetric matrix
        dfdr2_data[15] = dfdr2_data[11];
        dfdr2_data[16] = dfdr2_data[14];

        dfdr3_data[9] = -1/pow(d3,3) + 3*pow(r_p3(0),2)/pow(d3,5);        //dxdx3
        dfdr3_data[10] = 3*r_p3(0)*r_p3(1)/pow(d3,5);                  //dxdy3
        dfdr3_data[11] = 3*r_p3(0)*r_p3(2)/pow(d3,5);                  //dxdz3
        dfdr3_data[13] = -1/pow(d3,3) + 3*pow(r_p3(1),2)/pow(d3,5);       //dydy3
        dfdr3_data[14] = 3*r_p3(1)*r_p3(2)/pow(d3,5);                  //dydz3
        dfdr3_data[17] = -1/pow(d3,3) + 3*pow(r_p3(2),2)/pow(d3,5);       //dzdz3

        dfdr3_data[12] = dfdr3_data[10];      // Fill in symmetric matrix
        dfdr3_data[15] = dfdr3_data[11];
        dfdr3_data[16] = dfdr3_data[14];

        MatrixXRd dFdr2 = Eigen::Map<MatrixXRd>(dfdr2_data, 6, 3);
        MatrixXRd dFdr3 = Eigen::Map<MatrixXRd>(dfdr3_data, 6, 3);

        // scale matrices by constants
        dFdr2 *= -1*(mu - nu);
        dFdr3 *= -1*nu;

        // Compute partials of constraint function w.r.t. epoch time
        Eigen::VectorXd dFdT;
        dFdT.noalias() = dFdr2*(primVel.row(1).transpose()) + dFdr3*(primVel.row(2).transpose());

        DF[it->totalFree*row0 + epochVar.row0] = dFdT(3);
        DF[it->totalFree*(row0+1) + epochVar.row0] = dFdT(4);
        DF[it->totalFree*(row0+2) + epochVar.row0] = dFdT(5);
    }
}// End of SP Targeting ==============================

/**
 *  \brief Compute partials and constraint values for nodes constrained with <tt>SP_RANGE</tt>
 *  \details This function computes one constraint value and one row of partials for the Jacobian
 *  because the total acceleration magnitude is targeted rather than individual acceleration 
 *  vector components. The FX and DF matrices are update in place by editing their values 
 *  stored in <tt>it</tt>
 * 
 *  \param it the MultShootData object holding the current data for the corrections process
 *  \param con the constraint being applied
 *  \param c the index of the constraint within the constraint vector
 */
void DynamicsModel_bc4bp::multShoot_targetSP_mag(MultShootData* it, Constraint con, int c) const{

    // int row0 = it->conRows[c];
    // int n = con.getID();
    // double Amax = con.getData()[0];
    // int epochCol = it->bEqualArcTime ? 6*it->numNodes+1+n : 7*it->numNodes-1+n;
    // double t0 = it->bVarTime ? it->X[epochCol]/it->freeVarScale[3] : it->nodeset->getNodeByIx(n).getEpoch();

    // const SysData_bc4bp *bcSysData = static_cast<const SysData_bc4bp *> (it->sysData);
    // double sr = it->freeVarScale[0];

    // std::vector<double> primPosData = getPrimPos(t0, it->sysData);

    // // Get primary positions at the specified epoch time
    // Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    // double *X = &(it->X[0]);
    // Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+6*n, 3, 1);   // Position vector

    // // Create relative position vectors between s/c and primaries
    // Eigen::Vector3d r_p1 = r/sr - primPos.row(0).transpose();
    // Eigen::Vector3d r_p2 = r/sr - primPos.row(1).transpose();
    // Eigen::Vector3d r_p3 = r/sr - primPos.row(2).transpose();

    // double d1 = r_p1.norm();
    // double d2 = r_p2.norm();
    // double d3 = r_p3.norm();

    // double k = bcSysData->getK();
    // double mu = bcSysData->getMu();
    // double nu = bcSysData->getNu();

    // // Compute the acceleration vector at this point
    // Eigen::Vector3d A;
    // A.noalias() = -(1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2,3) - nu*r_p3/pow(d3, 3);

    // // Parials w.r.t. node position r
    // double dFdq_data[9] = {0};
    // dFdq_data[0] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(0),2)/pow(d1,5)) -
    //     (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
    //     3*pow(r_p3(0),2)/pow(d3,5));     //dxdx
    // dFdq_data[1] = (1/k - mu)*3*r_p1(0)*r_p1(1)/pow(d1,5) + 
    //     (mu - nu)*3*r_p2(0)*r_p2(1)/pow(d2,5) +
    //     nu*3*r_p3(0)*r_p3(1)/pow(d3,5);   //dxdy
    // dFdq_data[2] = (1/k - mu)*3*r_p1(0)*r_p1(2)/pow(d1,5) +
    //     (mu - nu)*3*r_p2(0)*r_p2(2)/pow(d2,5) +
    //     nu*3*r_p3(0)*r_p3(2)/pow(d3,5);   //dxdz
    // dFdq_data[3] = dFdq_data[1];    // dydx = dxdy
    // dFdq_data[4] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(1),2)/pow(d1,5)) -
    //     (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(1),2)/pow(d2,5)) - 
    //     nu*(1/pow(d3,3) - 3*pow(r_p3(1),2)/pow(d3,5));   //dydy
    // dFdq_data[5] = (1/k - mu)*3*r_p1(1)*r_p1(2)/pow(d1,5) +
    //     (mu - nu)*3*r_p2(1)*r_p2(2)/pow(d2,5) +
    //     nu*3*r_p3(1)*r_p3(2)/pow(d3,5);   //dydz
    // dFdq_data[6] = dFdq_data[2];    //dzdx = dxdz
    // dFdq_data[7] = dFdq_data[5];    //dzdy = dydz
    // dFdq_data[8] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(2),2)/pow(d1,5)) -
    //     (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
    //     3*pow(r_p3(2),2)/pow(d3,5)); //dzdz

    // Matrix3Rd dAdq = Eigen::Map<Matrix3Rd>(dFdq_data, 3, 3);
    // Eigen::Vector3d dFdq;
    // dFdq.noalias() = 2*dAdq*A/(sr*Amax*Amax);
    // // dFdq.noalias() = 2*dAdq*A/sr;

    // // Copy data into the correct vectors/matrices
    // double* dFdq_ptr = dFdq.data();

    // double *FX = &(it->FX[0]);
    // double *DF = &(it->DF[0]);

    // FX[row0] = A.squaredNorm()/(Amax*Amax) - 1;     // Found this one converges MUCH better
    // // FX[row0] = A.squaredNorm() - Amax*Amax;
    
    // std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + 6*n);

    // if(it->bVarTime){
    //     // Get primary velocities at the specified epoch time
    //     std::vector<double> primVelData = getPrimVel(t0, it->sysData);
    //     Matrix3Rd primVel = Eigen::Map<Matrix3Rd>(&(primVelData[0]), 3, 3);

    //     // Compute partials of state w.r.t. primary positions; dont' compute partials
    //     // for P1 because its velocity is zero in the rotating frame
    //     double dfdr2_data[18] = {0};   double dfdr3_data[18] = {0};

    //     dfdr2_data[9] = -1/pow(d2,3) + 3*pow(r_p2(0),2)/pow(d2,5);        //dxdx2
    //     dfdr2_data[10] = 3*r_p2(0)*r_p2(1)/pow(d2,5);                  //dxdy2
    //     dfdr2_data[11] = 3*r_p2(0)*r_p2(2)/pow(d2,5);                  //dxdz2
    //     dfdr2_data[13] = -1/pow(d2,3) + 3*pow(r_p2(1),2)/pow(d2,5);       //dydy2
    //     dfdr2_data[14] = 3*r_p2(1)*r_p2(2)/pow(d2,5);                  //dydz2
    //     dfdr2_data[17] = -1/pow(d2,3) + 3*pow(r_p2(2),2)/pow(d2,5);       //dzdz2

    //     dfdr2_data[12] = dfdr2_data[10];      // Fill in symmetric matrix
    //     dfdr2_data[15] = dfdr2_data[11];
    //     dfdr2_data[16] = dfdr2_data[14];

    //     dfdr3_data[9] = -1/pow(d3,3) + 3*pow(r_p3(0),2)/pow(d3,5);        //dxdx3
    //     dfdr3_data[10] = 3*r_p3(0)*r_p3(1)/pow(d3,5);                  //dxdy3
    //     dfdr3_data[11] = 3*r_p3(0)*r_p3(2)/pow(d3,5);                  //dxdz3
    //     dfdr3_data[13] = -1/pow(d3,3) + 3*pow(r_p3(1),2)/pow(d3,5);       //dydy3
    //     dfdr3_data[14] = 3*r_p3(1)*r_p3(2)/pow(d3,5);                  //dydz3
    //     dfdr3_data[17] = -1/pow(d3,3) + 3*pow(r_p3(2),2)/pow(d3,5);       //dzdz3

    //     dfdr3_data[12] = dfdr3_data[10];      // Fill in symmetric matrix
    //     dfdr3_data[15] = dfdr3_data[11];
    //     dfdr3_data[16] = dfdr3_data[14];

    //     Matrix3Rd dFdr2 = Eigen::Map<Matrix3Rd>(dfdr2_data+9, 3, 3);
    //     Matrix3Rd dFdr3 = Eigen::Map<Matrix3Rd>(dfdr3_data+9, 3, 3);

    //     // scale matrices by constants
    //     dFdr2 *= -1*(mu - nu)/sr;
    //     dFdr3 *= -1*nu/sr;

    //     Eigen::VectorXd dAdT;
    //     dAdT.noalias() = A.transpose()*dFdr2*primVel.row(1).transpose() + A.transpose()*dFdr3*primVel.row(2).transpose();
    //     dAdT *= 2*sr/(Amax*Amax*it->freeVarScale[3]);
    //     // dAdT *= 2*sr/it->freeVarScale[3];

    //     int epochCol = it->bEqualArcTime ? 6*it->numNodes+1+n : 7*it->numNodes-1+n;
    //     DF[it->totalFree*row0 + epochCol] = dAdT(0);
    // }

    // // figure out which of the slack variables correspond to this constraint
    // std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
    //     it->slackAssignCon.end(), c);

    // // which column of the DF matrix the slack variable is in
    // int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());

    // // Add squared slack variable from constraint
    // it->FX[row0] += it->X[slackCol]*it->X[slackCol];

    // // Partial with respect to slack variable
    // it->DF[it->totalFree*row0 + slackCol] = 2*it->X[slackCol];
    
    // // printf("Node %d %s Con: Slack Var = %.4e\n", n, con.getTypeStr(), it->X[slackCol]);

    int row0 = it->conRows[c];
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    // int nodeIx = it->nodeset->getNodeIx(con.getID());
    double Amax = con.getData()[0];

    // Get the node epoch either from the design vector or from the original set of nodes
    MSVarMap_Obj epochVar;
    double epoch = 0;
    if(it->bVarTime){
        epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        epoch = it->X[epochVar.row0];
    }else{
        epoch = it->nodeset->getEpoch(con.getID());
    }
    const SysData_bc4bp *bcSysData = static_cast<const SysData_bc4bp *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(epoch, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double *X = &(it->X[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+stateVar.row0, 3, 1);   // Position vector

    // Create relative position vectors between s/c and primaries
    Eigen::Vector3d r_p1 = r - primPos.row(0).transpose();
    Eigen::Vector3d r_p2 = r - primPos.row(1).transpose();
    Eigen::Vector3d r_p3 = r - primPos.row(2).transpose();

    double d1 = r_p1.norm();
    double d2 = r_p2.norm();
    double d3 = r_p3.norm();

    double k = bcSysData->getK();
    double mu = bcSysData->getMu();
    double nu = bcSysData->getNu();

    // Compute the acceleration vector at this point
    Eigen::Vector3d A;
    A.noalias() = -(1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2,3) - nu*r_p3/pow(d3, 3);

    // Parials w.r.t. node position r
    double dFdq_data[9] = {0};
    dFdq_data[0] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(0),2)/pow(d1,5)) -
        (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
        3*pow(r_p3(0),2)/pow(d3,5));     //dxdx
    dFdq_data[1] = (1/k - mu)*3*r_p1(0)*r_p1(1)/pow(d1,5) + 
        (mu - nu)*3*r_p2(0)*r_p2(1)/pow(d2,5) +
        nu*3*r_p3(0)*r_p3(1)/pow(d3,5);   //dxdy
    dFdq_data[2] = (1/k - mu)*3*r_p1(0)*r_p1(2)/pow(d1,5) +
        (mu - nu)*3*r_p2(0)*r_p2(2)/pow(d2,5) +
        nu*3*r_p3(0)*r_p3(2)/pow(d3,5);   //dxdz
    dFdq_data[3] = dFdq_data[1];    // dydx = dxdy
    dFdq_data[4] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(1),2)/pow(d1,5)) -
        (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(1),2)/pow(d2,5)) - 
        nu*(1/pow(d3,3) - 3*pow(r_p3(1),2)/pow(d3,5));   //dydy
    dFdq_data[5] = (1/k - mu)*3*r_p1(1)*r_p1(2)/pow(d1,5) +
        (mu - nu)*3*r_p2(1)*r_p2(2)/pow(d2,5) +
        nu*3*r_p3(1)*r_p3(2)/pow(d3,5);   //dydz
    dFdq_data[6] = dFdq_data[2];    //dzdx = dxdz
    dFdq_data[7] = dFdq_data[5];    //dzdy = dydz
    dFdq_data[8] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1(2),2)/pow(d1,5)) -
        (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
        3*pow(r_p3(2),2)/pow(d3,5)); //dzdz

    Matrix3Rd dAdq = Eigen::Map<Matrix3Rd>(dFdq_data, 3, 3);
    Eigen::Vector3d dFdq;
    dFdq.noalias() = 2*dAdq*A/(Amax*Amax);
    // dFdq.noalias() = 2*dAdq*A;

    // Get primary velocities at the specified epoch time
    std::vector<double> primVelData = getPrimVel(epoch, it->sysData);
    Matrix3Rd primVel = Eigen::Map<Matrix3Rd>(&(primVelData[0]), 3, 3);

    // Compute partials of state w.r.t. primary positions; dont' compute partials
    // for P1 because its velocity is zero in the rotating frame
    double dfdr2_data[18] = {0};   double dfdr3_data[18] = {0};

    dfdr2_data[9] = -1/pow(d2,3) + 3*pow(r_p2(0),2)/pow(d2,5);        //dxdx2
    dfdr2_data[10] = 3*r_p2(0)*r_p2(1)/pow(d2,5);                  //dxdy2
    dfdr2_data[11] = 3*r_p2(0)*r_p2(2)/pow(d2,5);                  //dxdz2
    dfdr2_data[13] = -1/pow(d2,3) + 3*pow(r_p2(1),2)/pow(d2,5);       //dydy2
    dfdr2_data[14] = 3*r_p2(1)*r_p2(2)/pow(d2,5);                  //dydz2
    dfdr2_data[17] = -1/pow(d2,3) + 3*pow(r_p2(2),2)/pow(d2,5);       //dzdz2

    dfdr2_data[12] = dfdr2_data[10];      // Fill in symmetric matrix
    dfdr2_data[15] = dfdr2_data[11];
    dfdr2_data[16] = dfdr2_data[14];

    dfdr3_data[9] = -1/pow(d3,3) + 3*pow(r_p3(0),2)/pow(d3,5);        //dxdx3
    dfdr3_data[10] = 3*r_p3(0)*r_p3(1)/pow(d3,5);                  //dxdy3
    dfdr3_data[11] = 3*r_p3(0)*r_p3(2)/pow(d3,5);                  //dxdz3
    dfdr3_data[13] = -1/pow(d3,3) + 3*pow(r_p3(1),2)/pow(d3,5);       //dydy3
    dfdr3_data[14] = 3*r_p3(1)*r_p3(2)/pow(d3,5);                  //dydz3
    dfdr3_data[17] = -1/pow(d3,3) + 3*pow(r_p3(2),2)/pow(d3,5);       //dzdz3

    dfdr3_data[12] = dfdr3_data[10];      // Fill in symmetric matrix
    dfdr3_data[15] = dfdr3_data[11];
    dfdr3_data[16] = dfdr3_data[14];

    Matrix3Rd dFdr2 = Eigen::Map<Matrix3Rd>(dfdr2_data+9, 3, 3);
    Matrix3Rd dFdr3 = Eigen::Map<Matrix3Rd>(dfdr3_data+9, 3, 3);

    // scale matrices by constants
    dFdr2 *= -1*(mu - nu);
    dFdr3 *= -1*nu;

    Eigen::VectorXd dAdT;
    dAdT.noalias() = A.transpose()*dFdr2*primVel.row(1).transpose() + A.transpose()*dFdr3*primVel.row(2).transpose();
    dAdT *= 2/(Amax*Amax);
    // dAdT *= 2;

    // Copy data into the correct vectors/matrices
    double* dFdq_ptr = dFdq.data();
    double* dFdT_ptr = dAdT.data();

    double *FX = &(it->FX[0]);
    double *DF = &(it->DF[0]);

    FX[row0] = A.squaredNorm()/(Amax*Amax) - 1;
    // FX[row0] = A.squaredNorm() - Amax*Amax;
    
    std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + stateVar.row0);

    if(it->bVarTime){
        std::copy(dFdT_ptr, dFdT_ptr+1, DF + it->totalFree*row0 + epochVar.row0);
    }

    // figure out which of the slack variables correspond to this constraint
    std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
        it->slackAssignCon.end(), c);

    // which column of the DF matrix the slack variable is in
    int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());

    // Add squared slack variable from constraint
    it->FX[row0] += it->X[slackCol]*it->X[slackCol];

    // Partial with respect to slack variable
    it->DF[it->totalFree*row0 + slackCol] = 2*it->X[slackCol];
}// End of SP Targeting (Magnitude) ==============================

/**
 *  \brief Compute an initial value for the slack variable for an SP_RANGE inequality constraint 
 *  \details If the constraint is satisified, the slack variable must be set such that the 
 *  constraint function evaluates to zero; that value is computed here. If the constraint is not
 *  satisfied, a small value is returned.
 * 
 *  \param it the MultShootData object holding the current data for the corrections process
 *  \param con the constraint being applied
 * 
 *  \return the initial value for the slack variable associated with an SP_RANGE constraint
 */
double DynamicsModel_bc4bp::multShoot_targetSPMag_compSlackVar(const MultShootData *it, Constraint con) const{
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    // int nodeIx = it->nodeset->getNodeIx(con.getID());
    double Amax = con.getData()[0];

    // Get the node epoch either from the design vector or from the original set of nodes
    double epoch = 0;
    if(it->bVarTime){
        MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        epoch = it->X[epochVar.row0];
    }else{
        epoch = it->nodeset->getEpoch(con.getID());
    }
    const SysData_bc4bp *bcSysData = static_cast<const SysData_bc4bp *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(epoch, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double rData[3];
    std::copy(&(it->X[0])+stateVar.row0, &(it->X[0])+stateVar.row0+3, rData);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(rData, 3, 1);   // Position vector

    // Create relative position vectors between s/c and primaries
    Eigen::Vector3d r_p1 = r - primPos.row(0).transpose();
    Eigen::Vector3d r_p2 = r - primPos.row(1).transpose();
    Eigen::Vector3d r_p3 = r - primPos.row(2).transpose();

    double d1 = r_p1.norm();
    double d2 = r_p2.norm();
    double d3 = r_p3.norm();

    double k = bcSysData->getK();
    double mu = bcSysData->getMu();
    double nu = bcSysData->getNu();

    // Compute the acceleration vector at this point
    Eigen::Vector3d A;
    A.noalias() = -(1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2,3) - nu*r_p3/pow(d3, 3);
    double diff = 1 - A.squaredNorm()/(Amax*Amax);
    // double diff = Amax*Amax - A.squaredNorm();

    /*  If diff is positive, then the constraint
     *  is satisfied, so compute the value of the slack variable that 
     *  sets the constraint function equal to zero. Otherwise, choose 
     *  a small value of the slack variable but don't set it to zero as 
     *  that will make the partials zero and will prevent the mulitple
     *  shooting algorithm from updating the slack variable
     */
    return diff > 0 ? sqrt(diff) : 1e-4;
}//======================================================

/**
 *  \brief Compute an initial value for the slack variable for an SP_DIST inequality constraint 
 *  \details If the constraint is satisified, the slack variable must be set such that the 
 *  constraint function evaluates to zero; that value is computed here. If the constraint is not
 *  satisfied, a small value is returned.
 * 
 *  \param it the MultShootData object holding the current data for the corrections process
 *  \param con the constraint being applied
 * 
 *  \return the initial value for the slack variable associated with an SP_DIST constraint
 */
double DynamicsModel_bc4bp::multShoot_targetSP_maxDist_compSlackVar(const MultShootData *it, Constraint con) const{
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    double T = 0;
    if(it->bVarTime){
        MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        T = it->X[epochVar.row0];
    }else{
        T = it->nodeset->getEpoch(con.getID());
    }

    std::vector<double> coeff = con.getData();

    // Compute SP position from polynomial approximation
    Eigen::Vector3d spPos;
    spPos(0) = T*T*coeff[1] + T*coeff[2] + coeff[3];
    spPos(1) = T*T*coeff[4] + T*coeff[5] + coeff[6];
    spPos(2) = T*T*coeff[7] + T*coeff[8] + coeff[9];

    // double *X = &(it->X[0]);
    double rData[3];
    std::copy(&(it->X[0])+stateVar.row0, &(it->X[0])+stateVar.row0+3, rData);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(rData, 3, 1);   // Position vector

    Eigen::Vector3d dist = r - spPos;
    // double diff = coeff[0] - dist.norm();
    double diff = coeff[0]*coeff[0] - dist.squaredNorm();
    
    return diff > 0 ? sqrt(diff) : 1e-4;
}//===================================================

/**
 *  \brief Compute partials and constraint values for node constrained with <tt>SP_DIST</tt>
 *  and <tt>SP_MAX_DIST</tt>. 
 *  \details One constraint value and one row of partials are computed. This constraint uses
 *  2nd-order polynomials to approximate the saddle point's location as a function of epoch 
 *  and then targets a node to be at or within a set radius of the saddle point.
 * 
 *  \param it the MultShootData object holding the current data for the corrections process
 *  \param con the constraint being applied
 *  \param c the index of the constraint within the storage vector
 */
void DynamicsModel_bc4bp::multShoot_targetSP_dist(MultShootData *it, Constraint con, int c) const{
    int row0 = it->conRows[c];
    MSVarMap_Obj stateVar = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    MSVarMap_Obj epochVar;
    double T = 0;
    if(it->bVarTime){
        epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, con.getID());
        T = it->X[epochVar.row0];
    }else{
        T = it->nodeset->getEpoch(con.getID());
    }

    // const SysData_bc4bp *bcSysData = static_cast<const SysData_bc4bp *> (it->sysData);
    std::vector<double> coeff = con.getData();

    // Compute SP position from polynomial approximation
    Eigen::Vector3d spPos;
    spPos(0) = T*T*coeff[1] + T*coeff[2] + coeff[3];
    spPos(1) = T*T*coeff[4] + T*coeff[5] + coeff[6];
    spPos(2) = T*T*coeff[7] + T*coeff[8] + coeff[9];

    double *X = &(it->X[0]);
    double *DF = &(it->DF[0]);
    double *FX = &(it->FX[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+stateVar.row0, 3, 1);   // Position vector

    Eigen::Vector3d dist = r - spPos;

    FX[row0] = dist.squaredNorm() - coeff[0]*coeff[0];    // This one is much more robust, so choose this one
    // double d = dist.norm();          
    // FX[row0] = d - coeff[0]*sr;    // This one converges much more quickly but is much less robust

    // Compute partials w.r.t. node states
    DF[it->totalFree*row0 + stateVar.row0+0] = 2*dist(0);
    DF[it->totalFree*row0 + stateVar.row0+1] = 2*dist(1);
    DF[it->totalFree*row0 + stateVar.row0+2] = 2*dist(2);
    // DF[it->totalFree*row0 + stateVar.row0+0] = dist(0)/d;     // dFdx
    // DF[it->totalFree*row0 + stateVar.row0+1] = dist(1)/d;     // dFdy
    // DF[it->totalFree*row0 + stateVar.row0+2] = dist(2)/d;     // dFdz

    // printf("Dist from Approx SP = [%.4f, %.4f, %.4f]\n", dist(0), dist(1), dist(2));
    // printf("  r = [%.4f, %.4f, %.4f]\n", r(0), r(1), r(2));
    // printf("  sp_pos = [%.4f, %.4f, %.4f]\n", spPos(0), spPos(1), spPos(2));

    // Compute partials w.r.t. epoch Time
    if(it->bVarTime){
        DF[it->totalFree*row0 + epochVar.row0] = -2*dist(0)*(2*coeff[1]*T + coeff[2]) -
            2*dist(1)*(2*coeff[4]*T + coeff[5]) - 2*dist(2)*(2*coeff[7]*T + coeff[8]);
        // DF[it->totalFree*row0 + epochVar.row0] = -1/d*dist(0)*(2*coeff[1]*T + coeff[2]) -
        //     1/d*dist(1)*(2*coeff[4]*T + coeff[5]) - 1/d*dist(2)*(2*coeff[7]*T + coeff[8]);
    }

    if(con.getType() == Constraint_tp::SP_MAX_DIST){
        // figure out which of the slack variables correspond to this constraint
        std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
            it->slackAssignCon.end(), c);

        // which column of the DF matrix the slack variable is in
        int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());

        // Add squared slack variable from constraint
        it->FX[row0] += it->X[slackCol]*it->X[slackCol];

        // Partial with respect to slack variable
        it->DF[it->totalFree*row0 + slackCol] = 2*it->X[slackCol];
        
        // printf("Node %d distance to SP = %.4f km\n", n, dist.norm()*bcSysData->getCharL());
        // if(FX[row0] < 0)     // I've tested this method in test_multShootCons and found no difference in convergence speed
        //     FX[row0]= 0;
    }

    // printf("SP_DIST constraint: dist = %.4e = %.4f km, ||F|| = %.4e\n", d, d*bcSysData->getCharL(), it->FX[row0]);

}//===================================================

/**
 *  \brief Take the final, corrected free variable vector <tt>X</tt> and create an output 
 *  nodeset
 *
 *  If <tt>findEvent</tt> is set to true, the
 *  output nodeset will contain extra information for the simulation engine to use. Rather than
 *  returning only the position and velocity states, the output nodeset will contain the STM 
 *  and dqdT values for the final node; this information will be appended to the extraParameter
 *  vector in the final node.
 *
 *  \param it an iteration data object containing all info from the corrections process
 *  \param nodes_in a pointer to the original, uncorrected nodeset
 *  \param findEvent whether or not this correction process is locating an event
 *  \param nodes_out pointer to the nodeset object that will contain the output of the
 *  shooting process
 *  
 *  \return a pointer to a nodeset containing the corrected nodes
 */
void DynamicsModel_bc4bp::multShoot_createOutput(const MultShootData *it, const Nodeset *nodes_in, bool findEvent, Nodeset *nodes_out) const{

    Nodeset_bc4bp *nodeset_out = static_cast<Nodeset_bc4bp *>(nodes_out);

    std::vector<int> newNodeIDs;
    for(int n = 0; n < it->numNodes; n++){
        MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getNodeByIx(n).getID());
        // double state[6];
        std::vector<double> state(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0 + coreStates);
        // std::copy(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0+6, state);

        Node node(state, 0);
        node.setConstraints(nodes_in->getNodeByIx(n).getConstraints());

        double T = 0;
        if(it->bVarTime){
            MSVarMap_Obj epochVar = it->getVarMap_obj(MSVar_tp::EPOCH, state_var.key.id);
            T = it->X[epochVar.row0];
        }else{
            T = it->nodeset->getEpoch(state_var.key.id);
        }
        node.setEpoch(T);

        if(n+1 == it->numNodes){
            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
                state for the last node. We do this by appending the STM elements and dqdT elements to the
                end of the node array. This output nodeset should have two "nodes": the first 6 elements
                are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
                information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                Traj lastSeg = it->propSegs.back();
                MatrixXRd lastSTM = lastSeg.getSTMByIx(-1);
                std::vector<double> stm_vec(lastSTM.data(), lastSTM.data() + lastSTM.rows()*lastSTM.cols());
                node.setExtraParamVec(PARAMKEY_STM, stm_vec);

                // Also append the last dqdT vector
                std::vector<double> dqdT = lastSeg.getExtraParamVecByIx(-1,PARAMKEY_STATE_EPOCH_DERIV);
                node.setExtraParamVec(PARAMKEY_STATE_EPOCH_DERIV, dqdT);
            }
        }

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(nodeset_out->addNode(node));
    }

    double tof;
    int newOrigID, newTermID;
    for(int s = 0; s < it->nodeset->getNumSegs(); s++){
        Segment seg = it->nodeset->getSegByIx(s);

        if(it->bVarTime){
            MSVarMap_Obj tofVar = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
                it->bEqualArcTime ? Linkable::INVALID_ID : seg.getID());
            // Get data
            tof = it->bEqualArcTime ? it->X[tofVar.row0]/(it->nodeset->getNumSegs()) : it->X[tofVar.row0];
        }else{
            tof = seg.getTOF();
        }

        newOrigID = newNodeIDs[it->nodeset->getNodeIx(seg.getOrigin())];
        int termID = seg.getTerminus();
        newTermID = termID == Linkable::INVALID_ID ? termID : newNodeIDs[it->nodeset->getNodeIx(termID)];
        
        Segment newSeg(newOrigID, newTermID, tof);
        newSeg.setConstraints(seg.getConstraints());
        newSeg.setVelCon(seg.getVelCon());
        newSeg.setSTM(it->propSegs[s].getSTMByIx(-1));
        newSeg.setCtrlLaw(seg.getCtrlLaw());
        nodeset_out->addSeg(newSeg);
    }

    std::vector<Constraint> arcCons = nodes_in->getArcConstraints();
    for(unsigned int i = 0; i < arcCons.size(); i++){
        nodeset_out->addConstraint(arcCons[i]);
    }
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *   \brief Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   \param t epoch at integration step
 *   \param s the 48-d state vector
 *   \param sdot the 48-d state derivative vector
 *   \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_bc4bp::fullEOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    // SysData_bc4bp *sysData = static_cast<SysData_bc4bp *>(params);
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_bc4bp *sysData = static_cast<const SysData_bc4bp *>(paramStruct->sysData);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    getPrimaryPos(t, sysData, primPosData);

    // Create relative position vectors between s/c and primaries
    double r_p1[] = {s[0] - primPosData[0], s[1] - primPosData[1], s[2] - primPosData[2]};
    double r_p2[] = {s[0] - primPosData[3], s[1] - primPosData[4], s[2] - primPosData[5]};
    double r_p3[] = {s[0] - primPosData[6], s[1] - primPosData[7], s[2] - primPosData[8]};
    double d1 = std::sqrt(r_p1[0]*r_p1[0] + r_p1[1]*r_p1[1] + r_p1[2]*r_p1[2]);
    double d2 = std::sqrt(r_p2[0]*r_p2[0] + r_p2[1]*r_p2[1] + r_p2[2]*r_p2[2]);
    double d3 = std::sqrt(r_p3[0]*r_p3[0] + r_p3[1]*r_p3[1] + r_p3[2]*r_p3[2]);
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Velocity
    std::copy(s+3, s+6, sdot);

    // Compute acceleration
    sdot[3] = 2*k*s[4] + k*k*(s[0] + 1/k - mu) - (1/k - mu)*(s[0] - primPosData[0])/pow(d1,3) -
        (mu-nu)*(s[0] - primPosData[3])/pow(d2,3) - nu*(s[0] - primPosData[6])/pow(d3,3);
    sdot[4] = -2*k*s[3] + k*k*s[1] - (1/k - mu)*(s[1] - primPosData[1])/pow(d1,3) -
        (mu-nu)*(s[1] - primPosData[4])/pow(d2,3) - nu*(s[1] - primPosData[7])/pow(d3,3);
    sdot[5] = -1*(1/k - mu)*(s[2] - primPosData[2])/pow(d1,3) -
        (mu-nu)*(s[2] - primPosData[5])/pow(d2,3) - nu*(s[2] - primPosData[8])/pow(d3,3);

    // Compute psuedo-potential
    double dxdx = k*k - (1/k - mu)*(1/pow(d1,3) - 3*r_p1[0]*r_p1[0]/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*r_p2[0]*r_p2[0]/pow(d2,5)) - nu*(1/pow(d3,3) -
                3*r_p3[0]*r_p3[0]/pow(d3,5));
    double dxdy = (1/k - mu)*3*r_p1[0]*r_p1[1]/pow(d1,5) +
            (mu - nu)*3*r_p2[0]*r_p2[1]/pow(d2,5) +
            nu*3*r_p3[0]*r_p3[1]/pow(d3,5);
    double dxdz = (1/k - mu)*3*r_p1[0]*r_p1[2]/pow(d1,5) +
            (mu - nu)*3*r_p2[0]*r_p2[2]/pow(d2,5) +
            nu*3*r_p3[0]*r_p3[2]/pow(d3,5);
    double dydy = k*k - (1/k - mu)*(1/pow(d1,3) - 3*r_p1[1]*r_p1[1]/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*r_p2[1]*r_p2[1]/pow(d2,5)) - nu*(1/pow(d3,3) -
            3*r_p3[1]*r_p3[1]/pow(d3,5));
    double dydz = (1/k - mu)*3*r_p1[1]*r_p1[2]/pow(d1,5) +
            (mu - nu)*3*r_p2[1]*r_p2[2]/pow(d2,5) +
            nu*3*r_p3[1]*r_p3[2]/pow(d3,5);
    double dzdz = -(1/k - mu)*(1/pow(d1,3) - 3*r_p1[2]*r_p1[2]/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*r_p2[2]*r_p2[2]/pow(d2,5)) - nu*(1/pow(d3,3) -
            3*r_p3[2]*r_p3[2]/pow(d3,5));

    /*  Compute the STM Derivative 
     *  PhiDot = A * Phi
     *  s[6] through s[42] represent the STM, Phi, in row-major order 
     *  sdot [6] through [42] is thus the derivative of the STM
     */
    std::copy(s+24, s+42, sdot+6); // First three rows are the last three rows of Phi
    for(int i = 0; i < 6; i++){
        sdot[24+i] = dxdx*s[6+i] + dxdy*s[12+i] + dxdz*s[18+i] + 2*k*s[30+i];
        sdot[30+i] = dxdy*s[6+i] + dydy*s[12+i] + dydz*s[18+i] - 2*k*s[24+i];
        sdot[36+i] = dxdz*s[6+i] + dydz*s[12+i] + dzdz*s[18+i];
    }   // Last three rows are a combo of A and Phi
    

    // Compute partials of state w.r.t. primary positions; dont' compute partials
    // for P1 because its velocity is zero in the rotating frame
    double dfdr2[18] = {0};   double dfdr3[18] = {0};

    dfdr2[9] = -1*(mu-nu) * (-1/pow(d2,3) + 3*r_p2[0]*r_p2[0]/pow(d2,5));       //dxdx2
    dfdr2[10] = -1*(mu-nu) * (3*r_p2[0]*r_p2[1]/pow(d2,5));                     //dxdy2
    dfdr2[11] = -1*(mu-nu) * (3*r_p2[0]*r_p2[2]/pow(d2,5));                     //dxdz2
    dfdr2[13] = -1*(mu-nu) * (-1/pow(d2,3) + 3*r_p2[1]*r_p2[1]/pow(d2,5));      //dydy2
    dfdr2[14] = -1*(mu-nu) * (3*r_p2[1]*r_p2[2]/pow(d2,5));                     //dydz2
    dfdr2[17] = -1*(mu-nu) * (-1/pow(d2,3) + 3*r_p2[2]*r_p2[2]/pow(d2,5));      //dzdz2

    dfdr2[12] = dfdr2[10];      // Fill in symmetric matrix
    dfdr2[15] = dfdr2[11];
    dfdr2[16] = dfdr2[14];

    dfdr3[9] = -nu * (-1/pow(d3,3) + 3*r_p3[0]*r_p3[0]/pow(d3,5));              //dxdx3
    dfdr3[10] = -nu * (3*r_p3[0]*r_p3[1]/pow(d3,5));                            //dxdy3
    dfdr3[11] = -nu * (3*r_p3[0]*r_p3[2]/pow(d3,5));                            //dxdz3
    dfdr3[13] = -nu * (-1/pow(d3,3) + 3*r_p3[1]*r_p3[1]/pow(d3,5));             //dydy3
    dfdr3[14] = -nu * (3*r_p3[1]*r_p3[2]/pow(d3,5));                            //dydz3
    dfdr3[17] = -nu * (-1/pow(d3,3) + 3*r_p3[2]*r_p3[2]/pow(d3,5));             //dzdz3

    dfdr3[12] = dfdr3[10];      // Fill in symmetric matrix
    dfdr3[15] = dfdr3[11];
    dfdr3[16] = dfdr3[14];
    
    // Compute time derivative of dqdT
    std::copy(s+45, s+48, sdot+42);     // First three rows are the same as the last three rows of dqdT
    double primVelData[9] = {0};
    getPrimaryVel(t, sysData, primVelData);

    // For the final three rows, first compute the product of A*dqdT
    sdot[45] = dxdx*s[42] + dxdy*s[43] + dxdz*s[44] + 2*k*s[46];
    sdot[46] = dxdy*s[42] + dydy*s[43] + dydz*s[44] - 2*k*s[45];
    sdot[47] = dxdz*s[42] + dydz*s[43] + dzdz*s[44];
    
    // Next, add product of df/dr_2 and v_2
    sdot[45] += dfdr2[9]*primVelData[3] + dfdr2[10]*primVelData[4] + dfdr2[11]*primVelData[5];
    sdot[46] += dfdr2[12]*primVelData[3] + dfdr2[13]*primVelData[4] + dfdr2[14]*primVelData[5];
    sdot[47] += dfdr2[15]*primVelData[3] + dfdr2[16]*primVelData[4] + dfdr2[17]*primVelData[5];
    
    // Finally, add product of df/dr_3 and v_3
    sdot[45] += dfdr3[9]*primVelData[6] + dfdr3[10]*primVelData[7] + dfdr3[11]*primVelData[8];
    sdot[46] += dfdr3[12]*primVelData[6] + dfdr3[13]*primVelData[7] + dfdr3[14]*primVelData[8];
    sdot[47] += dfdr3[15]*primVelData[6] + dfdr3[16]*primVelData[7] + dfdr3[17]*primVelData[8];

    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================

/**
 *   \brief Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   \param t epoch at integration step
 *   \param s the 6-d state vector
 *   \param sdot the 6-d state derivative vector
 *   \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_bc4bp::simpleEOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_bc4bp *sysData = static_cast<const SysData_bc4bp *>(paramStruct->sysData);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    getPrimaryPos(t, sysData, primPosData);

    // Create relative position vectors between s/c and primaries
    double r_p1[] = {s[0] - primPosData[0], s[1] - primPosData[1], s[2] - primPosData[2]};
    double r_p2[] = {s[0] - primPosData[3], s[1] - primPosData[4], s[2] - primPosData[5]};
    double r_p3[] = {s[0] - primPosData[6], s[1] - primPosData[7], s[2] - primPosData[8]};
    double d1 = std::sqrt(r_p1[0]*r_p1[0] + r_p1[1]*r_p1[1] + r_p1[2]*r_p1[2]);
    double d2 = std::sqrt(r_p2[0]*r_p2[0] + r_p2[1]*r_p2[1] + r_p2[2]*r_p2[2]);
    double d3 = std::sqrt(r_p3[0]*r_p3[0] + r_p3[1]*r_p3[1] + r_p3[2]*r_p3[2]);
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Velocity
    std::copy(s+3, s+6, sdot);

    // Compute acceleration
    sdot[3] = 2*k*s[4] + k*k*(s[0] + 1/k - mu) - (1/k - mu)*(s[0] - primPosData[0])/pow(d1,3) -
        (mu-nu)*(s[0] - primPosData[3])/pow(d2,3) - nu*(s[0] - primPosData[6])/pow(d3,3);
    sdot[4] = -2*k*s[3] + k*k*s[1] - (1/k - mu)*(s[1] - primPosData[1])/pow(d1,3) -
        (mu-nu)*(s[1] - primPosData[4])/pow(d2,3) - nu*(s[1] - primPosData[7])/pow(d3,3);
    sdot[5] = -1*(1/k - mu)*(s[2] - primPosData[2])/pow(d1,3) -
        (mu-nu)*(s[2] - primPosData[5])/pow(d2,3) - nu*(s[2] - primPosData[8])/pow(d3,3);

    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================

/**
 *  \brief Compute the location of the three primaries in the BCR4BP (rotating coord.)
 *
 *  \param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  \param sysData a system data object containing information about the BCR4BP primaries
 *  \param primPos a pointer to a 1x9 double array that will hold the positions of the three primaries in 
 *  row-major order. The first three elements are the position of P1, etc.
 */
void DynamicsModel_bc4bp::getPrimaryPos(double t, const SysData_bc4bp *sysData, double *primPos){
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();
    double theta0 = sysData->getTheta0();
    double phi0 = sysData->getPhi0();
    double gamma = sysData->getGamma();
    double ratio = sysData->getCharLRatio();

    // Compute the angles for the system at the specified time
    double theta = theta0 + k*t;
    double phi = phi0 + sqrt(mu/pow(ratio, 3)) * t;

    // P1 position
    // primPos[0] = -mu;    // original derivation
    primPos[0] = -1/k;        // new derivation
    primPos[1] = 0;
    primPos[2] = 0;

    // P2 position
    // primPos[3] = 1/k - mu - nu/mu*ratio * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[3] = -nu/mu*ratio * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[4] = -nu/mu*ratio * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[5] = nu/mu*ratio * cos(phi) * sin(gamma);

    // P3 position
    // primPos[6] = 1/k - mu + (1 - nu/mu)*ratio * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[6] = (1 - nu/mu)*ratio * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[7] = (1 - nu/mu)*ratio * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[8] = (nu/mu - 1)*ratio * cos(phi)*sin(gamma);
}//================================================

/**
 *  \brief Compute the velocity of the three primaries in the BCR4BP, rotating coordinates.
 *
 *  \param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  \param sysData a system data object containing information about the BCR4BP primaries
 *  \param primVel a pointer to a 3x3 double array that will hold the velocities of the three primaries in
 *  row-major order. The first three elements hold the velocity of P1, etc.
 */
void DynamicsModel_bc4bp::getPrimaryVel(double t, const SysData_bc4bp *sysData, double *primVel){

    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();
    double theta0 = sysData->getTheta0();
    double phi0 = sysData->getPhi0();
    double gamma = sysData->getGamma();
    double ratio = sysData->getCharLRatio();

    double thetaDot = k;
    double phiDot = sqrt(mu/pow(ratio, 3));

    double theta = theta0 + thetaDot*t;
    double phi = phi0 + phiDot * t;    

    // P1 is stationary in this coordinate system
    primVel[0] = 0;  primVel[1] = 0;  primVel[2] = 0;

    // Angular velocity of P2-P3 line
    double v_P2P3Line[3] = {0};
    v_P2P3Line[0] = thetaDot*(sin(phi)*cos(theta) - cos(phi)*sin(theta)*cos(gamma)) + 
        phiDot*(cos(phi)*sin(theta) - sin(phi)*cos(theta)*cos(gamma));
    v_P2P3Line[1] = (phiDot - thetaDot)*cos(phi - theta);
    v_P2P3Line[2] = phiDot*sin(phi)*sin(gamma);

    // Multiply by radii of P2 and P3 to get their velocities
    primVel[3] = v_P2P3Line[0] * (-nu/mu)*ratio;
    primVel[4] = v_P2P3Line[1] * (-nu/mu)*ratio;
    primVel[5] = v_P2P3Line[2] * (-nu/mu)*ratio;

    primVel[6] = v_P2P3Line[0] * (1-nu/mu)*ratio;
    primVel[7] = v_P2P3Line[1] * (1-nu/mu)*ratio;
    primVel[8] = v_P2P3Line[2] * (1-nu/mu)*ratio;
}//===================================================================

/**
 *  \brief Compute the acceleration of the primary bodies in the BCR4BP, rotating coordinates
 * 
 *  \param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  \param sysData a system data object containing information about the BCR4BP primaries
 *  \param primAccel a pointer to a 3x3 double array that will hold the accelerations of the three primaries in
 *  row-major order. The first three elements hold the acceleration of P1, etc.
 */
void DynamicsModel_bc4bp::getPrimaryAccel(double t, const SysData_bc4bp *sysData, double *primAccel){
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();
    double theta0 = sysData->getTheta0();
    double phi0 = sysData->getPhi0();
    double gamma = sysData->getGamma();
    double ratio = sysData->getCharLRatio();

    double thetaDot = k;
    double phiDot = sqrt(mu/pow(ratio, 3));

    double theta = theta0 + thetaDot*t;
    double phi = phi0 + phiDot * t;    

    // P1 is stationary in this coordinate system
    primAccel[0] = 0;  primAccel[1] = 0;  primAccel[2] = 0;

    // Acceleration of P2-P3 line
    double a_P2P3Line[3] = {0};
    a_P2P3Line[0] = (-thetaDot*thetaDot - phiDot*phiDot)*(cos(theta)*cos(phi)*cos(gamma) + sin(theta)*sin(phi)) + 
        2*thetaDot*phiDot*(cos(theta)*cos(phi) + sin(theta)*sin(phi)*cos(gamma));
    a_P2P3Line[1] = -1*(phiDot - thetaDot)*(phiDot - thetaDot)*sin(phi - theta);
    a_P2P3Line[2] = phiDot*phiDot*cos(phi)*sin(gamma);

    // Multiply by radii of P2 and P3 to get their velocities
    primAccel[3] = a_P2P3Line[0] * (-nu/mu)*ratio;
    primAccel[4] = a_P2P3Line[1] * (-nu/mu)*ratio;
    primAccel[5] = a_P2P3Line[2] * (-nu/mu)*ratio;

    primAccel[6] = a_P2P3Line[0] * (1-nu/mu)*ratio;
    primAccel[7] = a_P2P3Line[1] * (1-nu/mu)*ratio;
    primAccel[8] = a_P2P3Line[2] * (1-nu/mu)*ratio;
}//=====================================================================

/**
 *  \brief Orient a Sun-Earth-Moon BCR4BPR system so that T = 0 corresponds to the specified epoch time
 *  
 *  \param et epoch time (seconds, J2000, UTC)
 *  \param sysData pointer to system data object; A new theta0 and phi0 will be stored
 *  in this data object
 */
void DynamicsModel_bc4bp::orientAtEpoch(double et, SysData_bc4bp *sysData){
    // Both theta and phi are approximately equal to zero at REF_EPOCH
    double time_nonDim = (et - SysData_bc4bp::REF_EPOCH)/sysData->getCharT();
    
    // Compute theta and phi
    double theta = sysData->getK()*time_nonDim;
    double phi = sqrt(sysData->getMu()/pow(sysData->getCharLRatio(), 3))*time_nonDim;

    // Leverage SPICE data to compute the true angles between my approximation and ephemeris
    double totalTheta = 0, totalPhi = 0, lt = 0;
    double moonState[6] = {0}, sunState[6] = {0};
    // Loop through one year of data, centered at the desired epoch
    for(unsigned int i = 0; i < 365; i++){
        double t = et + 3600*24*(i - 182); 
        double t_nondim = (t - SysData_bc4bp::REF_EPOCH)/sysData->getCharT();

        // Get ephemeris states for Sun and Moon relative to EMB in inertial, Ecliptic J2000 coordinates
        spkezr_c("MOON", t, "ECLIPJ2000", "NONE", "EMB", moonState, &lt);
        checkAndReThrowSpiceErr("DynamicsModel_bc4bp::orientAtEpoch error");

        spkezr_c("SUN", t, "ECLIPJ2000", "NONE", "EMB", sunState, &lt);
        checkAndReThrowSpiceErr("DynamicsModel_bc4bp::orientAtEpoch error");

        // Compute the angle between the Sun-EMB line and the inertial x-axis (shift by 90 because my alignment is with the y-axis)
        double spice_theta = atan2(-sunState[1], -sunState[0]) + PI/2; 
        // Compute the angle between the EMB-Moon line, projected into ecliptic plane, and the inertial x-axis (shift by 90 deg as above)
        double spice_phi = atan2(moonState[1], moonState[0]) + PI/2;

        // Compute approximate angles for the specified time
        double approx_theta = sysData->getK()*t_nondim;
        double approx_phi = sqrt(sysData->getMu()/pow(sysData->getCharLRatio(), 3))*t_nondim;

        // Sum the sin of the diference between the two. Using sin(angle) avoids numerical issues like 361 - 2 deg instead of 1 - 2
        totalTheta += sin(spice_theta - approx_theta);
        totalPhi += sin(spice_phi - approx_phi);

        // printf("Day %03d: dTheta = %.4f deg  dPhi = %.4f deg\n", i - 182,
        //     (spice_theta - approx_theta)*180/PI, (spice_phi - approx_phi)*180/PI);
    }

    // Adjust approximate theta and phi to be between 0 and 2*PI
    theta -= floor(theta/(2*PI))*2*PI;
    phi -= floor(phi/(2*PI))*2*PI;

    printf("Theta from SPICE = %.4f deg\n", (theta + asin(totalTheta/365))*180/PI);
    printf("Theta from BC4BP = %.4f deg\n", theta*180/PI);
    printf("Phi from SPICE = %.4f deg\n", (phi + asin(totalPhi/365))*180/PI);
    printf("Phi from BC4BP = %.4f deg\n", phi*180/PI);

    // Shift approximated angles by the average amount from the year-long survey conducted above
    theta += asin(totalTheta/365);
    phi += asin(totalPhi/365);

    sysData->setTheta0(theta);
    sysData->setPhi0(phi);
    sysData->setEpoch0(et);
}//====================================================




}// END of Astrohelion namespace