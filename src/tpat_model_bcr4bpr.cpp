/**
 *  @file tpat_model_bcr4bpr.cpp
 *  @brief Derivative of tpat_model, specific to BCR4BPR
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

#include "tpat_model_bcr4bpr.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_nodeset_bcr4bp.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_traj_bcr4bp.hpp"
#include "tpat_traj_step.hpp"
#include "tpat_event.hpp"
#include "tpat_node.hpp"
#include "tpat_utilities.hpp"

/**
 *  @brief Construct a BCR4BP Dynamic Model
 */
tpat_model_bcr4bpr::tpat_model_bcr4bpr() : tpat_model(MODEL_CR3BP) {
    coreStates = 6;
    stmStates = 36;
    extraStates = 6;
    allowedCons.push_back(tpat_constraint::SP);
    allowedCons.push_back(tpat_constraint::SP_RANGE);
}//==============================================

/**
 *  @brief Copy constructor
 *  @param m a model reference
 */
tpat_model_bcr4bpr::tpat_model_bcr4bpr(const tpat_model_bcr4bpr &m) : tpat_model(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
tpat_model_bcr4bpr& tpat_model_bcr4bpr::operator =(const tpat_model_bcr4bpr &m){
	tpat_model::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
tpat_model::eom_fcn tpat_model_bcr4bpr::getSimpleEOM_fcn(){
	return &bcr4bpr_simple_EOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
tpat_model::eom_fcn tpat_model_bcr4bpr::getFullEOM_fcn(){
	return &bcr4bpr_EOMs;
}//==============================================

/**
 *  @brief Compute the positions of all primaries
 *
 *  @param t the epoch at which the computations occur
 *  @param sysData object describing the specific system
 *  @return an n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> tpat_model_bcr4bpr::getPrimPos(double t, tpat_sys_data *sysData){
    double primPos[9];
    tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
    bcr4bpr_getPrimaryPos(t, bcSys, primPos);

    return std::vector<double>(primPos, primPos+9);
}//==============================================

/**
 *  @brief Compute the velocities of all primaries
 *
 *  @param t the epoch at which the computations occur
 *  @param sysData object describing the specific system
 *  @return an n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one velocity vector in non-dimensional units
 */
std::vector<double> tpat_model_bcr4bpr::getPrimVel(double t, tpat_sys_data *sysData){
    double primVel[9];
    tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
    bcr4bpr_getPrimaryVel(t, bcSys, primVel);

    return std::vector<double>(primVel, primVel+9);
}//==============================================

/**
 *  @brief Takes an input state and time and saves the data to the trajectory
 *  @param y an array containing the core state and any extra states integrated
 *  by the EOM function, including STM elements.
 *  @param t the time at the current integration state
 *  @param traj a pointer to the trajectory we should store the data in
 */
void tpat_model_bcr4bpr::sim_saveIntegratedData(double* y, double t, tpat_traj* traj){
    // Save the position and velocity states
    double state[6];
    std::copy(y, y+6, state);

    // Save STM
    double stmElm[36];
    std::copy(y+6, y+42, stmElm);

	// Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr*>(traj->getSysData());

    // Compute acceleration (elements 3-5)
    double dsdt[6] = {0};
    bcr4bpr_simple_EOMs(0, y, dsdt, bcSys);

    double dqdT[6];
    std::copy(y+42, y+48, dqdT);
    
    tpat_traj_step step(state, t, dsdt+3, stmElm);

    traj->appendStep(step);
    
    tpat_traj_bcr4bp *bcTraj = static_cast<tpat_traj_bcr4bp*>(traj);
    bcTraj->set_dqdT(-1, dqdT);
}//=====================================================

/**
 *  @brief Use a correction algorithm to accurately locate an event crossing
 *
 *  The simulation engine calls this function if and when it determines that an event 
 *  has been crossed. To accurately locate the event, we employ differential corrections
 *  and find the exact event occurence in space and time.
 *
 *  @param event the event we're looking for
 *  @param traj a pointer to the trajectory the event should occur on
 *  @param ic the core state vector for this system
 *  @param t0 non-dimensional time at the beginning of the search arc
 *  @param tof the time-of-flight for the arc to search over
 *  @param verbose whether or not we should be verbose with output messages
 *
 *  @return wether or not the event has been located. If it has, a new point
 *  has been appended to the trajectory's data vectors.
 */
bool tpat_model_bcr4bpr::sim_locateEvent(tpat_event event, tpat_traj *traj,
    double *ic, double t0, double tof, verbosity_t verbose){

    // **** Make sure you fix the epoch of the first node as well as the states
    double IC[7] = {0};
    std::copy(ic, ic+6, IC);
    IC[6] = t0;

    // Recast system data pointer
    tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr*>(traj->getSysData());

    // Create a nodeset for this particular type of system
    printVerb(verbose == ALL_MSG, "  Creating nodeset for event location\n");
    tpat_nodeset_bcr4bp eventNodeset(IC, bcSys, t0,
        tof, 2, tpat_nodeset::DISTRO_TIME);

    // Constraint to keep first node unchanged
    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, IC, 7);

    // Constraint to enforce event
    tpat_constraint eventCon(event.getConType(), 1, event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose == ALL_MSG){ eventNodeset.print(); }

    printVerb(verbose == ALL_MSG, "  Applying corrections process to locate event\n");
    tpat_correction_engine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbose(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations
    try{
        corrector.multShoot(&eventNodeset);
    }catch(tpat_diverge &e){
        printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(tpat_linalg_err &e){
        printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    tpat_nodeset_bcr4bp correctedNodes = corrector.getBCR4BPR_Output();

    std::vector<double> state = correctedNodes.getNode(-1).getPosVelState();
    std::vector<double> extra = correctedNodes.getNode(-1).getExtraParams();
    extra.insert(extra.begin(), state.begin(), state.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOF(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    sim_saveIntegratedData(&(extra[0]), eventTime, traj);
    
    return true;
}//=========================================================

/**
 *  @brief Initialize the corrector's design vector with position and velocity states,
 *  and times-of-flight.
 *
 *  Derived models may replace this function or call it and then append more design 
 *  variables.
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param set a pointer to the nodeset being corrected
 */
void tpat_model_bcr4bpr::multShoot_initDesignVec(iterationData *it, tpat_nodeset *set){
    // Call base class to do most of the work
    tpat_model::multShoot_initDesignVec(it, set);

    if(it->equalArcTime)
        throw tpat_exception("tpat_model_bcr4bpr::multShoot_initDesignVec: Equal Arc Times have not been implemented in this Model (yet)!");

    // Append the Epoch for each node
    if(it->varTime){
        // epochs come after ALL the TOFs have been added
        tpat_nodeset_bcr4bp *bcSet = static_cast<tpat_nodeset_bcr4bp *>(set);
        for(int n = 0; n < bcSet->getNumNodes(); n++){
            it->X.push_back(bcSet->getEpoch(n));
        }
    }
}//============================================================

/**
 *  @brief Create continuity constraints for the correction algorithm; this function
 *  creates position and velocity constraints.
 *
 *  This function overrides the base model's to add time continuity
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param set a pointer to the nodeset being corrected
 */ 
void tpat_model_bcr4bpr::multShoot_createContCons(iterationData *it, tpat_nodeset *set){
    tpat_model::multShoot_createContCons(it, set);

    if(it->varTime){
        std::vector<double> zero {0};
        for(int n = 1; n < set->getNumNodes(); n++){
            tpat_constraint timeCont(tpat_constraint::CONT_EX, n, zero);   // 0th index is epoch time in extraParam
            it->allCons.push_back(timeCont);
        }
    }
}//============================================================

/**
 *  @brief Retrieve the initial conditions for a segment that the correction
 *  engine will integrate.
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param set a pointer to the nodeset being corrected
 *  @param n the index of the node that serves as the initial state
 *  @param ic a pointer to a 6-element initial state array
 *  @param t0 a pointer to a double representing the initial time (epoch)
 *  @param tof a pointer to a double the time-of-flight on the segment.
 */
void tpat_model_bcr4bpr::multShoot_getSimICs(iterationData *it, tpat_nodeset *set, int n,
    double *ic, double *t0, double *tof){
    
    std::copy(it->X.begin()+6*n, it->X.begin()+6*(n+1), ic);
    *tof = it->varTime ? it->X[6*it->numNodes+n] : set->getTOF(n);
    tpat_nodeset_bcr4bp *bcSet = static_cast<tpat_nodeset_bcr4bp *>(set);
    *t0 = it->varTime ? it->X[7*it->numNodes-1+n] : bcSet->getEpoch(n);
}//============================================================

/**
 *  @brief Compute the value of a slack variable for an inequality constraint.
 *  @details Computing the value of the slack variable can avoid unneccessary 
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
 *  @param it the iterationData object associated with the multiple shooting process
 *  @param con the inequality constraint for which the slack variable is being computed
 * 
 *  @return The value of the slack variable that minimizes the constraint function
 *  without setting the slack variable to zero
 */
double tpat_model_bcr4bpr::multShoot_getSlackVarVal(iterationData *it, tpat_constraint con){
    switch(con.getType()){
        case tpat_constraint::SP_RANGE:
            return multShoot_targetSPMag_compSlackVar(it, con);
        default:
            return tpat_model::multShoot_getSlackVarVal(it, con);
    }
}//===========================================================

/**
 *  @brief Compute constraint function and partial derivative values for a constraint
 *  
 *  This function calls its relative in the tpat_model base class and appends additional
 *  instructions specific to the BCR4BPR
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void tpat_model_bcr4bpr::multShoot_applyConstraint(iterationData *it, tpat_constraint con, int c){

    // Let the base class do its thing first
    tpat_model::multShoot_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case tpat_constraint::SP:
            multShoot_targetSP(it, con, row0);
            break;
        case tpat_constraint::SP_RANGE:
            multShoot_targetSP_mag(it, con, c);
            break;
        default: break;
    }
}//=========================================================

/**
 *  @brief Compute position and velocity constraint values and partial derivatives
 *
 *  This function computes and stores the default position continuity constraints as well
 *  as velocity constraints for all nodes marked continuous in velocity. The delta-Vs
 *  between arc segments and node states are recorded, and the partial derivatives of each
 *  node with respect to other nodes and integration time are all computed
 *  and placed in the appropriate spots in the Jacobian matrix.
 *
 *  This function replaces the one found in the base model
 *
 *  @param it a pointer to the correctors iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the first row this constraint applies to
 */
void tpat_model_bcr4bpr::multShoot_targetPosVelCons(iterationData* it, tpat_constraint con, int row0){
    // Call base function first to do most of the work
    tpat_model::multShoot_targetPosVelCons(it, con, row0);

    // Add epoch dependencies for this model
    if(it->varTime){
        int n = con.getNode();
        std::vector<double> conData = con.getData();
        std::vector<double> last_dqdT = it->allSegs.at(n-1).getExtraParam(-1, 1);

        // Loop through conData
        for(size_t s = 0; s < conData.size(); s++){
            if(!isnan(conData[s])){
                // Epoch dependencies
                it->DF[it->totalFree*(row0+s) + 7*it->numNodes-1+n-1] = last_dqdT[s];
            }
        }
    }
}//=========================================================

/**
 *  @brief Computes continuity constraints for constraints with the <tt>CONT_EX</tt> type.
 *
 *  This function overrides the base model function
 *
 *  @param it a pointer to the correctors iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the first row this constraint applies to
 */
void tpat_model_bcr4bpr::multShoot_targetExContCons(iterationData *it, tpat_constraint con, int row0){
    int n = con.getNode();
    if(con.getData()[0] == 0){
        /* Add time-continuity constraints if applicable; we need to match
        the epoch time of node n to the sum of node n-1's epoch and TOF */
        if(it->varTime){
            it->FX[row0] = it->X[7*it->numNodes-1+n] - (it->X[7*it->numNodes-1+n-1] +
                it->X[6*it->numNodes+n-1]);
            it->DF[it->totalFree*(row0) + 6*it->numNodes+n-1] = -1;
            it->DF[it->totalFree*(row0) + 7*it->numNodes-1+n-1] = -1;
            it->DF[it->totalFree*(row0) + 7*it->numNodes-1+n] = 1;
        }
    }else{
        throw tpat_exception("tpat_model_bcr4bpr::multShoot_createExContCons: Unrecognized extra constraint index");
    }
}//=========================================================

/**
 *  @brief Compute partials and constraint functions for nodes constrained with <tt>STATE</tt>.
 *
 *  This method replaces the base class function and allows the user to constrain epoch as well
 *  as the configuration space states
 *
 *  @param it a pointer to the class containing all the data relevant to the corrections process
 *  @param con the constraint being applied
 *  @param row0 the index of the row this constraint begins at
 */
void tpat_model_bcr4bpr::multShoot_targetState(iterationData* it, tpat_constraint con, int row0){
    std::vector<double> conData = con.getData();
    int n = con.getNode();
    // Allow user to constrain all 7 states
    
    int count = 0;  // Count # rows since some may be skipped (NAN)
    for(int s = 0; s < ((int)con.getData().size()); s++){
        if(!isnan(conData[s])){
            if(s < 6){
                it->FX[row0+count] = it->X[6*n+s] - conData[s];
                it->DF[it->totalFree*(row0 + count) + 6*n + s] = 1;
                count++;
            }else if(s == 6){
                // Allow constraining epoch
                it->FX[row0+count] = it->X[7*it->numNodes-1+n] - conData[s];
                it->DF[it->totalFree*(row0 + count) + 7*it->numNodes-1+n] = 1;
                count++;
            }else{
                throw tpat_exception("State constraints must have <= 6 elements");
            }
        }
    }
}//=================================================

/**
 *  @brief Compute partials and constraint functions for nodes constrained with <tt>DIST</tt>, 
 *  <tt>MIN_DIST</tt>, or <tt>MAX_DIST</tt>
 *
 *  This method overrides the base class function to add functionality for epoch-time dependencies
 *
 *  @param it a pointer to the class containing all the data relevant to the corrections process
 *  @param con a copy of the constraint object
 *  @param c the index of this constraint in the constraint vector object
 */
void tpat_model_bcr4bpr::multShoot_targetDist(iterationData* it, tpat_constraint con, int c){

    std::vector<double> conData = con.getData();
    int n = con.getNode();
    int Pix = (int)(conData[0]);    // index of primary
    int row0 = it->conRows[c];
    
    // Get the node epoch either from the design vector or from the original set of nodes
    double t0 = it->varTime ? it->X[7*it->numNodes-1+n] : it->origNodes.at(n).getExtraParam(1);

    tpat_sys_data *sysData = it->sysData;

    // Get the primary position
    std::vector<double> primPos = getPrimPos(t0, sysData);

    // Get distance between node and primary in x, y, and z-coordinates
    double dx = it->X[6*n+0] - primPos[Pix*3+0];
    double dy = it->X[6*n+1] - primPos[Pix*3+1];
    double dz = it->X[6*n+2] - primPos[Pix*3+2];

    double h = sqrt(dx*dx + dy*dy + dz*dz);     // true distance

    // Compute difference between desired distance and true distance
    it->FX[row0] = h - conData[1];

    // Partials with respect to node position states
    it->DF[it->totalFree*row0 + 6*n + 0] = dx/h;
    it->DF[it->totalFree*row0 + 6*n + 1] = dy/h;
    it->DF[it->totalFree*row0 + 6*n + 2] = dz/h;

    if(it->varTime){
        // Epoch dependencies from primary positions
        double dhdr_data[] = {-dx/h, -dy/h, -dz/h};
        std::vector<double> primVel = getPrimVel(t0, sysData);
    
        Eigen::RowVector3d dhdr = Eigen::Map<Eigen::RowVector3d>(dhdr_data, 1, 3);
        Eigen::Vector3d drdT = Eigen::Map<Eigen::Vector3d>(&(primVel[Pix*3]), 3, 1);
        
        double prod = dhdr*drdT;
        it->DF[it->totalFree*row0 + 7*it->numNodes - 1 + n] = prod;
    }

    // Extra stuff for inequality constraints
    if(con.getType() == tpat_constraint::MIN_DIST || 
        con.getType() == tpat_constraint::MAX_DIST ){
        // figure out which of the slack variables correspond to this constraint
        std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
            it->slackAssignCon.end(), c);

        // which column of the DF matrix the slack variable is in
        int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
        int sign = con.getType() == tpat_constraint::MAX_DIST ? 1 : -1;

        // printf("Dist from P%d is %f (%s %f)\n", Pix, h,
        //     con.getType() == tpat_constraint::MIN_DIST ? "Min" : "Max", conData[1]);
        // printf("  Slack Var^2 = %e\n", it->X[slackCol]*it->X[slackCol]);
        // Subtract squared slack variable from constraint
        it->FX[row0] += sign*it->X[slackCol]*it->X[slackCol];

        // Partial with respect to slack variable
        it->DF[it->totalFree*row0 + slackCol] = sign*2*it->X[slackCol];
    }
}// End of targetDist() =========================================

/**
 *  @brief  Compute the value of the slack variable for inequality distance constraints
 *  @details This function computes a value for the slack variable in an
 *  inequality distance constraint. If the constraint is already met by the initial
 *  design, using this value will prevent the multiple shooting algorithm from
 *  searching all over for the propper value.
 *  
 *  Note: This method overrides the base class function to add functionality for non-zero epochs
 * 
 *  @param it the iteration data object for the multiple shooting process
 *  @param con the constraint the slack variable applies to
 *  @return the value of the slack variable that minimizes the constraint function
 *  without setting the slack variable equal to zero
 */
double tpat_model_bcr4bpr::multShoot_targetDist_compSlackVar(iterationData* it, tpat_constraint con){
    std::vector<double> conData = con.getData();
    int n = con.getNode();
    int Pix = (int)(conData[0]);    // index of primary 
    // Get the node epoch either from the design vector or from the original set of nodes
    double t0 = it->varTime ? it->X[7*it->numNodes-1+n] : it->origNodes.at(n).getExtraParam(1);
    tpat_sys_data *sysData = it->sysData;

    // Get the primary position
    std::vector<double> primPos = getPrimPos(t0, sysData);

    // Get distance between node and primary in x, y, and z-coordinates
    double dx = it->X[6*n+0] - primPos[Pix*3+0];
    double dy = it->X[6*n+1] - primPos[Pix*3+1];
    double dz = it->X[6*n+2] - primPos[Pix*3+2];

    double h = sqrt(dx*dx + dy*dy + dz*dz);     // true distance
    int sign = con.getType() == tpat_constraint::MAX_DIST ? 1 : -1;
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
 *  @brief Compute partials and constraints for all nodes constrained with <tt>DELTA_V</tt> or
 *  <tt>MIN_DELTA_V</tt>
 *
 *  Because the delta-V constraint applies to the entire trajectory, the constraint function values
 *  and partial derivatives must be computed for each node along the trajectory. This function
 *  takes care of all of them at once.
 *
 *  This function overrides the base targeting function to add support for epoch dependencies
 *
 *  @param it a pointer to the class containing all the data relevant to the corrections process
 *  @param con the constraint being applied
 *  @param c the index of the first row for this constraint
 */
void tpat_model_bcr4bpr::multShoot_targetDeltaV(iterationData* it, tpat_constraint con, int c){
    // Call base function to take care of most of the constraint computations and partials
    tpat_model::multShoot_targetDeltaV(it, con, c);

    if(it->varTime){
        // Add partials w.r.t. epoch time
        int row0 = it->conRows[c];

        // Compute total deltaV magnitude
        for(int n = 0; n < it->numNodes-1; n++){
            // Compute parial w.r.t. node n+1 (where velocity is discontinuous)
            double dFdq_n2_data[] = {0, 0, 0, -2*it->deltaVs[n*3], -2*it->deltaVs[n*3+1], -2*it->deltaVs[n*3+2]};
            Eigen::RowVectorXd dFdq_n2 = Eigen::Map<Eigen::RowVectorXd>(dFdq_n2_data, 1, 6);

            // Compute partial w.r.t. epoch time n
            std::vector<double> last_dqdT = it->allSegs.at(n).getExtraParam(-1, 1);

            Eigen::VectorXd dqdT = Eigen::Map<Eigen::VectorXd>(&(last_dqdT[0]), 6, 1);
            MatrixXRd dFdT_n;
            dFdT_n.noalias() = -1*dFdq_n2*dqdT;
            it->DF[it->totalFree*row0 + 7*it->numNodes-1+n] = dFdT_n(0);
        }
    }
}//==============================================

/**
 *  @brief Compute partials and constraint values for nodes constrained with <tt>SP</tt>
 *
 *  This function computes three constraint values and three rows of partials for the Jacobian.
 *  Each row/function corresponds to one position state. The FX and DF matrices are updated
 *  in place by editing their values stored in <tt>it</tt>
 *
 *  @param it the iterationData object holding the current data for the corrections process
 *  @param con the constraint being applied
 *  @param row0 the index of the first row for this constraint
 */
void tpat_model_bcr4bpr::multShoot_targetSP(iterationData* it, tpat_constraint con, int row0){

    int n = con.getNode();
    double epoch = it->varTime ? it->X[7*it->numNodes-1+n] : it->origNodes.at(n).getExtraParam(1);
    tpat_sys_data_bcr4bpr *bcSysData = static_cast<tpat_sys_data_bcr4bpr *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(epoch, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double *X = &(it->X[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+6*n, 3, 1);   // Position vector

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

    MatrixXRd dFdr2 = Eigen::Map<MatrixXRd>(dfdr2_data, 6, 3);
    MatrixXRd dFdr3 = Eigen::Map<MatrixXRd>(dfdr3_data, 6, 3);

    // scale matrices by constants
    dFdr2 *= -1*(mu - nu);
    dFdr3 *= -1*nu;

    // Compute partials of constraint function w.r.t. epoch time
    MatrixXRd dFdT;
    dFdT.noalias() = dFdr2*(primVel.row(1).transpose()) + dFdr3*(primVel.row(2).transpose());

    // Copy data into the correct vectors/matrices
    double* conEvalPtr = conEval.data();
    double* dFdq_ptr = dFdq.data();
    double* dFdT_ptr = dFdT.data();

    double *FX = &(it->FX[0]);
    double *DF = &(it->DF[0]);

    std::copy(conEvalPtr, conEvalPtr+3, FX+row0);
    std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + 6*n);
    std::copy(dFdq_ptr+3, dFdq_ptr+6, DF + it->totalFree*(row0+1) + 6*n);
    std::copy(dFdq_ptr+6, dFdq_ptr+9, DF + it->totalFree*(row0+2) + 6*n);

    if(it->varTime){
        std::copy(dFdT_ptr, dFdT_ptr+1, DF + it->totalFree*row0 + 7*it->numNodes-1+n);
        std::copy(dFdT_ptr+1, dFdT_ptr+2, DF + it->totalFree*(row0+1) + 7*it->numNodes-1+n);
        std::copy(dFdT_ptr+2, dFdT_ptr+3, DF + it->totalFree*(row0+2) + 7*it->numNodes-1+n);
    }
}// End of SP Targeting ==============================

void tpat_model_bcr4bpr::multShoot_targetSP_mag(iterationData* it, tpat_constraint con, int c){

    int row0 = it->conRows[c];
    int n = con.getNode();
    double Amax = con.getData()[0];
    double epoch = it->varTime ? it->X[7*it->numNodes-1+n] : it->origNodes.at(n).getExtraParam(1);
    tpat_sys_data_bcr4bpr *bcSysData = static_cast<tpat_sys_data_bcr4bpr *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(epoch, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double *X = &(it->X[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+6*n, 3, 1);   // Position vector

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
    
    std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + 6*n);

    if(it->varTime){
        std::copy(dFdT_ptr, dFdT_ptr+1, DF + it->totalFree*row0 + 7*it->numNodes-1+n);
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

double tpat_model_bcr4bpr::multShoot_targetSPMag_compSlackVar(iterationData *it, tpat_constraint con){
    int n = con.getNode();
    double Amax = con.getData()[0];
    double epoch = it->varTime ? it->X[7*it->numNodes-1+n] : it->origNodes.at(n).getExtraParam(1);
    tpat_sys_data_bcr4bpr *bcSysData = static_cast<tpat_sys_data_bcr4bpr *> (it->sysData);

    std::vector<double> primPosData = getPrimPos(epoch, it->sysData);

    // Get primary positions at the specified epoch time
    Matrix3Rd primPos = Eigen::Map<Matrix3Rd>(&(primPosData[0]), 3, 3);

    double *X = &(it->X[0]);
    Eigen::Vector3d r = Eigen::Map<Eigen::Vector3d>(X+6*n, 3, 1);   // Position vector

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
    return diff > 0 ? sqrt(std::abs(diff)) : 1e-4;
}//======================================================

/**
 *  @brief Take the final, corrected free variable vector <tt>X</tt> and create an output 
 *  nodeset
 *
 *  If <tt>findEvent</tt> is set to true, the
 *  output nodeset will contain extra information for the simulation engine to use. Rather than
 *  returning only the position and velocity states, the output nodeset will contain the STM 
 *  and dqdT values for the final node; this information will be appended to the extraParameter
 *  vector in the final node.
 *
 *  @param it an iteration data object containing all info from the corrections process
 *  @param nodes_in a pointer to the original, uncorrected nodeset
 *  @param findEvent whether or not this correction process is locating an event
 *
 *  @return a pointer to a nodeset containing the corrected nodes
 */
tpat_nodeset* tpat_model_bcr4bpr::multShoot_createOutput(iterationData *it, tpat_nodeset *nodes_in, bool findEvent){

    // Create a nodeset with the same system data as the input
    tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(it->sysData);
    tpat_nodeset_bcr4bp *nodeset_out = new tpat_nodeset_bcr4bp(bcSys);

    int numNodes = (int)(it->origNodes.size());
    for(int i = 0; i < numNodes; i++){
        tpat_node node(&(it->X[i*6]), NAN);
        node.setExtraParam(1, it->X[7*numNodes-1 + i]);     // save epoch time
        node.setVelCon(nodes_in->getNode(i).getVelCon());   // save velocity continuity info
        node.setConstraints(nodes_in->getNode(i).getConstraints());    // save constraints

        if(i + 1 < numNodes){
            node.setTOF(it->X[numNodes*6 + i]);
        }else{
            node.setTOF(NAN);
            
            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
                state for the last node. We do this by appending the STM elements and dqdT elements to the
                end of the node array. This output nodeset should have two "nodes": the first 6 elements
                are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
                information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                tpat_traj lastSeg = it->allSegs.back();
                MatrixXRd lastSTM = lastSeg.getSTM(-1);
                
                // Create a vector of extra parameters from existing extraParam vector
                std::vector<double> extraParams = nodeset_out->getNode(-1).getExtraParams();
                // append the STM elements at the end
                extraParams.insert(extraParams.end(), lastSTM.data(), lastSTM.data()+36);
                // append the last dqdT vector
                std::vector<double> dqdT = lastSeg.getExtraParam(-1,1);
                extraParams.insert(extraParams.end(), dqdT.begin(), dqdT.end());
                node.setExtraParams(extraParams);
            }
        }

        nodeset_out->appendNode(node);
    }

    return nodeset_out;
}//====================================================




