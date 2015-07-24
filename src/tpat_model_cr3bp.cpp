/**
 *  @file tpat_model_cr3bp.cpp
 *
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

#include "tpat_model_cr3bp.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_event.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"

tpat_model_cr3bp::tpat_model_cr3bp() : tpat_model(MODEL_CR3BP) {
    // Allow a few more constraints than the default
    allowedCons.push_back(tpat_constraint::JC);
}//==============================================

tpat_model_cr3bp::tpat_model_cr3bp(const tpat_model_cr3bp &m) : tpat_model(m) {}

tpat_model_cr3bp& tpat_model_cr3bp::operator =(const tpat_model_cr3bp &m){
	tpat_model::operator =(m);
	return *this;
}//==============================================

tpat_model::eom_fcn tpat_model_cr3bp::getSimpleEOM_fcn(){
	return &cr3bp_simple_EOMs;
}//==============================================

tpat_model::eom_fcn tpat_model_cr3bp::getFullEOM_fcn(){
	return &cr3bp_EOMs;
}//==============================================

std::vector<double> tpat_model_cr3bp::getPrimPos(double t, tpat_sys_data *sysData){
    double primPos[6] = {0};
    tpat_sys_data_cr3bp crSys(*static_cast<tpat_sys_data_cr3bp *>(sysData));
    
    primPos[0] = -1*crSys.getMu();
    primPos[3] = 1 - crSys.getMu();

    return std::vector<double>(primPos, primPos+6);
}//==============================================

std::vector<double> tpat_model_cr3bp::getPrimVel(double t, tpat_sys_data *sysData){
    double primVel[6] = {0};
    
    return std::vector<double>(primVel, primVel+6);
}//==============================================

void tpat_model_cr3bp::saveIntegratedData(double* y, double t, tpat_traj* traj){
    // Save the position and velocity states
    for(int i = 0; i < 6; i++){
        traj->getState()->push_back(y[i]);
    }

    // Save time
    traj->getTime()->push_back(t);

    // Save STM
    double stmElm[36];
    std::copy(y+6, y+42, stmElm);
    traj->getSTM()->push_back(tpat_matrix(6,6,stmElm));

	// Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    tpat_traj_cr3bp *cr3bpTraj = static_cast<tpat_traj_cr3bp*>(traj);
    tpat_sys_data_cr3bp sysData = cr3bpTraj->getSysData();

    // Compute acceleration
    double dsdt[6] = {0};
    cr3bp_simple_EOMs(0, y, dsdt, &sysData);

    // Save the accelerations
    traj->getAccel()->push_back(dsdt[3]);
    traj->getAccel()->push_back(dsdt[4]);
    traj->getAccel()->push_back(dsdt[5]);
    
    cr3bpTraj->getJacobi()->push_back(cr3bp_getJacobi(y, sysData.getMu()));
}//=====================================================

/**
 *  @brief Use a correction algorithm to check and see if an event has occurred
 *
 *  The simulation engine calls this function if and when it determines that an event 
 *  has been crossed. To accurately locate the event, we employ differential corrections
 *  and find the exact event occurence in space and time.
 *
 *  @param event the event we're looking for
 *  @param sysData a system data object describing the system we're working in
 *  @param model the dynamical model we're working in
 *  @param ic the core state vector for this system
 *  @param t0 non-dimensional time at the beginning of the search arc
 *  @param tof the time-of-flight for the arc to search over
 *  @param verbose whether or not we should be verbose with output messages
 *
 *  @param return wether or not the event has been located. If it has, a new point
 *  has been appended to the trajectory's data vectors.
 */
bool tpat_model_cr3bp::locateEvent(tpat_event event, tpat_traj* traj, tpat_model* model,
    double *ic, double t0, double tof, bool verbose){

    // Copy system data object
    tpat_traj_cr3bp *crTraj = static_cast<tpat_traj_cr3bp *>(traj);

    // Create a nodeset for this particular type of system
    printVerb(verbose, "  Creating nodeset for event location\n");
    tpat_nodeset_cr3bp eventNodeset(ic, crTraj->getSysData(), tof, 2, tpat_nodeset::DISTRO_TIME);

    // Constraint to keep first node unchanged
    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, ic, 6);

    // Constraint to enforce event
    tpat_constraint eventCon(event.getConType(), event.getConNode(), event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose){ eventNodeset.print(); }

    printVerb(verbose, "  Applying corrections process to locate event\n");
    tpat_correction_engine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbose(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations
    try{
        corrector.correct_cr3bp(&eventNodeset);
    }catch(tpat_diverge &e){
        printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(tpat_linalg_err &e){
        printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    tpat_nodeset_cr3bp correctedNodes = corrector.getCR3BP_Output();

    std::vector<double> state = correctedNodes.getNode(-1).getPosVelState();
    std::vector<double> extra = correctedNodes.getNode(-1).getExtraParams();
    extra.insert(extra.begin(), state.begin(), state.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOF(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    model->saveIntegratedData(&(extra[0]), eventTime, traj);

    return true;
}//======================================================


/**
 *  @brief Compute constraint function and partial derivative values for a constraint
 *  
 *  This function calls its relative in the tpat_model base class and appends additional
 *  instructions specific to the CR3BP
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void tpat_model_cr3bp::corrector_applyConstraint(iterationData *it, tpat_constraint con, int c){

    // Let the base class do its thing first
    tpat_model::corrector_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case tpat_constraint::JC:
            corrector_targetJC(it, con, row0);
            break;
        default: break;
    }
}//=========================================================

void tpat_model_cr3bp::corrector_targetJC(iterationData* it, tpat_constraint con, int row0){
    std::vector<double> conData = con.getData();
    int n = con.getNode();
    tpat_sys_data_cr3bp *crSys = static_cast<tpat_sys_data_cr3bp *> (it->sysData);

    // Compute the value of Jacobi at this node
    double mu = crSys->getMu();
    double *nodeState = &(it->X[6*n]);
    double nodeJC = cr3bp_getJacobi(nodeState, mu);
    
    // temp variables to make equations more readable; compute partials w.r.t. node state
    double x = nodeState[0];
    double y = nodeState[1];
    double z = nodeState[2];
    double vx = nodeState[3];
    double vy = nodeState[4];
    double vz = nodeState[5];

    double d = sqrt((x + mu)*(x + mu) + y*y + z*z);
    double r = sqrt((x + mu - 1)*(x + mu - 1) + y*y + z*z);

    it->FX[row0] = nodeJC - conData[0];
    it->DF[it->totalFree*row0 + 6*n + 0] = -2*(x + mu)*(1 - mu)/pow(d,3) - 2*(x + mu - 1)*mu/pow(r,3) + 2*x;    //dFdx
    it->DF[it->totalFree*row0 + 6*n + 1] = -2*y*(1 - mu)/pow(d,3) - 2*y*mu/pow(r,3) - 2*y;                      //dFdy
    it->DF[it->totalFree*row0 + 6*n + 2] = -2*z*(1 - mu)/pow(d,3) - 2*z*mu/pow(r,3);                            //dFdz
    it->DF[it->totalFree*row0 + 6*n + 3] = -2*vx;   //dFdx_dot
    it->DF[it->totalFree*row0 + 6*n + 4] = -2*vy;   //dFdy_dot
    it->DF[it->totalFree*row0 + 6*n + 5] = -2*vz;   //dFdz_dot
}//=============================================

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
 *  @param nodeset_out a pointer to the output nodeset object
 *  @param findEvent whether or not this correction process is locating an event
 */
tpat_nodeset* tpat_model_cr3bp::corrector_createOutput(iterationData *it, bool findEvent){

    // Create a nodeset with the same system data as the input
    tpat_sys_data_cr3bp *crSys = static_cast<tpat_sys_data_cr3bp *>(it->sysData);
    tpat_nodeset_cr3bp *nodeset_out = new tpat_nodeset_cr3bp(*crSys);
    nodeset_out->setNodeDistro(tpat_nodeset::DISTRO_NONE);

    int numNodes = (int)(it->origNodes.size());
    for(int i = 0; i < numNodes; i++){
        if(i + 1 < numNodes){
            tpat_node node(&(it->X[i*6]), it->X[numNodes*6 + i]);   // create node with state and TOF
            nodeset_out->appendNode(node);
        }else{
            tpat_node node(&(it->X[i*6]), NAN);                 // create node with state and fake TOF (last node)
            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
            state for the last node. We do this by appending the STM elements and dqdT elements to the
            end of the node array. This output nodeset should have two "nodes": the first 6 elements
            are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
            information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                tpat_traj lastSeg = it->allSegs.back();
                tpat_matrix stm = lastSeg.getSTM(-1);
                std::vector<double> extraParam(stm.getDataPtr(), stm.getDataPtr()+36);
                
                node.setExtraParams(extraParam);
            }
            nodeset_out->appendNode(node);
        }
    }


    return nodeset_out;
}//====================================================



