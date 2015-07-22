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

tpat_model_cr3bp::tpat_model_cr3bp() : tpat_model(MODEL_CR3BP) {}
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
    std::vector<double> *nodes = correctedNodes.getNodes();

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOF(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    model->saveIntegratedData(&(nodes->at(6)), eventTime, traj);

    return true;
}//======================================================