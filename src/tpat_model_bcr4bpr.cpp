/**
 *  @file tpat_model_bcr4bpr.cpp
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

#include "tpat_model_bcr4bpr.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_traj_bcr4bpr.hpp"
#include "tpat_event.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"

tpat_model_bcr4bpr::tpat_model_bcr4bpr() : tpat_model(MODEL_CR3BP) {
    coreStates = 6;
    stmStates = 36;
    extraStates = 6;
}//==============================================

tpat_model_bcr4bpr::tpat_model_bcr4bpr(const tpat_model_bcr4bpr &m) : tpat_model(m) {}

tpat_model_bcr4bpr& tpat_model_bcr4bpr::operator =(const tpat_model_bcr4bpr &m){
	tpat_model::operator =(m);
	return *this;
}//==============================================

tpat_model::eom_fcn tpat_model_bcr4bpr::getSimpleEOM_fcn(){
	return &bcr4bpr_simple_EOMs;
}//==============================================

tpat_model::eom_fcn tpat_model_bcr4bpr::getFullEOM_fcn(){
	return &bcr4bpr_EOMs;
}//==============================================


void tpat_model_bcr4bpr::saveIntegratedData(double* y, double t, tpat_traj* traj){
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
    tpat_traj_bcr4bpr *bcr4bprTraj = static_cast<tpat_traj_bcr4bpr*>(traj);
    tpat_sys_data_bcr4bpr sysData = bcr4bprTraj->getSysData();

    // Compute acceleration
    double dsdt[6] = {0};
    bcr4bpr_simple_EOMs(0, y, dsdt, &sysData);

    // Save the accelerations
    traj->getAccel()->push_back(dsdt[3]);
    traj->getAccel()->push_back(dsdt[4]);
    traj->getAccel()->push_back(dsdt[5]);
    
    // Save dqdT
    bcr4bprTraj->get_dqdT()->push_back(y[42]);
    bcr4bprTraj->get_dqdT()->push_back(y[43]);
    bcr4bprTraj->get_dqdT()->push_back(y[44]);
    bcr4bprTraj->get_dqdT()->push_back(y[45]);
    bcr4bprTraj->get_dqdT()->push_back(y[46]);
    bcr4bprTraj->get_dqdT()->push_back(y[47]);
}//=====================================================


bool tpat_model_bcr4bpr::locateEvent(tpat_event event, tpat_traj *traj, tpat_model* model,
    double *ic, double t0, double tof, bool verbose){

    // **** Make sure you fix the epoch of the first node as well as the states
    double IC[7] = {0};
    std::copy(ic, ic+6, IC);
    IC[6] = t0;

    // Copy system data object
    tpat_traj_bcr4bpr *bcTraj = static_cast<tpat_traj_bcr4bpr *>(traj);

    // Create a nodeset for this particular type of system
    printVerb(verbose, "  Creating nodeset for event location\n");
    tpat_nodeset_bcr4bpr eventNodeset(IC, bcTraj->getSysData(), t0,
        tof, 2, tpat_nodeset::DISTRO_TIME);

    // Constraint to keep first node unchanged
    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, IC, 7);

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
        corrector.correct_bcr4bpr(&eventNodeset);
    }catch(tpat_diverge &e){
        printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(tpat_linalg_err &e){
        printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    tpat_nodeset_bcr4bpr correctedNodes = corrector.getBCR4BPR_Output();
    std::vector<double> *nodes = correctedNodes.getNodes();

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOF(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    model->saveIntegratedData(&(nodes->at(6)), eventTime, traj);
    
    return true;
}