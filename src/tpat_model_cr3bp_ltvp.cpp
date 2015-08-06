/**
 *  @file tpat_model_cr3bp_ltvp.cpp
 *  @brief Derivative of tpat_model, specific to CR3BP-LTVP
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

#include "tpat_model_cr3bp_ltvp.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"
#include "tpat_event.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"

/**
 *  @brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic Model
 */
tpat_model_cr3bp_ltvp::tpat_model_cr3bp_ltvp() : tpat_model(MODEL_CR3BP_LTVP) {
    coreStates = 6;
    stmStates = 36;
    extraStates = 0;
    allowedCons.push_back(tpat_constraint::JC);
    allowedEvents.push_back(tpat_event::JC);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */ 
tpat_model_cr3bp_ltvp::tpat_model_cr3bp_ltvp(const tpat_model_cr3bp_ltvp &m) : tpat_model(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
tpat_model_cr3bp_ltvp& tpat_model_cr3bp_ltvp::operator =(const tpat_model_cr3bp_ltvp &m){
	tpat_model::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
tpat_model::eom_fcn tpat_model_cr3bp_ltvp::getSimpleEOM_fcn(){
	return &cr3bp_ltvp_simple_EOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
tpat_model::eom_fcn tpat_model_cr3bp_ltvp::getFullEOM_fcn(){
	return &cr3bp_ltvp_EOMs;
}//==============================================

/**
 *  @brief Compute the positions of all primaries
 *
 *  @param t the epoch at which the computations occur (unused for this system)
 *  @param sysData object describing the specific system
 *  @return an n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> tpat_model_cr3bp_ltvp::getPrimPos(double t, tpat_sys_data *sysData){
    (void)t;
    double primPos[6] = {0};
    tpat_sys_data_cr3bp_ltvp crSys(*static_cast<tpat_sys_data_cr3bp_ltvp *>(sysData));

    primPos[0] = -1*crSys.getMu();
    primPos[3] = 1 - crSys.getMu();

    return std::vector<double>(primPos, primPos+6);
}//==============================================

/**
 *  @brief Compute the velocities of all primaries
 *
 *  @param t the epoch at which the computations occur (unused for this system)
 *  @param sysData object describing the specific system (unused for this system)
 *  @return an n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one velocity vector in non-dimensional units
 */
std::vector<double> tpat_model_cr3bp_ltvp::getPrimVel(double t, tpat_sys_data *sysData){
    (void)t;
    (void)sysData;
    double primVel[6] = {0};
    
    return std::vector<double>(primVel, primVel+6);
}//==============================================

/**
 *  @brief Takes an input state and time and saves the data to the trajectory
 *  @param y an array containing the core state and any extra states integrated
 *  by the EOM function, including STM elements.
 *  @param t the time at the current integration state
 *  @param traj a pointer to the trajectory we should store the data in
 */
void tpat_model_cr3bp_ltvp::sim_saveIntegratedData(double* y, double t, tpat_traj* traj){
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
    tpat_traj_cr3bp_ltvp *cr3bpTraj = static_cast<tpat_traj_cr3bp_ltvp*>(traj);
    tpat_sys_data_cr3bp_ltvp sysData = cr3bpTraj->getSysData();

    // Compute acceleration
    double dsdt[6] = {0};
    cr3bp_ltvp_simple_EOMs(t, y, dsdt, &sysData);

    // Save the accelerations
    traj->getAccel()->push_back(dsdt[3]);
    traj->getAccel()->push_back(dsdt[4]);
    traj->getAccel()->push_back(dsdt[5]);

    // Save Jacobi for CR3BP - it won't be constant any more, but is definitely useful to have
    cr3bpTraj->getJacobi()->push_back(cr3bp_getJacobi(y, sysData.getMu()));

    // Compute and save mass of s/c; assumes t began at 0
    double g0_nonDim = G_GRAV_0*sysData.getCharT()*sysData.getCharT()/sysData.getCharL();
    cr3bpTraj->getMass()->push_back(sysData.getM0() - sysData.getThrust()/(sysData.getIsp()*g0_nonDim) * t);
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
 *  @param model the dynamical model we're working in
 *  @param ic the core state vector for this system
 *  @param t0 non-dimensional time at the beginning of the search arc
 *  @param tof the time-of-flight for the arc to search over
 *  @param verbose whether or not we should be verbose with output messages
 *
 *  @return wether or not the event has been located. If it has, a new point
 *  has been appended to the trajectory's data vectors.
 */
bool tpat_model_cr3bp_ltvp::sim_locateEvent(tpat_event event, tpat_traj* traj, tpat_model* model,
    double *ic, double t0, double tof, bool verbose){

    return true;
}//=======================================================

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
 */
tpat_nodeset* tpat_model_cr3bp_ltvp::corrector_createOutput(iterationData *it, tpat_nodeset *nodes_in, bool findEvent){

    return new tpat_nodeset_cr3bp;
}//====================================================