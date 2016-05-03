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

#include "tpat_model_cr3bp_ltvp.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_event.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"
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
tpat_model::eom_fcn tpat_model_cr3bp_ltvp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
tpat_model::eom_fcn tpat_model_cr3bp_ltvp::getFullEOM_fcn() const{
	return &fullEOMs;
}//==============================================

/**
 *  @brief Compute the positions of all primaries
 *
 *  @param t the epoch at which the computations occur (unused for this system)
 *  @param sysData object describing the specific system
 *  @return an n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> tpat_model_cr3bp_ltvp::getPrimPos(double t, const tpat_sys_data *sysData) const{
    (void)t;
    double primPos[6] = {0};
    const tpat_sys_data_cr3bp_ltvp crSys(*static_cast<const tpat_sys_data_cr3bp_ltvp *>(sysData));

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
std::vector<double> tpat_model_cr3bp_ltvp::getPrimVel(double t, const tpat_sys_data *sysData) const{
    (void)t;
    (void)sysData;
    double primVel[6] = {0};
    
    return std::vector<double>(primVel, primVel+6);
}//==============================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Takes an input state and time and saves the data to the trajectory
 *  @param y an array containing the core state and any extra states integrated
 *  by the EOM function, including STM elements.
 *  @param t the time at the current integration state
 *  @param traj a pointer to the trajectory we should store the data in
 */
void tpat_model_cr3bp_ltvp::sim_saveIntegratedData(const double* y, double t, tpat_traj* traj) const{
    // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const tpat_sys_data_cr3bp_ltvp *ltSys = static_cast<const tpat_sys_data_cr3bp_ltvp*>(traj->getSysData());

    // Compute acceleration (elements 3-5)
    double dsdt[6] = {0};
    eomParamStruct paramStruct(ltSys);
    simpleEOMs(t, y, dsdt, &paramStruct);
    
    // node(state, accel, epoch) - y(0:5) holds the state, y(6:41) holds the STM
    int id = traj->addNode(tpat_node(y, dsdt+3, t));

    if(id > 0){
        double tof = t - traj->getNode(id-1).getEpoch();
        traj->addSeg(tpat_segment(id-1, id, tof, y+6));
    }

    tpat_traj_cr3bp_ltvp *cr3bpTraj = static_cast<tpat_traj_cr3bp_ltvp*>(traj);

    // Save Jacobi for CR3BP - it won't be constant any more, but is definitely useful to have
    cr3bpTraj->setJacobiByIx(-1, tpat_model_cr3bp::getJacobi(y, ltSys->getMu()));

    // Compute and save mass of s/c; assumes t began at 0
    double g0_nonDim = G_GRAV_0*ltSys->getCharT()*ltSys->getCharT()/ltSys->getCharL();
    cr3bpTraj->setMassByIx(-1, ltSys->getM0() - ltSys->getThrust()/(ltSys->getIsp()*g0_nonDim) * t);
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
bool tpat_model_cr3bp_ltvp::sim_locateEvent(tpat_event event, tpat_traj* traj,
    const double *ic, double t0, double tof, tpat_verbosity_tp verbose) const{

    (void) event;
    (void) traj;
    (void) ic;
    (void) t0;
    (void) tof;
    (void) verbose;

    return true;
}//=======================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

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
tpat_nodeset* tpat_model_cr3bp_ltvp::multShoot_createOutput(const iterationData *it, const tpat_nodeset *nodes_in, bool findEvent) const{
    (void) it;
    (void) findEvent;
    const tpat_sys_data_cr3bp_ltvp *sys = static_cast<const tpat_sys_data_cr3bp_ltvp*>(nodes_in->getSysData());
    return new tpat_nodeset_cr3bp(sys);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP
 *  @param t the current time of the integration
 *  @param s the 43-d state vector. The first 6 elements are position and velocity,
 *  the 7th is mass, and the final 36 are STM elements
 *  @param sdot the 43-d state derivative vector
 *  @param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an eomParamStruct object
 */
int tpat_model_cr3bp_ltvp::fullEOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params
    // tpat_sys_data_cr3bp_ltvp *sysData = static_cast<tpat_sys_data_cr3bp_ltvp *>(params);
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const tpat_sys_data_cr3bp_ltvp *sysData = static_cast<const tpat_sys_data_cr3bp_ltvp *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    
    double g0_nonDim = (G_GRAV_0/charL)*charT*charT;
    double m = sysData->getM0() - T/(Isp*g0_nonDim)*t;    // assumes t began at 0

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );
    double v = sqrt( (xdot - y)*(xdot - y) + (ydot + x)*(ydot + x) + zdot*zdot);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3) - T/(m*v)*(xdot - y);
    sdot[4] = -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3) - T/(m*v)*(ydot + x);
    sdot[5] = -(1-mu)*z/pow(d,3) - mu*z/pow(r,3) - T/(m*v)*zdot; 

    /*
     * Next step, compute STM
     */
    
    // Partials of x_ddot w.r.t. state variables
    double dxdx = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((x + mu),2)/pow(d,5) + 
        3*mu*pow((x + mu - 1), 2)/pow(r,5) + (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    double dxdy = 3*(1-mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) +
        (T/m)*(v*v + (xdot - y)*(xdot - y))/pow(v,3);
    double dxdz = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5);
    double dxdxdot = -1*(T/m)*(v*v - (xdot - y)*(xdot - y))/pow(v,3);
    double dxdydot = (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    double dxdzdot = (T/m)*(xdot - y)*zdot/pow(v,3);

    // Partials of y_ddot w.r.t. state variables
    double dydx = 3*(1 - mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) - 
        (T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    double dydy = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*y*y/pow(d,5) + 3*mu*y*y/pow(r,5) +
        +(T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    double dydz = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5);
    double dydxdot = -2 + (T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    double dydydot = -1*(T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    double dydzdot = (T/m)*(ydot + x)*zdot/pow(v,3);

    // Partials of z_ddot w.r.t. state variables
    double dzdx = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5) + (T/m)*(ydot + x)*zdot/pow(v,3);
    double dzdy = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5) - (T/m)*(xdot - y)*zdot/pow(v,3);
    double dzdz = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*z*z/pow(d,5) + 3*mu*z*z/pow(r,5);
    double dzdxdot = (T/m)*zdot*(xdot - y)/pow(v,3);
    double dzdydot = (T/m)*zdot*(ydot + x)/pow(v,3);
    double dzdzdot = -1*(T/m)*(v*v - zdot*zdot)/pow(v,3);

    // Create A Matrix
    double a_data[] = { 0,    0,    0,    1,       0,       0,
                        0,    0,    0,    0,       1,       0,
                        0,    0,    0,    0,       0,       1,
                        dxdx, dxdy, dxdz, dxdxdot, dxdydot, dxdzdot,
                        dydx, dydy, dydz, dydxdot, dydydot, dydzdot,
                        dzdx, dzdy, dzdz, dzdxdot, dzdydot, dzdzdot};
    MatrixXRd A = Eigen::Map<MatrixXRd>(a_data, 6, 6);

    // Copy the STM states into a sub-array
    double stmElements[36];
    std::copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    MatrixXRd phi = Eigen::Map<MatrixXRd>(stmElements, 6, 6);

    // Compute derivative of STM
    MatrixXRd phiDot(6,6);
    phiDot.noalias() = A*phi;     // use noalias() to avoid creating an unnecessary temporary matrix in Eigen library

    // Copy the elements of phiDot into the derivative array
    double *phiDotData = phiDot.data();
    std::copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 7-d state vector
 *  @param sdot the 7-d state derivative vector
 *  @param params points to an eomParamStruct object
 */
int tpat_model_cr3bp_ltvp::simpleEOMs(double t, const double s[], double sdot[], void *params){
    // tpat_sys_data_cr3bp_ltvp *sysData = static_cast<tpat_sys_data_cr3bp_ltvp *>(params);
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const tpat_sys_data_cr3bp_ltvp *sysData = static_cast<const tpat_sys_data_cr3bp_ltvp *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    
    double g0_nonDim = G_GRAV_0/charL*charT*charT;
    double m = sysData->getM0() - T/(Isp*g0_nonDim)*t;    // assumes t began at 0

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );
    double v = sqrt( (xdot - y)*(xdot - y) + (ydot + x)*(ydot + x) + zdot*zdot);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3) - T/(m*v)*(xdot - y);
    sdot[4] = -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3) - T/(m*v)*(ydot + x);
    sdot[5] = -(1-mu)*z/pow(d,3) - mu*z/pow(r,3) - T/(m*v)*zdot; 

    return GSL_SUCCESS;
}//=====================================================
