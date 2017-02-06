/**
 *  @file DynamicsModel_cr3bp_ltvp.cpp
 *  @brief Derivative of DynamicsModel, specific to CR3BP-LTVP
 *  
 *  @author Andrew Cox
 *  @version May 25, 2016
 *  @copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "DynamicsModel_cr3bp_ltvp.hpp"

#include "Calculations.hpp"
#include "CorrectionEngine.hpp"
#include "Event.hpp"
#include "MultShootData.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SysData_cr3bp_ltvp.hpp"
#include "Traj_cr3bp_ltvp.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{
/**
 *  @brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic DynamicsModel
 */
DynamicsModel_cr3bp_ltvp::DynamicsModel_cr3bp_ltvp() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP_LTVP) {
    coreStates = 6;
    stmStates = 36;
    extraStates = 0;
    allowedCons.push_back(Constraint_tp::JC);
    allowedEvents.push_back(Event_tp::JC);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */ 
DynamicsModel_cr3bp_ltvp::DynamicsModel_cr3bp_ltvp(const DynamicsModel_cr3bp_ltvp &m) : DynamicsModel(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
DynamicsModel_cr3bp_ltvp& DynamicsModel_cr3bp_ltvp::operator =(const DynamicsModel_cr3bp_ltvp &m){
	DynamicsModel::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_ltvp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_ltvp::getFullEOM_fcn() const{
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
std::vector<double> DynamicsModel_cr3bp_ltvp::getPrimPos(double t, const SysData *sysData) const{
    (void)t;
    double primPos[6] = {0};
    const SysData_cr3bp_ltvp *crSys = static_cast<const SysData_cr3bp_ltvp *>(sysData);

    primPos[0] = -1*crSys->getMu();
    primPos[3] = 1 - crSys->getMu();

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
std::vector<double> DynamicsModel_cr3bp_ltvp::getPrimVel(double t, const SysData *sysData) const{
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
void DynamicsModel_cr3bp_ltvp::sim_saveIntegratedData(const double* y, double t, Traj* traj) const{
    // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const SysData_cr3bp_ltvp *ltSys = static_cast<const SysData_cr3bp_ltvp*>(traj->getSysData());

    // Compute acceleration (elements 3-5)
    double dsdt[6] = {0};
    EOM_ParamStruct paramStruct(ltSys);
    simpleEOMs(t, y, dsdt, &paramStruct);
    
    // node(state, accel, epoch) - y(0:5) holds the state, y(6:41) holds the STM
    int id = traj->addNode(Node(y, dsdt+3, t));

    if(id > 0){
        double tof = t - traj->getNode(id-1).getEpoch();
        traj->addSeg(Segment(id-1, id, tof, y+6));
    }

    Traj_cr3bp_ltvp *cr3bpTraj = static_cast<Traj_cr3bp_ltvp*>(traj);

    // Save Jacobi for CR3BP - it won't be constant any more, but is definitely useful to have
    cr3bpTraj->setJacobiByIx(-1, DynamicsModel_cr3bp::getJacobi(y, ltSys->getMu()));

    // Compute and save mass of s/c; assumes t began at 0
    double g0_nonDim = G_GRAV_0*ltSys->getCharT()*ltSys->getCharT()/ltSys->getCharL();
    // cr3bpTraj->setMassByIx(-1, ltSys->getM0() - ltSys->getThrust()/(ltSys->getIsp()*g0_nonDim) * t);
    throw Exception("Need to update mass in CR3BP LTVP saveIntegratedData!");
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
bool DynamicsModel_cr3bp_ltvp::sim_locateEvent(Event event, Traj* traj,
    const double *ic, double t0, double tof, Verbosity_tp verbose) const{

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
 *  @param nodes_out pointer to the nodeset object that will contain the output of the
 *  shooting process
 */
void DynamicsModel_cr3bp_ltvp::multShoot_createOutput(const MultShootData *it, const Nodeset *nodes_in, bool findEvent, Nodeset *nodes_out) const{
    (void) it;
    (void) findEvent;
    (void) nodes_in;
    (void) nodes_out;
}//====================================================

/**
 *  @brief Perform model-specific initializations on the MultShootData object
 *  @param it pointer to the object to be initialized
 */
void DynamicsModel_cr3bp_ltvp::multShoot_initIterData(MultShootData *it) const{
    Traj_cr3bp_ltvp traj(static_cast<const SysData_cr3bp_ltvp *>(it->sysData));
    it->propSegs.assign(it->nodeset->getNumSegs(), traj);
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
 *  function, the pointer points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_ltvp::fullEOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params
    // SysData_cr3bp_ltvp *sysData = static_cast<SysData_cr3bp_ltvp *>(params);
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_ltvp *sysData = static_cast<const SysData_cr3bp_ltvp *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    double m = s[6];    // Mass

    double g0_nonDim = (G_GRAV_0/charL)*charT*charT;

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );
    double v = sqrt( (xdot - y)*(xdot - y) + (ydot + x)*(ydot + x) + zdot*zdot);    // Velocity in an inertial frame

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
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_ltvp::simpleEOMs(double t, const double s[], double sdot[], void *params){
    // SysData_cr3bp_ltvp *sysData = static_cast<SysData_cr3bp_ltvp *>(params);
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_ltvp *sysData = static_cast<const SysData_cr3bp_ltvp *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    
    double g0_nonDim = G_GRAV_0/charL*charT*charT;
    // double m = sysData->getM0() - T/(Isp*g0_nonDim)*t;    // assumes t began at 0
    double m = 0;
    throw Exception("Need to update mass in CR3BP LTVP EOMs!");

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


}// END of Astrohelion namespace