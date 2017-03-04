/**
 *  @file DynamicsModel_cr3bp_lt.cpp
 *  @brief Derivative of DynamicsModel, specific to CR3BP-LTVP
 *  
 *  @author Andrew Cox
 *  @version May 25, 2016
 *  @copyright GNU GPL v3.0
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

#include "DynamicsModel_cr3bp_lt.hpp"

#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "CorrectionEngine.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj_cr3bp_lt.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{
/**
 *  @brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic DynamicsModel
 */
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP_LT) {
    coreStates = 7;
    extraStates = 0;
    allowedCons.push_back(Constraint_tp::JC);
    allowedEvents.push_back(Event_tp::JC);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */ 
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt &m) : DynamicsModel(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
DynamicsModel_cr3bp_lt& DynamicsModel_cr3bp_lt::operator =(const DynamicsModel_cr3bp_lt &m){
	DynamicsModel::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_lt::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_lt::getFullEOM_fcn() const{
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
std::vector<double> DynamicsModel_cr3bp_lt::getPrimPos(double t, const SysData *sysData) const{
    (void)t;
    double primPos[6] = {0};
    const SysData_cr3bp_lt *crSys = static_cast<const SysData_cr3bp_lt *>(sysData);

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
std::vector<double> DynamicsModel_cr3bp_lt::getPrimVel(double t, const SysData *sysData) const{
    (void)t;
    (void)sysData;
    double primVel[6] = {0};
    
    return std::vector<double>(primVel, primVel+6);
}//==============================================

std::vector<double> DynamicsModel_cr3bp_lt::getAccel(const SysData *pSys, double t, std::vector<double> state) const{
    if(state.size() != coreStates)
        throw Exception("DynamicsModel_cr3bp_lt::getAccel: State size does not match the core state size specified by the dynamical model");

    // Compute the acceleration
    std::vector<double> dsdt(coreStates, 0);
    EOM_ParamStruct paramStruct(pSys, ControlLaw::NO_CTRL);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), &paramStruct);
    
    return std::vector<double>(dsdt.begin()+3, dsdt.end());
}//==================================================

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
void DynamicsModel_cr3bp_lt::sim_saveIntegratedData(const double* y, double t, Traj* traj) const{
    
    DynamicsModel::sim_saveIntegratedData(y, t, traj);

    // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const SysData_cr3bp_lt *ltSys = static_cast<const SysData_cr3bp_lt*>(traj->getSysData());
    Traj_cr3bp_lt *ltTraj = static_cast<Traj_cr3bp_lt*>(traj);

    // Save Jacobi for CR3BP - it won't be constant any more, but is definitely useful to have
    ltTraj->setJacobiByIx(-1, DynamicsModel_cr3bp::getJacobi(y, ltSys->getMu()));
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
bool DynamicsModel_cr3bp_lt::sim_locateEvent(Event event, Traj* traj,
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
void DynamicsModel_cr3bp_lt::multShoot_createOutput(const MultShootData *it, const Nodeset *nodes_in, bool findEvent, Nodeset *nodes_out) const{
    (void) it;
    (void) findEvent;
    (void) nodes_in;
    (void) nodes_out;
}//====================================================

/**
 *  @brief Perform model-specific initializations on the MultShootData object
 *  @param it pointer to the object to be initialized
 */
void DynamicsModel_cr3bp_lt::multShoot_initIterData(MultShootData *it) const{
    Traj_cr3bp_lt traj(static_cast<const SysData_cr3bp_lt *>(it->sysData));
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
int DynamicsModel_cr3bp_lt::fullEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();
    double f_nondim = T*charT*charT/charL;

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );
    double v = sqrt( s[3]*s[3] + s[4]*s[4] + s[5]*s[5] );   // Rotating Velocity

    double thrust_dir[3];
    const ControlLaw_cr3bp_lt* law = static_cast<const ControlLaw_cr3bp_lt *>(sysData->getControlLaw());
    law->getLaw(t, s, sysData, paramStruct->ctrlLawID, thrust_dir, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[2];
    sdot[6] = -T*charT/(Isp*G_GRAV_0);

    /*
     * Next step, compute STM
     */
    
    // Partials of x_ddot w.r.t. state variables
    // double dxdx = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((x + mu),2)/pow(d,5) + 
    //     3*mu*pow((x + mu - 1), 2)/pow(r,5) + (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    // double dxdy = 3*(1-mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) +
    //     (T/m)*(v*v + (xdot - y)*(xdot - y))/pow(v,3);
    // double dxdz = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5);
    // double dxdxdot = -1*(T/m)*(v*v - (xdot - y)*(xdot - y))/pow(v,3);
    // double dxdydot = (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    // double dxdzdot = (T/m)*(xdot - y)*zdot/pow(v,3);

    // // Partials of y_ddot w.r.t. state variables
    // double dydx = 3*(1 - mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) - 
    //     (T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    // double dydy = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*y*y/pow(d,5) + 3*mu*y*y/pow(r,5) +
    //     +(T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    // double dydz = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5);
    // double dydxdot = -2 + (T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    // double dydydot = -1*(T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    // double dydzdot = (T/m)*(ydot + x)*zdot/pow(v,3);

    // // Partials of z_ddot w.r.t. state variables
    // double dzdx = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5) + (T/m)*(ydot + x)*zdot/pow(v,3);
    // double dzdy = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5) - (T/m)*(xdot - y)*zdot/pow(v,3);
    // double dzdz = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*z*z/pow(d,5) + 3*mu*z*z/pow(r,5);
    // double dzdxdot = (T/m)*zdot*(xdot - y)/pow(v,3);
    // double dzdydot = (T/m)*zdot*(ydot + x)/pow(v,3);
    // double dzdzdot = -1*(T/m)*(v*v - zdot*zdot)/pow(v,3);

    // // Create A Matrix
    // double a_data[] = { 0,    0,    0,    1,       0,       0,
    //                     0,    0,    0,    0,       1,       0,
    //                     0,    0,    0,    0,       0,       1,
    //                     dxdx, dxdy, dxdz, dxdxdot, dxdydot, dxdzdot,
    //                     dydx, dydy, dydz, dydxdot, dydydot, dydzdot,
    //                     dzdx, dzdy, dzdz, dzdxdot, dzdydot, dzdzdot};
    // MatrixXRd A = Eigen::Map<MatrixXRd>(a_data, 6, 6);
    MatrixXRd A = MatrixXRd::Identity(7,7);

    // Copy the STM states into a sub-array
    double stmElements[49];
    std::copy(s+7, s+56, stmElements);

    // Turn sub-array into matrix object for math stuffs
    MatrixXRd phi = Eigen::Map<MatrixXRd>(stmElements, 7, 7);

    // Compute derivative of STM
    MatrixXRd phiDot(7,7);
    phiDot.noalias() = A*phi;     // use noalias() to avoid creating an unnecessary temporary matrix in Eigen library

    // Copy the elements of phiDot into the derivative array
    double *phiDotData = phiDot.data();
    std::copy(phiDotData, phiDotData+49, sdot+7);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 7-d state vector
 *  @param sdot the 7-d state derivative vector
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::simpleEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->sysData);

    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();
    double f_nondim = T*charT*charT/charL;

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );
    double v = sqrt( s[3]*s[3] + s[4]*s[4] + s[5]*s[5] );   // Rotating Velocity

    double thrust_dir[3];
    const ControlLaw_cr3bp_lt* law = static_cast<const ControlLaw_cr3bp_lt *>(sysData->getControlLaw());
    law->getLaw(t, s, sysData, paramStruct->ctrlLawID, thrust_dir, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + f_nondim/(s[6]*v)*thrust_dir[2];
    sdot[6] = -T*charT/(Isp*G_GRAV_0);

    return GSL_SUCCESS;
}//=====================================================


}// END of Astrohelion namespace