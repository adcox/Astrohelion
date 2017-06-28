/**
 *  \file DynamicsModel_cr3bp_lt.cpp
 *  \brief Derivative of DynamicsModel, specific to CR3BP-LTVP
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

#include "DynamicsModel_cr3bp_lt.hpp"

#include "Arcset_cr3bp_lt.hpp"
#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "MultShootEngine.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{

/**
 *  \brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic DynamicsModel
 */
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt() : DynamicsModel_cr3bp() {
    coreDim = 7;
    extraDim = 0;
    allowedEvents.push_back(Event_tp::MASS);
}//==============================================

/**
 *  \brief Copy Constructor
 *  \param m a model reference
 */ 
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt &m) : DynamicsModel_cr3bp(m) {}

/**
 *  \brief Assignment operator
 *  \param m a model reference
 */
DynamicsModel_cr3bp_lt& DynamicsModel_cr3bp_lt::operator =(const DynamicsModel_cr3bp_lt &m){
	DynamicsModel::operator =(m);
	return *this;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_lt::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp_lt::getFullEOM_fcn() const{
	return &fullEOMs;
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
std::vector<double> DynamicsModel_cr3bp_lt::getStateDeriv(double t, std::vector<double> state, EOM_ParamStruct *params) const{
    const unsigned int ctrlDim = params->pCtrlLaw ? params->pCtrlLaw->getNumStates() : 0;

    if(state.size() != coreDim + ctrlDim)
        throw Exception("DynamicsModel_cr3bp_lt::getStateDeriv: State size does not match the state size specified by the dynamical model and control law");

    // Compute the acceleration
    std::vector<double> dsdt(coreDim + ctrlDim, 0);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), params);
    
    return dsdt;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Create default events for a simulation run
 *  \details These events are intended to prevent numerical issues, e.g., to avoid singularities.
 * 
 *  \param pSys pointer to system data object
 *  \return A vector of events to use in the simulation
 */
std::vector<Event> DynamicsModel_cr3bp_lt::sim_makeDefaultEvents(const SysData *pSys) const{
    // Create crash events from base dynamics model
    std::vector<Event> events = DynamicsModel::sim_makeDefaultEvents(pSys);

    // Add event to keep mass greater than 0.01 (1% of spacecraft reference mass)
    std::vector<double> minMass {0.01};
    events.push_back(Event(Event_tp::MASS, -1, true, minMass));

    return events;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Perform model-specific initializations on the MultShootData object
 *  \param it pointer to the object to be initialized
 */
void DynamicsModel_cr3bp_lt::multShoot_initIterData(MultShootData *it) const{
    Arcset_cr3bp_lt traj(static_cast<const SysData_cr3bp_lt *>(it->nodesIn->getSysData()));
    it->propSegs.assign(it->nodesIn->getNumSegs(), traj);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Integrate the equations of motion for the CR3BP LTVP
 *  \param t the current time of the integration
 *  \param s the state vector passed in from the SimEngine. This vector includes
 *  the core states, STM states, extra states, and control states, in that order.
 *  \param sdot the 43-d state derivative vector
 *  \param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::fullEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->pSysData);
    const ControlLaw_cr3bp_lt *law = static_cast<const ControlLaw_cr3bp_lt *>(paramStruct->pCtrlLaw);
    const unsigned int coreDim = sysData->getDynamicsModel()->getCoreStateSize();

    double mu = sysData->getMu();           // nondimensional mass ratio

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Retrieve the control law acceleration values
    double control_accel[3] = {0};
    if(law){
        law->getLaw_Output(t, s, sysData, control_accel, 3);
    }

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + control_accel[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + control_accel[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + control_accel[2];

    // nondimensional mass flow rate; simplified a bit by cancelling some of the constants
    if(law)
        sdot[6] = -(law->getThrust()/1000)*sysData->getCharT()/(sysData->getRefMass()*law->getIsp()*G_GRAV_0);
    else
        sdot[6] = 0;

    // Save any time-derivatives of the control states
    const unsigned int ctrlDim = law ? law->getNumStates() : 0;
    const unsigned int ctrlOutDim = law ? law->getNumOutputs() : 0;
    if(ctrlDim > 0){
        std::vector<double> control_stateDeriv(ctrlDim, 0);
        law->getLaw_StateDeriv(t, s, sysData, &(control_stateDeriv.front()), ctrlDim);

        std::copy(control_stateDeriv.begin(), control_stateDeriv.begin() + ctrlDim, sdot + coreDim);
    }

    
    const unsigned int stmSide = (coreDim + ctrlDim);
    MatrixXRd A = MatrixXRd::Zero(stmSide, stmSide);

    // Velocity relationships
    A(0, 3) = 1;    // d/dvx (dx/dt)
    A(1, 4) = 1;    // d/dvy (dy/dt)
    A(2, 5) = 1;    // d/dvz (dz/dt)

    // Uxx = d/dx (dvx/dt)
    A(3, 0) = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((s[0] + mu),2)/pow(d,5) + 
        3*mu*pow((s[0] + mu - 1), 2)/pow(r,5);
    // Uxy = d/dy (dvx/dt)
    A(3, 1) = 3*(1-mu)*(s[0] + mu)*s[1]/pow(d,5) + 3*mu*(s[0] + mu - 1)*s[1]/pow(r,5);
    // Uxz = d/dz (dvx/dt)
    A(3, 2) = 3*(1-mu)*(s[0] + mu)*s[2]/pow(d,5) + 3*mu*(s[0] + mu - 1)*s[2]/pow(r,5);

    // Uyy = d/dy (dvy/dt)
    A(4, 1) = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*s[1]*s[1]/pow(d,5) + 3*mu*s[1]*s[1]/pow(r,5);

    // Uyz = d/dz (dvy/dt)
    A(4, 2) = 3*(1-mu)*s[1]*s[2]/pow(d,5) + 3*mu*s[1]*s[2]/pow(r,5);

    // Uzz = d/dz (dvz/dt)
    A(5, 2) = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*s[2]*s[2]/pow(d,5) + 3*mu*s[2]*s[2]/pow(r,5);

    // Symmetry
    A(4,0) = A(3,1);
    A(5,0) = A(3,2);
    A(5,1) = A(4,2);

    A(3, 4) = 2;    // d/dvy (dvx/dt)
    A(4, 3) = -2;   // d/dvx (dvy/dt)

    if(law){
        // Get partial derivatives of acceleration terms (which are part of EOMS 3, 4, 5) w.r.t. all state variables
        std::vector<double> law_accelPartials(ctrlOutDim*coreDim, 0);
        
        law->getLaw_OutputPartials(t, s, sysData, &(law_accelPartials.front()), law_accelPartials.size());

        // Add the control output partials to the existing partials (control outputs are part of core state EOMs)
        for(unsigned int r = 3; r < 3+ctrlOutDim; r++){
            for(unsigned int c = 0; c < coreDim; c++){
                A(r, c) += law_accelPartials.at((r - 3)*coreDim + c);
            }
        }

        if(ctrlDim > 0){
            std::vector<double> law_eomPartials(coreDim*ctrlDim, 0), law_stateDerivPartials(ctrlDim*(coreDim+ctrlDim), 0);

            law->getLaw_StateDerivPartials(t, s, sysData, &(law_stateDerivPartials.front()), law_stateDerivPartials.size());
            law->getLaw_EOMPartials(t, s, sysData, &(law_eomPartials.front()), law_eomPartials.size());

            // Assign the partial derivatives of the time derivatives of the control states w.r.t. all core and control states
            for(unsigned int r = coreDim; r < coreDim + ctrlDim; r++){
                for(unsigned int c = 0; c < coreDim + ctrlDim; c++){
                    A(r, c) = law_stateDerivPartials.at((r - coreDim)*(coreDim + ctrlDim) + c);
                }
            }

            // Assign the partial derivatives of the EOMs w.r.t. control states
            for(unsigned int r = 0; r < coreDim; r++){
                for(unsigned int c = coreDim; c < coreDim + ctrlDim; c++){
                    A(r, c) = law_eomPartials.at(r*ctrlDim + c - coreDim);
                }
            }
        }
    }

    // Turn sub-array into matrix object for math stuffs
    const MatrixXRd phi = Eigen::Map<const MatrixXRd>(s + coreDim + ctrlDim, stmSide, stmSide);

    // Compute derivative of STM
    MatrixXRd phiDot(stmSide, stmSide);
    phiDot.noalias() = A*phi;     // use noalias() to avoid creating an unnecessary temporary matrix in Eigen library

    // Copy the elements of phiDot into the derivative array
    const double *phiDotData = phiDot.data();
    std::copy(phiDotData, phiDotData + stmSide*stmSide, sdot + coreDim + ctrlDim);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  \brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  \param t time at integration step (unused)
 *  \param s the state vector passed in from the SimEngine. This vector includes
 *  the core states and control states, in that order.
 *  \param sdot the state derivative vector
 *  \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::simpleEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->pSysData);
    const ControlLaw_cr3bp_lt *law = static_cast<const ControlLaw_cr3bp_lt *>(paramStruct->pCtrlLaw);

    double mu = sysData->getMu();           // nondimensional mass ratio

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Retrieve the control law acceleration values
    double control_accel[3] = {0};
    if(law)
        law->getLaw_Output(t, s, sysData, control_accel, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + control_accel[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + control_accel[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + control_accel[2];
    
    // nondimensional mass flow rate; simplified a bit by cancelling some of the constants
    sdot[6] = -(law->getThrust()/1000)*sysData->getCharT()/(sysData->getRefMass()*law->getIsp()*G_GRAV_0);

    // Save any derivatives of the control states
    if(unsigned int ctrl_dim = law->getNumStates() > 0){
        std::vector<double> control_stateDeriv(ctrl_dim, 0);
        law->getLaw_StateDeriv(t, s, sysData, &(control_stateDeriv.front()), ctrl_dim);

        std::copy(control_stateDeriv.begin(), control_stateDeriv.begin() + ctrl_dim, sdot+7);
    }

    return GSL_SUCCESS;
}//=====================================================


/**
 *  \brief Construct a new control law and allocated it on the stack.
 *  \details Each dynamic model will return a pointer to the specific control
 *  law applicable to the system / model
 *  \return A pointer to a control law object. The object has been allocated
 *  on the stack so the delete() function must be employed to free the memory
 */
ControlLaw* DynamicsModel_cr3bp_lt::createControlLaw() const{ return new ControlLaw_cr3bp_lt; }

/**
 *  \brief Determine whether or not the model supports a specific control law
 *  \details By default, a DynamicsModel does not support any control laws.
 * 
 *  \param pLaw Pointer to a control law
 *  \return Whether or not the Dynamics model supports the specified control law
 */
bool DynamicsModel_cr3bp_lt::supportsControl(const ControlLaw *pLaw) const{
    // Does support control laws
    (void) pLaw;
    return true;
}//================================================





}// END of Astrohelion namespace


