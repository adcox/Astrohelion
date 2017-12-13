/**
 *  @file DynamicsModel_2bp.cpp
 *  @brief Derivative of DynamicsModel, specific to 2BP
 *  
 *  @author Andrew Cox
 *  @version August 24, 2016
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

#include "DynamicsModel_2bp.hpp"

#include "Calculations.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "Event.hpp"
#include "Node.hpp"
#include "Segment.hpp"
#include "SysData_2bp.hpp"
#include "Arcset_2bp.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{

DynamicsModel_2bp::DynamicsModel_2bp() : DynamicsModel() {}

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */
DynamicsModel_2bp::DynamicsModel_2bp(const DynamicsModel_2bp &m) : DynamicsModel(m) {

}//====================================================

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
DynamicsModel_2bp& DynamicsModel_2bp::operator =(const DynamicsModel_2bp &m){
	DynamicsModel::operator =(m);
	return *this;
}//====================================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_2bp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//====================================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_2bp::getFullEOM_fcn() const{
	return &fullEOMs;
}//====================================================

/**
 *  @brief Compute the position of the primary body
 * 
 *  @param t the epoch at which the computations occur (unused for this system)
 *  @param pSysData object describing the specific system (unused for this system)
 * 
 *  @return n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> DynamicsModel_2bp::getPrimPos(double t, const SysData *pSysData) const{
	(void)t;
	(void)pSysData;
	return std::vector<double>(3,0);
}//====================================================

/**
 *  @brief Compute the position of a specified primary
 *  @details This is the faster alternative to getPrimPos(t, pSysData).
 * 
 *  @param t Nondimensional time
 *  @param pSysData pointer to system data object
 *  @param pIx Index of the primary; a value of -1 will return the positions of all primaries,
 *  in order of largest to smallest mass
 *  @param pos An array to store the primary position(s) in with all elements initialized to zero.
 *  For a single primary position, the array must have at least three elements allocated. For all 
 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
 *  of primaries.
 */
void DynamicsModel_2bp::getPrimPos(double t, const SysData *pSysData, int pIx, double *pos) const{
    (void) t;
    (void) pSysData;
    (void) pIx;
    (void) pos; // Always zero, nothing to do here...
}//====================================================

/**
 *  @brief Compute the velocity of the primary body
 * 
 *  @param t the epoch at which the computations occur (unused for this system)
 *  @param pSysData object describing the specific system (unused for this system)
 * 
 *  @return n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> DynamicsModel_2bp::getPrimVel(double t, const SysData *pSysData) const{
	(void)t;
	(void)pSysData;
	return std::vector<double>(3,0);
}//====================================================

/**
 *  @brief Compute the velocity of a specified primary
 *  @details This is the faster alternative to getPrimVel(t, pSysData).
 * 
 *  @param t Nondimensional time
 *  @param pSysData pointer to system data object
 *  @param pIx Index of the primary; a value of -1 will return the velocities of all primaries,
 *  in order of largest to smallest mass
 *  @param vel An array to store the primary velocity(s) in with all elements initialized to zero. 
 *  For a single primary velocity, the array must have at least three elements allocated. For all 
 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
 *  of primaries.
 */
void DynamicsModel_2bp::getPrimVel(double t, const SysData *pSysData, int pIx, double *vel) const{
    (void) t;
    (void) pSysData;
    (void) pIx;
    (void) vel;
}//====================================================
/**
 *  @brief Retrieve the state derivative
 *  @details Evaluate the equations of motion to compute the state time-derivative at 
 *  the specified time and state
 * 
 *  @param t time parameter
 *  @param state state vector
 *  @param params structure containing parameters relevant to the integration
 *  @return the time-derivative of the state vector
 */
std::vector<double> DynamicsModel_2bp::getStateDeriv(double t, std::vector<double> state, EOM_ParamStruct *params) const{
    const unsigned int ctrlDim = params->pCtrlLaw ? params->pCtrlLaw->getNumStates() : 0;

    if(state.size() != coreDim + ctrlDim)
        throw Exception("DynamicsModel_2bp::getStateDeriv: State size does not match the state size specified by the dynamical model and control law");

    // Compute the acceleration
    std::vector<double> dsdt(coreDim + ctrlDim, 0);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), params);
    
    return dsdt;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

    // Use the default Base class functions

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Perform model-specific initializations on the MultShootData object
 *  @details NOT YET IMPLEMENTED
 *  
 *  @param it reference to the object to be initialized
 */
void DynamicsModel_2bp::multShoot_initIterData(MultShootData& it) const{
	(void)it;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the 2BP
 *  @param t time at integration step (unused)
 *  @param s the 42-d state vector
 *  @param sdot the 42-d state derivative vector
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_2bp::fullEOMs(double t, const double s[], double sdot[], void *params){
	(void)t;
	
	EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_2bp *pSys = static_cast<const SysData_2bp *>(paramStruct->pSysData);

    double r = sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] );
    double mult = -pSys->getMu()/pow(r,3);	// G*(mass_primary)/(r^3)

    // Position derivatives = velocity
    std::copy(s+3, s+6, sdot);

    // Velocity derivatives = acceleration
    sdot[3] = mult*s[0];
    sdot[4] = mult*s[1];
    sdot[5] = mult*s[2];

    /*
     * Next step, compute STM derivatives
     */
    double ddots[6];    // {dxdx, dydy, dzdz, dxdy, dxdz, dydz}
    getUDDots(pSys->getMu(), s[0], s[1], s[2], ddots);

    /*
     *  PhiDot = A * Phi
     *  s[6] through s[42] represent the STM, Phi, in row-major order 
     *  sdot [6] through [42] is thus the derivative of the STM
     */
    std::copy(s+24, s+42, sdot+6); // First three rows of PhiDot are the last three rows of Phi
    for(int i = 0; i < 6; i++){
        sdot[24+i] = ddots[0]*s[6+i] + ddots[3]*s[12+i] + ddots[4]*s[18+i];
        sdot[30+i] = ddots[3]*s[6+i] + ddots[1]*s[12+i] + ddots[5]*s[18+i];
        sdot[36+i] = ddots[4]*s[6+i] + ddots[5]*s[12+i] + ddots[2]*s[18+i];
    }   // Last three rows are a combo of A and Phi

    return GSL_SUCCESS;
}//====================================================

/**
 *  @brief Integrate the equations of motion for the 2BP; currently the same as the full EOMs
 *  @param t time at integration step (unused)
 *  @param s the 6-d state vector
 *  @param sdot the 6-d state derivative vector
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_2bp::simpleEOMs(double t, const double s[], double sdot[], void *params){
	(void)t;
	
	EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_2bp *pSys = static_cast<const SysData_2bp *>(paramStruct->pSysData);

    double r = sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] );
    double mult = -pSys->getMu()/pow(r,3);	// G*(mass_primary)/(r^3)

    // Position derivatives = velocity
    std::copy(s+3, s+6, sdot);

    // Velocity derivatives = acceleration
    sdot[3] = mult*s[0];
    sdot[4] = mult*s[1];
    sdot[5] = mult*s[2];

    return GSL_SUCCESS;
}//====================================================

/**
 *  @brief Compute the second derivatives of the potential function
 *
 *  @param mu the mass parameter of the system, km^3/s^2
 *  @param x coordinate, km
 *  @param y coordinate, km
 *  @param z coordinate, km
 *  @param ddots a pointer to a 6-element double array where the function will store 
 *  values for {Uxx, Uyy, Uzz, Uxy, Uxz, Uyz}. Note that Uyx = Uxy, etc.
 */
void DynamicsModel_2bp::getUDDots(double mu, double x, double y, double z, double* ddots){
    // compute distance to primaries
    double r = sqrt(x*x + y*y + z*z);

    ddots[0] = 3*mu*x*x/pow(r,5) - mu/pow(r,3);	// Uxx
    ddots[1] = 3*mu*y*y/pow(r,5) - mu/pow(r,3);	// Uyy
    ddots[2] = 3*mu*z*z/pow(r,5) - mu/pow(r,3);	// Uzz
    ddots[3] = 3*mu*x*y/pow(r,5);	// Uxy
    ddots[4] = 3*mu*x*z/pow(r,5);	// Uxz
    ddots[5] = 3*mu*y*z/pow(r,5);	// Uyz
}//========================================================




}