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

#include "DynamicsModel_2bp.hpp"

#include "Calculations.hpp"
// #include "CorrectionEngine.hpp"
#include "EigenDefs.hpp"
#include "Event.hpp"
// #include "MultShootData.hpp"
#include "Node.hpp"
// #include "Nodeset_2bp.hpp"
#include "Segment.hpp"
#include "SysData_2bp.hpp"
#include "Traj_2bp.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{

DynamicsModel_2bp::DynamicsModel_2bp() : DynamicsModel(DynamicsModel_tp::MODEL_2BP) {
	// stmStates = 0;
}//====================================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */
DynamicsModel_2bp::DynamicsModel_2bp(const DynamicsModel_2bp &m) : DynamicsModel(m) {
	// stmStates = 0;
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

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

void DynamicsModel_2bp::sim_saveIntegratedData(const double *y, double t, Traj* traj) const{
	const SysData_2bp *sys = static_cast<const SysData_2bp *>(traj->getSysData());

	// Compute acceleration
	double dsdt[6] = {0};
	EOM_ParamStruct paramStruct(sys);
	simpleEOMs(t, y, dsdt, &paramStruct);

	// Save the state (pos, vel, accel)
	int id = traj->addNode(Node(y, dsdt+3, t));

	if(id > 0){
		double tof = t - traj->getNode(id-1).getEpoch();
		traj->addSeg(Segment(id-1, id, tof, y+6));
	}
}//====================================================

bool DynamicsModel_2bp::sim_locateEvent(Event event, Traj *traj, const double *ic, double t0, double tof, Verbosity_tp verbose) const{
	
	// const SysData_2bp *sys = static_cast<const SysData_2bp*>(traj->getSysData());

 //    // Create a nodeset for this particular type of system
 //    astrohelion::printVerb(verbose == Verbosity_tp::ALL_MSG, "  Creating nodeset for event location\n");
 //    Nodeset_2bp eventNodeset(ic, sys, tof, 2, Nodeset::DISTRO_TIME);

 //    // Constraint to keep first node unchanged
 //    Constraint fixFirstCon(Constraint_tp::STATE, 0, ic, 6);

 //    // Constraint to enforce event
 //    Constraint eventCon(event.getConType(), 1, event.getConData());

 //    eventNodeset.addConstraint(fixFirstCon);
 //    eventNodeset.addConstraint(eventCon);

 //    if(verbose == Verbosity_tp::ALL_MSG){ eventNodeset.print(); }

 //    astrohelion::printVerb(verbose == Verbosity_tp::ALL_MSG, "  Applying corrections process to locate event\n");
 //    CorrectionEngine corrector;
 //    corrector.setVarTime(true);
 //    corrector.setTol(traj->getTol());
 //    corrector.setVerbose(verbose);
 //    corrector.setFindEvent(true);   // apply special settings to minimize computations

 //    // Because we set findEvent to true, this output nodeset should contain
 //    // the full (42 or 48 element) final state
 //    Nodeset_2bp correctedNodes(sys);
 //    try{
 //        corrector.multShoot(&eventNodeset, &correctedNodes);
 //    }catch(DivergeException &e){
 //        astrohelion::printErr("Unable to locate event; corrector diverged\n");
 //        return false;
 //    }catch(LinAlgException &e){
 //        astrohelion::printErr("LinAlg Err while locating event; bug in corrector!\n");
 //        return false;
 //    }

 //    std::vector<double> state = correctedNodes.getNodeByIx(-1).getState();
 //    std::vector<double> extra = correctedNodes.getNodeByIx(-1).getExtraParams();
 //    extra.insert(extra.begin(), state.begin(), state.end());

 //    // event time is the TOF of corrected path + time at the state we integrated from
 //    double eventTime = correctedNodes.getTOFByIx(0) + t0;

 //    // Use the data stored in nodes and save the state and time of the event occurence
 //    sim_saveIntegratedData(&(extra[0]), eventTime, traj);

    return true;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

void DynamicsModel_2bp::multShoot_initIterData(MultShootData *it) const{
	(void)it;
}//====================================================

void DynamicsModel_2bp::multShoot_createOutput(const MultShootData* it, const Nodeset *nodes_in, bool findEvent, Nodeset *nodesOut) const{
	(void)it;
	(void)nodes_in;
	(void)findEvent;
	(void)nodesOut;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the 2BP
 *  @param t time at integration step (unused)
 *  @param s the 6-d state vector
 *  @param sdot the 6-d state derivative vector
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_2bp::fullEOMs(double t, const double s[], double sdot[], void *params){
	(void)t;
	
	EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_2bp *pSys = static_cast<const SysData_2bp *>(paramStruct->sysData);

    double r = sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2] );
    double mult = -pSys->getMu()/pow(r,3);	// G*(mass_primary)/(r^3)

    // Position derivatives = velocity
    std::copy(s+3, s+6, sdot);

    // Velocity derivatives = acceleration
    sdot[3] = mult*s[0];
    sdot[4] = mult*s[1];
    sdot[5] = mult*s[2];

    //STM Derivatives are 0 for now...
    for(int i = 6; i < 42; i++){
    	sdot[i] = 0;
    }

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
    const SysData_2bp *pSys = static_cast<const SysData_2bp *>(paramStruct->sysData);

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






}