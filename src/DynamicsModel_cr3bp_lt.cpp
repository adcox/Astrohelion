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

#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "CorrectionEngine.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Nodeset_cr3bp_lt.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj_cr3bp_lt.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{

/**
 *  \brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic DynamicsModel
 */
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP_LT) {
    coreStates = 7;
    extraStates = 0;
    allowedCons.push_back(Constraint_tp::JC);
    allowedEvents.push_back(Event_tp::JC);
    allowedEvents.push_back(Event_tp::MASS);
}//==============================================

/**
 *  \brief Copy Constructor
 *  \param m a model reference
 */ 
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt &m) : DynamicsModel(m) {}

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
 *  \brief Compute the positions of all primaries
 *
 *  \param t the epoch at which the computations occur (unused for this system)
 *  \param sysData object describing the specific system
 *  \return an n x 3 vector (row-major order) containing the positions of
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
 *  \brief Compute the velocities of all primaries
 *
 *  \param t the epoch at which the computations occur (unused for this system)
 *  \param sysData object describing the specific system (unused for this system)
 *  \return an n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one velocity vector in non-dimensional units
 */
std::vector<double> DynamicsModel_cr3bp_lt::getPrimVel(double t, const SysData *sysData) const{
    (void)t;
    (void)sysData;
    double primVel[6] = {0};
    
    return std::vector<double>(primVel, primVel+6);
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
    if(state.size() != coreStates)
        throw Exception("DynamicsModel_cr3bp_lt::getStateDeriv: State size does not match the core state size specified by the dynamical model");

    // Compute the acceleration
    std::vector<double> dsdt(coreStates, 0);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), params);
    
    return dsdt;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Takes an input state and time and saves the data to the trajectory
 *  \param y an array containing the core state and any extra states integrated
 *  by the EOM function, including STM elements.
 *  \param t the time at the current integration state
 *  \param traj a pointer to the trajectory we should store the data in
 *  \param params structure containing parameters required by the EOMs
 */
void DynamicsModel_cr3bp_lt::sim_saveIntegratedData(const double* y, double t, Traj* traj, EOM_ParamStruct *params) const{
    
    DynamicsModel::sim_saveIntegratedData(y, t, traj, params);

    // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const SysData_cr3bp_lt *ltSys = static_cast<const SysData_cr3bp_lt*>(params->sysData);
    Traj_cr3bp_lt *ltTraj = static_cast<Traj_cr3bp_lt*>(traj);

    // Save Jacobi for CR3BP - it won't be constant any more, but is definitely useful to have
    ltTraj->setJacobiByIx(-1, DynamicsModel_cr3bp::getJacobi(y, ltSys->getMu()));
}//=====================================================

/**
 *  \brief Use a correction algorithm to accurately locate an event crossing
 *
 *  The simulation engine calls this function if and when it determines that an event 
 *  has been crossed. To accurately locate the event, we employ differential corrections
 *  and find the exact event occurence in space and time.
 *
 *  \param event the event we're looking for
 *  \param traj a pointer to the trajectory the event should occur on
 *  \param ic the core state vector for this system
 *  \param t0 non-dimensional time at the beginning of the search arc
 *  \param tof the time-of-flight for the arc to search over
 *  \param params structure containing parameters relevant to the integration
 *  \param verbose whether or not we should be verbose with output messages
 *
 *  \return wether or not the event has been located. If it has, a new point
 *  has been appended to the trajectory's data vectors.
 */
bool DynamicsModel_cr3bp_lt::sim_locateEvent(Event event, Traj* traj,
    const double *ic, double t0, double tof, EOM_ParamStruct *params, Verbosity_tp verbose) const{

    const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt *>(params->sysData);

    astrohelion::printVerb(verbose >= Verbosity_tp::ALL_MSG, "  Creating nodeset for event location\n");
    Nodeset_cr3bp_lt eventNodeset(pSys, ic, tof, 2, Nodeset::TIME, params->ctrlLawID);

    Constraint fixFirstCon(Constraint_tp::STATE, 0, ic, 7);
    Constraint eventCon(event.getConType(), 1, event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose == Verbosity_tp::ALL_MSG){ eventNodeset.print(); }

    astrohelion::printVerb(verbose >= Verbosity_tp::ALL_MSG, "  Applying corrections process to locate event\n");
    CorrectionEngine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbosity(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    Nodeset_cr3bp_lt correctedNodes(pSys);
    try{
        corrector.multShoot(&eventNodeset, &correctedNodes);
    }catch(DivergeException &e){
        if(verbose >= Verbosity_tp::SOME_MSG)
            astrohelion::printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(LinAlgException &e){
        if(verbose >= Verbosity_tp::SOME_MSG)
            astrohelion::printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    std::vector<double> state = correctedNodes.getStateByIx(-1);
    std::vector<double> stm = correctedNodes.getExtraParamVecByIx(-1, PARAMKEY_STM);
    state.insert(state.end(), stm.begin(), stm.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOFByIx(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    sim_saveIntegratedData(&(state[0]), eventTime, traj, params);

    return true;
}//=======================================================

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
 *  \brief Take the final, corrected free variable vector <tt>X</tt> and create an output 
 *  nodeset
 *
 *  If <tt>findEvent</tt> is set to true, the
 *  output nodeset will contain extra information for the simulation engine to use. Rather than
 *  returning only the position and velocity states, the output nodeset will contain the STM 
 *  and dqdT values for the final node; this information will be appended to the extraParameter
 *  vector in the final node.
 *
 *  \param it an iteration data object containing all info from the corrections process
 *  \param nodes_in a pointer to the original, uncorrected nodeset
 *  \param findEvent whether or not this correction process is locating an event
 *  \param nodes_out pointer to the nodeset object that will contain the output of the
 *  shooting process
 */
void DynamicsModel_cr3bp_lt::multShoot_createOutput(const MultShootData *it, const Nodeset *nodes_in, bool findEvent, Nodeset *nodes_out) const{
    
    const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt *>(it->sysData);
    Nodeset_cr3bp_lt *nodeset_out = static_cast<Nodeset_cr3bp_lt *>(nodes_out);

    std::vector<int> newNodeIDs;
    for(int n = 0; n < it->numNodes; n++){
        MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodeset->getNodeByIx(n).getID());
        std::vector<double> state(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0 + coreStates);

        Node node(state, 0);
        node.setConstraints(it->nodeset->getNodeByIx(n).getConstraints());

        if(n+1 == it->numNodes){
            // Set Jacobi Constant
            // node.setExtraParam("J", getJacobi(&(state[0]), pSys->getMu()));

            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
            state for the last node. We do this by appending the STM elements and dqdT elements to the
            end of the node array. This output nodeset should have two "nodes": the first 6 elements
            are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
            information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                Traj lastSeg = it->propSegs.back();
                MatrixXRd stm = lastSeg.getSTMByIx(-1);
                std::vector<double> stm_vec(stm.data(), stm.data() + stm.rows()*stm.cols());
                
                node.setExtraParamVec(PARAMKEY_STM, stm_vec);
            }
        }

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(nodeset_out->addNode(node));
        // nodeset_out->setJacobi(newNodeIDs.back(), getJacobi(&(state[0]), pSys->getMu()));
    }

    double tof;
    int newOrigID, newTermID;
    for(int s = 0; s < it->nodeset->getNumSegs(); s++){
        Segment seg = it->nodeset->getSegByIx(s);

        if(it->bVarTime){
            MSVarMap_Obj tofVar = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
                it->bEqualArcTime ? Linkable::INVALID_ID : seg.getID());
            // Get data
            tof = it->bEqualArcTime ? it->X[tofVar.row0]/(it->nodeset->getNumSegs()) : it->X[tofVar.row0];
        }else{
            tof = seg.getTOF();
        }

        newOrigID = newNodeIDs[it->nodeset->getNodeIx(seg.getOrigin())];
        int termID = seg.getTerminus();
        newTermID = termID == Linkable::INVALID_ID ? termID : newNodeIDs[it->nodeset->getNodeIx(termID)];
        
        Segment newSeg(newOrigID, newTermID, tof);
        newSeg.setConstraints(seg.getConstraints());
        newSeg.setVelCon(seg.getVelCon());
        newSeg.setSTM(it->propSegs[s].getSTMByIx(-1));
        newSeg.setCtrlLaw(seg.getCtrlLaw());
        nodeset_out->addSeg(newSeg);
    }

    // Determine the chronological order of the nodeset
    // nodeset_out->print();
    // nodeset_out->printInChrono();
    std::vector<ArcPiece> order = nodeset_out->getChronoOrder();
    // Set the epoch of each node based on the time of flight from
    // the first node
    double epoch = NAN;
    for(unsigned int i = 0; i < order.size(); i++){
        if(order[i].type == ArcPiece::Piece_tp::NODE){
            if(std::isnan(epoch)){
                // Copy the epoch value of the first node
                epoch = nodeset_out->getNode(order[i].id).getEpoch();
            }else{
                // Set the epoch value of all other nodes
                nodeset_out->getNodeRef(order[i].id).setEpoch(epoch);       
            }
        }
        if(order[i].type == ArcPiece::Piece_tp::SEG){
            if(!std::isnan(epoch)){
                // When stepping through in chronological order, every step is
                // forward in time; negative TOFs are associated with segments that
                // flow opposite the chronological order; ignore sign here.
                epoch += std::abs(nodeset_out->getSeg(order[i].id).getTOF());
            }
        }
    }

    // nodeset_out->print();
    std::vector<Constraint> arcCons = it->nodeset->getArcConstraints();
    for(unsigned int i = 0; i < arcCons.size(); i++){
        nodeset_out->addConstraint(arcCons[i]);
    }
}//======================================================

/**
 *  \brief Perform model-specific initializations on the MultShootData object
 *  \param it pointer to the object to be initialized
 */
void DynamicsModel_cr3bp_lt::multShoot_initIterData(MultShootData *it) const{
    Traj_cr3bp_lt traj(static_cast<const SysData_cr3bp_lt *>(it->sysData));
    it->propSegs.assign(it->nodeset->getNumSegs(), traj);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Integrate the equations of motion for the CR3BP LTVP
 *  \param t the current time of the integration
 *  \param s the 43-d state vector. The first 6 elements are position and velocity,
 *  the 7th is mass, and the final 36 are STM elements
 *  \param sdot the 43-d state derivative vector
 *  \param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::fullEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->sysData);

    double mu = sysData->getMu();           // nondimensional mass ratio
    double f = sysData->getThrust();        // nondimensional thrust magnitude
    double Isp = sysData->getIsp();         // specific impulse (sec)
    double charT = sysData->getCharT();     // characteristic time (sec)
    double charL = sysData->getCharL();     // characteristic length (km)

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    double thrust_dir[3];
    const ControlLaw_cr3bp_lt* law = static_cast<const ControlLaw_cr3bp_lt *>(sysData->getControlLaw());
    law->getLaw(t, s, sysData, paramStruct->ctrlLawID, thrust_dir, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + f/(s[6])*thrust_dir[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + f/(s[6])*thrust_dir[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + f/(s[6])*thrust_dir[2];
    sdot[6] = -charL*f/(charT*Isp*G_GRAV_0);   // nondimensional mass flow rate

    MatrixXRd A = MatrixXRd::Zero(7, 7);

    // Velocity relationships
    A(0, 3) = 1;
    A(1, 4) = 1;
    A(2, 5) = 1;

    // Uxx
    A(3, 0) = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((s[0] + mu),2)/pow(d,5) + 
        3*mu*pow((s[0] + mu - 1), 2)/pow(r,5);
    // Uxy
    A(3, 1) = 3*(1-mu)*(s[0] + mu)*s[1]/pow(d,5) + 3*mu*(s[0] + mu - 1)*s[1]/pow(r,5);
    // Uxz
    A(3, 2) = 3*(1-mu)*(s[0] + mu)*s[2]/pow(d,5) + 3*mu*(s[0] + mu - 1)*s[2]/pow(r,5);

    // Uyy
    A(4, 1) = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*s[1]*s[1]/pow(d,5) + 3*mu*s[1]*s[1]/pow(r,5);

    // Uyz
    A(4, 2) = 3*(1-mu)*s[1]*s[2]/pow(d,5) + 3*mu*s[1]*s[2]/pow(r,5);

    // Uzz
    A(5, 2) = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*s[2]*s[2]/pow(d,5) + 3*mu*s[2]*s[2]/pow(r,5);

    // Symmetry
    A(4,0) = A(3,1);
    A(5,0) = A(3,2);
    A(5,1) = A(4,2);

    A(3, 4) = 2;
    A(4, 3) = -2;

    // Get partials of control law and add them to the linear relationship matrix, A
    double law_partials[3*7];
    law->getPartials_State(t, s, sysData, paramStruct->ctrlLawID, law_partials, 21);

    for(unsigned int r = 0; r < 3; r++){
        for(unsigned int c = 0; c < 7; c++){
            A(3+r, c) += law_partials[r*7 + c];
        }
    }

    // toCSV(A, "LT_A.csv");
    // waitForUser();

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
 *  \brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  \param t time at integration step (unused)
 *  \param s the 7-d state vector
 *  \param sdot the 7-d state derivative vector
 *  \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::simpleEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->sysData);

    double mu = sysData->getMu();           // nondimensional mass ratio
    double f = sysData->getThrust();        // nondimensional thrust magnitude
    double Isp = sysData->getIsp();         // specific impulse (sec)
    double charT = sysData->getCharT();     // characteristic time (sec)
    double charL = sysData->getCharL();     // characteristic length (km)

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    double thrust_dir[3];
    const ControlLaw_cr3bp_lt* law = static_cast<const ControlLaw_cr3bp_lt *>(sysData->getControlLaw());
    law->getLaw(t, s, sysData, paramStruct->ctrlLawID, thrust_dir, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + f/(s[6])*thrust_dir[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + f/(s[6])*thrust_dir[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + f/(s[6])*thrust_dir[2];
    sdot[6] = -charL*f/(charT*Isp*G_GRAV_0);   // nondimensional mass flow rate (G_GRAV_0 is in km/s^2)

    return GSL_SUCCESS;
}//=====================================================

/**
 *  \brief Compute the Jacobi Constant for the CR3BP
 *
 *  \param s the state vector; only the position and velocity states are required
 *  \param mu the non-dimensional system mass ratio
 *
 *  \return the Jacobi Constant at this specific state and system
 */
double DynamicsModel_cr3bp_lt::getJacobi(const double s[], double mu){
    double v_squared = s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);
    return 2*U - v_squared;
}//================================================


}// END of Astrohelion namespace