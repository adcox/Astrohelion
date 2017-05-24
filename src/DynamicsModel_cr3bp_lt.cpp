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
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP_LT) {
    coreDim = 7;
    extraDim = 0;
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
 *  \param pSysData object describing the specific system
 *  \return an n x 3 vector (row-major order) containing the positions of
 *  n primaries; each row is one position vector in non-dimensional units
 */
std::vector<double> DynamicsModel_cr3bp_lt::getPrimPos(double t, const SysData *pSysData) const{
    std::vector<double> primPos(6,0);
    getPrimPos(t, pSysData, -1, &(primPos.front()));
    return primPos;
}//====================================================

/**
 *  \brief Compute the position of a specified primary
 *  \details This is the faster alternative to getPrimPos(t, pSysData).
 * 
 *  \param t Nondimensional time
 *  \param pSysData pointer to system data object
 *  \param pIx Index of the primary; a value of -1 will return the positions of all primaries,
 *  in order of largest to smallest mass
 *  \param pos An array to store the primary position(s) in with all elements initialized to zero.
 *  For a single primary position, the array must have at least three elements allocated. For all 
 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
 *  of primaries.
 */
void DynamicsModel_cr3bp_lt::getPrimPos(double t, const SysData *pSysData, int pIx, double *pos) const{
    (void) t;
    const SysData_cr3bp *pCrSys = static_cast<const SysData_cr3bp *>(pSysData);

    switch(pIx){
        case -1:
            pos[0] = -1*pCrSys->getMu();
            pos[3] = 1 - pCrSys->getMu();
            break;
        case 0:
            pos[0] = -1*pCrSys->getMu();
            break;
        case 1:
            pos[0] = 1 - pCrSys->getMu();
            break;
        default:
            throw Exception("DynamicsModel_cr3bp::getPrimPos: primary index out of bounds.");
    }
}//====================================================

/**
 *  \brief Compute the velocities of all primaries
 *
 *  \param t the epoch at which the computations occur (unused for this system)
 *  \param pSysData object describing the specific system (unused for this system)
 *  \return an n x 3 vector (row-major order) containing the velocities of
 *  n primaries; each row is one velocity vector in non-dimensional units
 */
std::vector<double> DynamicsModel_cr3bp_lt::getPrimVel(double t, const SysData *pSysData) const{
    std::vector<double> vel(6,0);
    getPrimVel(t, pSysData, -1, &(vel.front()));
    return vel;
}//====================================================

/**
 *  \brief Compute the velocity of a specified primary
 *  \details This is the faster alternative to getPrimVel(t, pSysData).
 * 
 *  \param t Nondimensional time
 *  \param pSysData pointer to system data object
 *  \param pIx Index of the primary; a value of -1 will return the velocities of all primaries,
 *  in order of largest to smallest mass
 *  \param vel An array to store the primary velocity(s) in with all elements initialized to zero. 
 *  For a single primary velocity, the array must have at least three elements allocated. For all 
 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
 *  of primaries.
 */
void DynamicsModel_cr3bp_lt::getPrimVel(double t, const SysData *pSysData, int pIx, double *vel) const{
    (void) t;
    (void) pSysData;
    (void) pIx;
    (void) vel;
}//====================================================

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
    if(state.size() != coreDim)
        throw Exception("DynamicsModel_cr3bp_lt::getStateDeriv: State size does not match the core state size specified by the dynamical model");

    // Compute the acceleration
    std::vector<double> dsdt(coreDim, 0);
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
void DynamicsModel_cr3bp_lt::multShoot_createOutput(const MultShootData *it) const{

    std::vector<int> newNodeIDs;
    newNodeIDs.reserve(it->numNodes);
    
    for(int n = 0; n < it->numNodes; n++){
        MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodesIn->getNodeByIx(n).getID());
        std::vector<double> state;

        if(state_var.row0 == -1){
            state = it->nodesIn->getState(state_var.key.id);
        }else{
            state = std::vector<double>(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0 + coreDim);
        }

        Node node(state, 0);
        node.setConstraints(it->nodesIn->getNodeByIx(n).getConstraints());

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(it->nodesOut->addNode(node));
    }

    double tof;
    int newOrigID, newTermID;
    for(unsigned int s = 0; s < it->nodesIn->getNumSegs(); s++){
        Segment seg = it->nodesIn->getSegByIx(s);

        if(it->bVarTime){
            MSVarMap_Obj tofVar = it->getVarMap_obj(it->bEqualArcTime ? MSVar_tp::TOF_TOTAL : MSVar_tp::TOF,
                it->bEqualArcTime ? Linkable::INVALID_ID : seg.getID());
            // Get data
            tof = it->bEqualArcTime ? it->X[tofVar.row0]/(it->nodesIn->getNumSegs()) : it->X[tofVar.row0];
        }else{
            tof = seg.getTOF();
        }

        newOrigID = newNodeIDs[it->nodesIn->getNodeIx(seg.getOrigin())];
        int termID = seg.getTerminus();
        newTermID = termID == Linkable::INVALID_ID ? termID : newNodeIDs[it->nodesIn->getNodeIx(termID)];
        
        Segment newSeg(newOrigID, newTermID, tof);
        newSeg.setConstraints(seg.getConstraints());
        newSeg.setVelCon(seg.getVelCon());
        newSeg.setSTM(it->propSegs[s].getSTMByIx(-1));
        newSeg.setCtrlLaw(seg.getCtrlLaw());
        newSeg.setStateVector(it->propSegs[s].getSegRef_const(0).getStateVector());
        newSeg.setTimeVector(it->propSegs[s].getSegRef_const(0).getTimeVector());
        it->nodesOut->addSeg(newSeg);
    }

    // Determine the chronological order of the nodeset
    // it->nodesOut->print();
    // it->nodesOut->printInChrono();
    std::vector<ArcPiece> order = it->nodesOut->getChronoOrder();
    // Set the epoch of each node based on the time of flight from
    // the first node
    double epoch = NAN;
    for(unsigned int i = 0; i < order.size(); i++){
        if(order[i].type == ArcPiece::Piece_tp::NODE){
            if(std::isnan(epoch)){
                // Copy the epoch value of the first node
                epoch = it->nodesOut->getNode(order[i].id).getEpoch();
            }else{
                // Set the epoch value of all other nodes
                it->nodesOut->getNodeRef(order[i].id).setEpoch(epoch);       
            }
        }
        if(order[i].type == ArcPiece::Piece_tp::SEG){
            if(!std::isnan(epoch)){
                // When stepping through in chronological order, every step is
                // forward in time; negative TOFs are associated with segments that
                // flow opposite the chronological order; ignore sign here.
                epoch += std::abs(it->nodesOut->getSeg(order[i].id).getTOF());
            }
        }
    }

    // it->nodesOut->print();
    std::vector<Constraint> arcCons = it->nodesIn->getArcConstraints();
    for(unsigned int i = 0; i < arcCons.size(); i++){
        it->nodesOut->addConstraint(arcCons[i]);
    }
}//======================================================

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

    double mu = sysData->getMu();           // nondimensional mass ratio
    double charT = sysData->getCharT();     // characteristic time (sec)
    double charL = sysData->getCharL();     // characteristic length (km)

    double T = law->getThrust();            // thrust magnitude (Newtons)
    double Isp = law->getIsp();             // specific impulse (sec)

    double f = (T/1000)*charT*charT/charL/sysData->getRefMass();    // nondimensional thrust

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    double thrust_dir[3];
    
    law->getLaw(t, s, sysData, thrust_dir, 3);

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
    double law_partials[21];
    law->getPartials_State(t, s, sysData, law_partials, 21);

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
 *  \param s the state vector passed in from the SimEngine. This vector includes
 *  the core states and control states, in that order.
 *  \param sdot the 7-d state derivative vector
 *  \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::simpleEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *sysData = static_cast<const SysData_cr3bp_lt *>(paramStruct->pSysData);
    const ControlLaw_cr3bp_lt *law = static_cast<const ControlLaw_cr3bp_lt *>(paramStruct->pCtrlLaw);

    double mu = sysData->getMu();           // nondimensional mass ratio
    double charT = sysData->getCharT();     // characteristic time (sec)
    double charL = sysData->getCharL();     // characteristic length (km)

    double T = law->getThrust();            // thrust magnitude (Newtons)
    double Isp = law->getIsp();             // specific impulse (sec)

    double f = (T/1000)*charT*charT/charL/sysData->getRefMass();    // nondimensional thrust

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    double thrust_dir[3];
    
    law->getLaw(t, s, sysData, thrust_dir, 3);

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