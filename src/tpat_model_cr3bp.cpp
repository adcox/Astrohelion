/**
 *  @file tpat_model_cr3bp.cpp
 *  @brief Derivative of TPAT_Model, specific to CR3BP
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

#include "tpat_model_cr3bp.hpp"

#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_event.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <gsl/gsl_errno.h>
/**
 *  @brief Construct a CR3BP Dynamic Model
 */
TPAT_Model_CR3BP::TPAT_Model_CR3BP() : TPAT_Model(MODEL_CR3BP) {
    // Allow a few more constraints than the default
    allowedCons.push_back(TPAT_Constraint_Tp::JC);
    allowedCons.push_back(TPAT_Constraint_Tp::PSEUDOARC);
    allowedEvents.push_back(TPAT_Event_Tp::JC);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */
TPAT_Model_CR3BP::TPAT_Model_CR3BP(const TPAT_Model_CR3BP &m) : TPAT_Model(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
TPAT_Model_CR3BP& TPAT_Model_CR3BP::operator =(const TPAT_Model_CR3BP &m){
	TPAT_Model::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
TPAT_Model::eom_fcn TPAT_Model_CR3BP::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
TPAT_Model::eom_fcn TPAT_Model_CR3BP::getFullEOM_fcn() const{
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
std::vector<double> TPAT_Model_CR3BP::getPrimPos(double t, const TPAT_Sys_Data *sysData) const{
    (void)t;
    double primPos[6] = {0};
    const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP *>(sysData);
    
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
std::vector<double> TPAT_Model_CR3BP::getPrimVel(double t, const TPAT_Sys_Data *sysData) const{
    (void)t;
    (void)sysData;
    
    return std::vector<double>(6, 0);
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
void TPAT_Model_CR3BP::sim_saveIntegratedData(const double* y, double t, TPAT_Traj* traj) const{

	// Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP*>(traj->getSysData());

    // Compute acceleration (elements 3-5)
    double dsdt[6] = {0};
    eomParamStruct paramStruct(crSys);
    simpleEOMs(t, y, dsdt, &paramStruct);
    
    // node(state, accel, epoch) - y(0:5) holds the state, y(6:41) holds the STM
    int id = traj->addNode(TPAT_Node(y, dsdt+3, t));

    if(id > 0){
        double tof = t - traj->getNode(id-1).getEpoch();
        traj->addSeg(TPAT_Segment(id-1, id, tof, y+6));
    }

    TPAT_Traj_CR3BP *cr3bpTraj = static_cast<TPAT_Traj_CR3BP*>(traj);    
    cr3bpTraj->setJacobiByIx(-1, getJacobi(y, crSys->getMu()));
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
bool TPAT_Model_CR3BP::sim_locateEvent(TPAT_Event event, TPAT_Traj* traj,
    const double *ic, double t0, double tof, TPAT_Verbosity_Tp verbose) const{

    const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP*>(traj->getSysData());

    // Create a nodeset for this particular type of system
    printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  Creating nodeset for event location\n");
    TPAT_Nodeset_CR3BP eventNodeset(ic, crSys, tof, 2, TPAT_Nodeset::DISTRO_TIME);

    // Constraint to keep first node unchanged
    TPAT_Constraint fixFirstCon(TPAT_Constraint_Tp::STATE, 0, ic, 6);

    // Constraint to enforce event
    TPAT_Constraint eventCon(event.getConType(), 1, event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose == TPAT_Verbosity_Tp::ALL_MSG){ eventNodeset.print(); }

    printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  Applying corrections process to locate event\n");
    TPAT_Correction_Engine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbose(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    TPAT_Nodeset_CR3BP correctedNodes(crSys);
    try{
        corrector.multShoot(&eventNodeset, &correctedNodes);
    }catch(TPAT_Diverge &e){
        printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(TPAT_LinAlg_Err &e){
        printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    std::vector<double> state = correctedNodes.getNodeByIx(-1).getState();
    std::vector<double> extra = correctedNodes.getNodeByIx(-1).getExtraParams();
    extra.insert(extra.begin(), state.begin(), state.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOFByIx(0) + t0;

    // Use the data stored in nodes and save the state and time of the event occurence
    sim_saveIntegratedData(&(extra[0]), eventTime, traj);

    return true;
}//======================================================

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Compute constraint function and partial derivative values for a constraint
 *  
 *  This function calls its relative in the TPAT_Model base class and appends additional
 *  instructions specific to the CR3BP
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void TPAT_Model_CR3BP::multShoot_applyConstraint(TPAT_MultShoot_Data *it, TPAT_Constraint con, int c) const{

    // Let the base class do its thing first
    TPAT_Model::multShoot_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case TPAT_Constraint_Tp::JC:
            multShoot_targetJC(it, con, row0);
            break;
        case TPAT_Constraint_Tp::PSEUDOARC:
            multShoot_targetPseudoArc(it, con, row0);
        default: break;
    }
}//=========================================================

/**
 *  @brief Perform model-specific initializations on the TPAT_MultShoot_Data object
 *  @param it pointer to the object to be initialized
 */
void TPAT_Model_CR3BP::multShoot_initIterData(TPAT_MultShoot_Data *it) const{
    TPAT_Traj_CR3BP traj(static_cast<const TPAT_Sys_Data_CR3BP *>(it->sysData));
    it->propSegs.assign(it->nodeset->getNumSegs(), traj);
}//====================================================

/**
 *  @brief Compute constraint function and partial derivative values for a Jacobi Constraint
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the row this constraint begins on
 */
void TPAT_Model_CR3BP::multShoot_targetJC(TPAT_MultShoot_Data* it, TPAT_Constraint con, int row0) const{
    std::vector<double> conData = con.getData();
    MSVarMap_Obj state_var = it->getVarMap_obj(MSVarType::STATE, con.getID());
    // int nodeIx = it->nodeset->getNodeIx(con.getID());
    const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP *> (it->sysData);

    // Compute the value of Jacobi at this node
    double mu = crSys->getMu();
    double nodeState[6];
    std::copy(&(it->X[state_var.row0]), &(it->X[state_var.row0])+6, nodeState);
    
    double sr = it->freeVarScale[0];
    double sv = it->freeVarScale[1];    

    // Reverse scaling to compute Jacobi at the node
    for(int i = 0; i < 6; i++){
        double scale = i < 3 ? sr : sv;
        nodeState[i] /= scale;
    }

    // printf("Node State = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n", nodeState[0],
    //     nodeState[1], nodeState[2], nodeState[3], nodeState[4], nodeState[5]);
    double nodeJC = getJacobi(nodeState, mu);
    
    // temp variables to make equations more readable; compute partials w.r.t. node state
    double x = nodeState[0];
    double y = nodeState[1];
    double z = nodeState[2];
    double vx = nodeState[3];
    double vy = nodeState[4];
    double vz = nodeState[5];

    double d = sqrt((x + mu)*(x + mu) + y*y + z*z);
    double r = sqrt((x + mu - 1)*(x + mu - 1) + y*y + z*z);

    it->FX[row0] = nodeJC - conData[0];
    // printf("Targeting JC = %.4f, value is %.4f\n", conData[0], nodeJC);

    it->DF[it->totalFree*row0 + state_var.row0 + 0] = (-2*(x + mu)*(1 - mu)/pow(d,3) - 2*(x + mu - 1)*mu/pow(r,3) + 2*x)/sr;    //dFdx
    it->DF[it->totalFree*row0 + state_var.row0 + 1] = (-2*y*(1 - mu)/pow(d,3) - 2*y*mu/pow(r,3) + 2*y)/sr;                      //dFdy
    it->DF[it->totalFree*row0 + state_var.row0 + 2] = (-2*z*(1 - mu)/pow(d,3) - 2*z*mu/pow(r,3))/sr;                            //dFdz
    it->DF[it->totalFree*row0 + state_var.row0 + 3] = -2*vx/sv;   //dFdx_dot
    it->DF[it->totalFree*row0 + state_var.row0 + 4] = -2*vy/sv;   //dFdy_dot
    it->DF[it->totalFree*row0 + state_var.row0 + 5] = -2*vz/sv;   //dFdz_dot
}//=============================================

/**
 *  @brief Compute constraint function and partial derivative values for Pseudo Arc-Length
 *  
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the row this constraint begins on
 *  @throw TPAT_Exception if the pseudo arclength constraint is not listed as the final constraint
 *  @throw TPAT_Exception if the Jacobian matrix (w/o the PAL constraint) is nonsquare.
 */
void TPAT_Model_CR3BP::multShoot_targetPseudoArc(TPAT_MultShoot_Data *it, TPAT_Constraint con, int row0) const{
    std::vector<double> conData = con.getData();

    if(row0 != it->totalCons-1)
        throw TPAT_Exception("TPAT_Model_CR3BP::multShoot_targetPseudoArc: Pseudo Arc-Length constraint must be the final constraint; please re-create the nodeset accordingly");

    if(it->totalCons != it->totalFree)
        throw TPAT_Exception("TPAT_Model_CR3BP::multShoot_targetPseudoArc: Jacobian matrix is not square; cannot apply pseudo arc-length");

    // All elements except the last are the free-variable vector for a converged family member
    std::vector<double> famFreeVec(conData.begin(), conData.begin()+it->totalFree);
    std::vector<double> nullspace(conData.begin()+it->totalFree, conData.end()-1);
    double stepSize = conData.back();   // The last element is the step size

    Eigen::RowVectorXd X = Eigen::Map<Eigen::RowVectorXd>(&(it->X[0]), 1, it->totalFree);
    Eigen::RowVectorXd X_fam = Eigen::Map<Eigen::RowVectorXd>(&(famFreeVec[0]), 1, famFreeVec.size());
    Eigen::VectorXd N = Eigen::Map<Eigen::VectorXd>(&(nullspace[0]), nullspace.size(), 1);
    
    MatrixXRd dotProd;
    dotProd.noalias() = (X - X_fam)*N;

    it->FX[row0] = dotProd(0) - stepSize;

    for(int i = 0; i < it->totalFree; i++){
        it->DF[it->totalFree*row0 + i] = N(i);
    }
}//=============================================

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
 *  @return a pointer to a nodeset containing the corrected nodes
 */
void TPAT_Model_CR3BP::multShoot_createOutput(const TPAT_MultShoot_Data *it, const TPAT_Nodeset *nodes_in, bool findEvent, TPAT_Nodeset *nodes_out) const{

    // Create a nodeset with the same system data as the input
    const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP *>(it->sysData);
    TPAT_Nodeset_CR3BP *nodeset_out = static_cast<TPAT_Nodeset_CR3BP *>(nodes_out);

    std::vector<int> newNodeIDs;
    for(int n = 0; n < it->numNodes; n++){
        MSVarMap_Obj state_var = it->getVarMap_obj(MSVarType::STATE, it->nodeset->getNodeByIx(n).getID());
        double state[6];
        std::copy(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0+6, state);

        // Reverse scaling
        for(int i = 0; i < 6; i++){
            state[i] /= i < 3 ? it->freeVarScale[0] : it->freeVarScale[1];
        }

        TPAT_Node node(state, 0);
        node.setConstraints(nodes_in->getNodeByIx(n).getConstraints());

        if(n+1 == it->numNodes){
            // Set Jacobi Constant
            node.setExtraParam(0, getJacobi(state, crSys->getMu()));

            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
            state for the last node. We do this by appending the STM elements and dqdT elements to the
            end of the node array. This output nodeset should have two "nodes": the first 6 elements
            are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
            information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                TPAT_Traj lastSeg = it->propSegs.back();
                MatrixXRd stm = lastSeg.getSTMByIx(-1);
                std::vector<double> extraParam(stm.data(), stm.data()+36);
                
                node.setExtraParams(extraParam);
            }
        }

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(nodeset_out->addNode(node));
        nodeset_out->setJacobi(newNodeIDs.back(), getJacobi(state, crSys->getMu()));
    }

    double tof;
    int newOrigID, newTermID;
    for(int s = 0; s < it->nodeset->getNumSegs(); s++){
        TPAT_Segment seg = it->nodeset->getSegByIx(s);

        if(it->varTime){
            MSVarMap_Obj tofVar = it->getVarMap_obj(it->equalArcTime ? MSVarType::TOF_TOTAL : MSVarType::TOF,
                it->equalArcTime ? TPAT_Linkable::INVALID_ID : seg.getID());
            // Get data
            tof = it->equalArcTime ? it->X[tofVar.row0]/(it->nodeset->getNumSegs()) : it->X[tofVar.row0];
            // Reverse scaling
            tof /= it->freeVarScale[2];     // TOF scaling
        }else{
            tof = seg.getTOF();
        }

        newOrigID = newNodeIDs[it->nodeset->getNodeIx(seg.getOrigin())];
        int termID = seg.getTerminus();
        newTermID = termID == TPAT_Linkable::INVALID_ID ? termID : newNodeIDs[it->nodeset->getNodeIx(termID)];
        
        TPAT_Segment newSeg(newOrigID, newTermID, tof);
        newSeg.setConstraints(seg.getConstraints());
        newSeg.setVelCon(seg.getVelCon());
        nodeset_out->addSeg(newSeg);
    }

    std::vector<TPAT_Constraint> arcCons = nodes_in->getArcConstraints();
    for(size_t i = 0; i < arcCons.size(); i++){
        nodeset_out->addConstraint(arcCons[i]);
    }
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the CR3BP
 *  @param t the current time of the integration; not used for this system
 *  @param s the 42-d state vector
 *  @param sdot the 42-d state derivative vector
 *  @param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an eomParamStruct object
 */
int TPAT_Model_CR3BP::fullEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;

    // Extract mu from params
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const TPAT_Sys_Data_CR3BP *sysData = static_cast<const TPAT_Sys_Data_CR3BP *>(paramStruct->sysData);
    
    double mu = sysData->getMu();

    // double x = s[0];    double y = s[1];    double z = s[2];
    // double xdot = s[3]; double ydot = s[4];

    // compute distance to primaries
    double d = sqrt( (s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0] - 1+mu)*(s[0] - 1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Position derivatives = velocity
    std::copy(s+3, s+6, sdot);

    // Velocity derivatives = acceleraiton
    sdot[3] =   2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3);
    sdot[4] =  -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3);
    sdot[5] =  -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3); 

    /*
     * Next step, compute STM
     */
    double ddots[6];    // {dxdx, dydy, dzdz, dxdy, dxdz, dydz}
    getUDDots(mu, s[0], s[1], s[2], ddots);

    /*  Compute the STM Derivative 
     *  PhiDot = A * Phi
     *  s[6] through s[42] represent the STM, Phi, in row-major order 
     *  sdot [6] through [42] is thus the derivative of the STM
     */
    std::copy(s+24, s+42, sdot+6); // First three rows are the last three rows of Phi
    for(int i = 0; i < 6; i++){
        sdot[24+i] = ddots[0]*s[6+i] + ddots[3]*s[12+i] + ddots[4]*s[18+i] + 2*s[30+i];
        sdot[30+i] = ddots[3]*s[6+i] + ddots[1]*s[12+i] + ddots[5]*s[18+i] - 2*s[24+i];
        sdot[36+i] = ddots[4]*s[6+i] + ddots[5]*s[12+i] + ddots[2]*s[18+i];
    }   // Last three rows are a combo of A and Phi

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 6-d state vector
 *  @param sdot the 6-d state derivative vector
 *  @param params points to an eomParamStruct object
 */
int TPAT_Model_CR3BP::simpleEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;
    
    // Extract mu from params
    // TPAT_Sys_Data_CR3BP *sysData = static_cast<TPAT_Sys_Data_CR3BP *>(params);
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const TPAT_Sys_Data_CR3BP *sysData = static_cast<const TPAT_Sys_Data_CR3BP *>(paramStruct->sysData);
    double mu = sysData->getMu();

    // double x = s[0];    double y = s[1];    double z = s[2];
    // double xdot = s[3]; double ydot = s[4];

    // compute distance to primaries
    double d = sqrt( (s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0] - 1+mu)*(s[0] - 1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Position derivatives = velocity
    std::copy(s+3, s+6, sdot);

    // Velocity derivatives = acceleraiton
    sdot[3] =   2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3);
    sdot[4] =  -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3);
    sdot[5] =  -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3); 

    return GSL_SUCCESS;
}//=====================================================

/**
 *  @brief Compute the location of a Lagrange point in the CR3BP
 *
 *  @param sysData pointer to an object describing the particular CR3BP
 *  @param L the Lagrange point number, 1 to 5
 *  @param tol the tolerance to use; if NAN is input, then a default value of 1e-14 will
 *  be used.
 *  @param pos a 3-element array to store the position of the Lagrange point
 *  @throws TPAT_Diverge if the Newton-Raphson process fails to converge on the 
 *  Lagrange point location
 */
void TPAT_Model_CR3BP::getEquilibPt(const TPAT_Sys_Data_CR3BP *sysData, int L, double tol, double pos[3]){
    if(L < 1 || L > 5){
        throw TPAT_Exception("Invalid Lagrange Point");
    }

    if(tol == NAN){
        tol = 1e-14;
    }

    double mu = sysData->getMu();
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;

    double gamma;
    double gamma_prev = -999;
    int count = 0;
    const int maxCount = 20;

    switch(L){
        case 1:
            gamma = 0.1;    // Initial guess is 10% of orbital radius
            while(std::abs(gamma - gamma_prev) > tol && count < maxCount){   // Newton-Raphson for L1
                gamma_prev = gamma;
                gamma = gamma - ( mu/(gamma*gamma) - (1-mu)/pow(1-gamma, 2) - gamma - mu + 1)/
                    ( -2*mu/pow(gamma,3) - 2*(1-mu)/pow(1-gamma,3) - 1 );
                count++;
            }
            pos[0] = 1 - mu - gamma;
            break;
        case 2:
            gamma = 0.1;    // Initial guess is 10% of orbital radius
            while(std::abs(gamma - gamma_prev) > tol && count < maxCount){
                gamma_prev = gamma;
                gamma = gamma - ( -1*mu/(gamma*gamma) - (1-mu)/pow(1+gamma, 2) - mu + 1 + gamma)/
                    ( 2*mu/pow(gamma, 3) + 2*(1-mu)/pow(1+gamma, 3) + 1 );
                count++;
            }
            pos[0] = 1 - mu + gamma;
            break;
        case 3:
            gamma = 1;  // Initial guess is 100% of orbital radius
            while(std::abs(gamma - gamma_prev) > tol && count < maxCount){
                gamma_prev = gamma;
                gamma = gamma - ( mu/pow(-1 - gamma, 2) + (1-mu)/(gamma*gamma) - mu - gamma)/
                    ( -2*mu/pow(1+gamma,3) - 2*(1-mu)/pow(gamma, 3) - 1);
                count++;
            }
            pos[0] = -1*mu - gamma;
            break;
        case 4:
        case 5:
            pos[0] = 0.5 - mu;
            pos[1] = L == 4 ? sin(PI/3) : -1*sin(PI/3);
            break;
    }

    if(L < 4 && std::abs(gamma - gamma_prev) > tol){
        throw TPAT_Diverge("TPAT_Model_CR3BP::getEquilibPt: Could not converge on Lagrange point.");
    }
}//========================================

/**
 *  @brief Compute the Jacobi Constant for the CR3BP
 *
 *  @param s the state vector; only the position and velocity states are required
 *  @param mu the non-dimensional system mass ratio
 *
 *  @return the Jacobi Constant at this specific state and system
 */
double TPAT_Model_CR3BP::getJacobi(const double s[], double mu){
    double v_squared = s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);
    return 2*U - v_squared;
}//================================================

/**
 *  @brief Compute the second derivatives of the pseudo-potential function
 *
 *  @param mu the mass ratio of the system, non-dimensional
 *  @param x coordinate, non-dimensional units 
 *  @param y coordinate, non-dimensional units 
 *  @param z coordinate, non-dimensional units 
 *  @param ddots a pointer to a 6-element double array where the function will store 
 *  values for {Uxx, Uyy, Uzz, Uxy, Uxz, Uyz}. Note that Uyx = Uxy, etc.
 */
void TPAT_Model_CR3BP::getUDDots(double mu, double x, double y, double z, double* ddots){
    // compute distance to primaries
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );

    // Uxx
    ddots[0] = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((x + mu),2)/pow(d,5) + 
        3*mu*pow((x + mu - 1), 2)/pow(r,5);
    // Uyy
    ddots[1] = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*y*y/pow(d,5) + 3*mu*y*y/pow(r,5);
    // Uzz
    ddots[2] = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*z*z/pow(d,5) + 3*mu*z*z/pow(r,5);

    // Uxy
    ddots[3] = 3*(1-mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5);
    // Uxz
    ddots[4] = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5);
    // Uyz
    ddots[5] = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5);
}//========================================================

