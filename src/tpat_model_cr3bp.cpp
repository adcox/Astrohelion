/**
 *  @file tpat_model_cr3bp.cpp
 *  @brief Derivative of tpat_model, specific to CR3BP
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
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_traj_step.hpp"
#include "tpat_utilities.hpp"

/**
 *  @brief Construct a CR3BP Dynamic Model
 */
tpat_model_cr3bp::tpat_model_cr3bp() : tpat_model(MODEL_CR3BP) {
    // Allow a few more constraints than the default
    allowedCons.push_back(tpat_constraint::JC);
    allowedCons.push_back(tpat_constraint::PSEUDOARC);
    allowedEvents.push_back(tpat_event::JC);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */
tpat_model_cr3bp::tpat_model_cr3bp(const tpat_model_cr3bp &m) : tpat_model(m) {}

/**
 *  @brief Assignment operator
 *  @param m a model reference
 */
tpat_model_cr3bp& tpat_model_cr3bp::operator =(const tpat_model_cr3bp &m){
	tpat_model::operator =(m);
	return *this;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
tpat_model::eom_fcn tpat_model_cr3bp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  @brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
tpat_model::eom_fcn tpat_model_cr3bp::getFullEOM_fcn() const{
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
std::vector<double> tpat_model_cr3bp::getPrimPos(double t, const tpat_sys_data *sysData) const{
    (void)t;
    double primPos[6] = {0};
    const tpat_sys_data_cr3bp crSys(*static_cast<const tpat_sys_data_cr3bp *>(sysData));
    
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
std::vector<double> tpat_model_cr3bp::getPrimVel(double t, const tpat_sys_data *sysData) const{
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
void tpat_model_cr3bp::sim_saveIntegratedData(const double* y, double t, tpat_traj* traj) const{

	// Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
    const tpat_sys_data_cr3bp *crSys = static_cast<const tpat_sys_data_cr3bp*>(traj->getSysData());

    // Compute acceleration (elements 3 - 5)
    double dsdt[6] = {0};
    eomParamStruct paramStruct(crSys);
    simpleEOMs(t, y, dsdt, &paramStruct);

    // step(state, time, accel, stm) - y(0:5) holds the state, y(6:41) holds the STM
    tpat_traj_step step(y, t, dsdt+3, y+6);
    traj->appendStep(step);

    tpat_traj_cr3bp *cr3bpTraj = static_cast<tpat_traj_cr3bp*>(traj);    
    cr3bpTraj->setJacobi(-1, getJacobi(y, crSys->getMu()));
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
bool tpat_model_cr3bp::sim_locateEvent(tpat_event event, tpat_traj* traj,
    const double *ic, double t0, double tof, tpat_verbosity_tp verbose) const{

    const tpat_sys_data_cr3bp *crSys = static_cast<const tpat_sys_data_cr3bp*>(traj->getSysData());

    // Create a nodeset for this particular type of system
    printVerb(verbose == ALL_MSG, "  Creating nodeset for event location\n");
    tpat_nodeset_cr3bp eventNodeset(ic, crSys, tof, 2, tpat_nodeset::DISTRO_TIME);

    // Constraint to keep first node unchanged
    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, ic, 6);

    // Constraint to enforce event
    tpat_constraint eventCon(event.getConType(), 1, event.getConData());

    eventNodeset.addConstraint(fixFirstCon);
    eventNodeset.addConstraint(eventCon);

    if(verbose == ALL_MSG){ eventNodeset.print(); }

    printVerb(verbose == ALL_MSG, "  Applying corrections process to locate event\n");
    tpat_correction_engine corrector;
    corrector.setVarTime(true);
    corrector.setTol(traj->getTol());
    corrector.setVerbose(verbose);
    corrector.setFindEvent(true);   // apply special settings to minimize computations
    try{
        corrector.multShoot(&eventNodeset);
    }catch(tpat_diverge &e){
        printErr("Unable to locate event; corrector diverged\n");
        return false;
    }catch(tpat_linalg_err &e){
        printErr("LinAlg Err while locating event; bug in corrector!\n");
        return false;
    }

    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    tpat_nodeset_cr3bp correctedNodes = corrector.getCR3BP_Output();

    std::vector<double> state = correctedNodes.getNode(-1).getPosVelState();
    std::vector<double> extra = correctedNodes.getNode(-1).getExtraParams();
    extra.insert(extra.begin(), state.begin(), state.end());

    // event time is the TOF of corrected path + time at the state we integrated from
    double eventTime = correctedNodes.getTOF(0) + t0;

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
 *  This function calls its relative in the tpat_model base class and appends additional
 *  instructions specific to the CR3BP
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void tpat_model_cr3bp::multShoot_applyConstraint(iterationData *it, tpat_constraint con, int c) const{

    // Let the base class do its thing first
    tpat_model::multShoot_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case tpat_constraint::JC:
            multShoot_targetJC(it, con, row0);
            break;
        case tpat_constraint::PSEUDOARC:
            multShoot_targetPseudoArc(it, con, row0);
        default: break;
    }
}//=========================================================

/**
 *  @brief Compute constraint function and partial derivative values for a Jacobi Constraint
 *
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the row this constraint begins on
 */
void tpat_model_cr3bp::multShoot_targetJC(iterationData* it, tpat_constraint con, int row0) const{
    std::vector<double> conData = con.getData();
    int n = con.getNode();
    const tpat_sys_data_cr3bp *crSys = static_cast<const tpat_sys_data_cr3bp *> (it->sysData);

    // Compute the value of Jacobi at this node
    double mu = crSys->getMu();
    double nodeState[6];
    std::copy(&(it->X[6*n]), &(it->X[6*n])+6, nodeState);
    
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

    it->DF[it->totalFree*row0 + 6*n + 0] = (-2*(x + mu)*(1 - mu)/pow(d,3) - 2*(x + mu - 1)*mu/pow(r,3) + 2*x)/sr;    //dFdx
    it->DF[it->totalFree*row0 + 6*n + 1] = (-2*y*(1 - mu)/pow(d,3) - 2*y*mu/pow(r,3) + 2*y)/sr;                      //dFdy
    it->DF[it->totalFree*row0 + 6*n + 2] = (-2*z*(1 - mu)/pow(d,3) - 2*z*mu/pow(r,3))/sr;                            //dFdz
    it->DF[it->totalFree*row0 + 6*n + 3] = -2*vx/sv;   //dFdx_dot
    it->DF[it->totalFree*row0 + 6*n + 4] = -2*vy/sv;   //dFdy_dot
    it->DF[it->totalFree*row0 + 6*n + 5] = -2*vz/sv;   //dFdz_dot
}//=============================================

/**
 *  @brief Compute constraint function and partial derivative values for Pseudo Arc-Length
 *  
 *  @param it a pointer to the corrector's iteration data structure
 *  @param con the constraint being applied
 *  @param row0 the row this constraint begins on
 */
void tpat_model_cr3bp::multShoot_targetPseudoArc(iterationData *it, tpat_constraint con, int row0) const{
    std::vector<double> conData = con.getData();

    if(row0 != it->totalCons-1)
        throw tpat_exception("tpat_model_cr3bp::multShoot_targetPseudoArc: Pseudo Arc-Length constraint must be the final constraint; please re-create the nodeset accordingly");

    if(it->totalCons != it->totalFree)
        throw tpat_exception("tpat_model_cr3bp::multShoot_targetPseudoArc: Jacobian matrix is not square; cannot apply pseudo arc-length");

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
 *
 *  @return a pointer to a nodeset containing the corrected nodes
 */
tpat_nodeset* tpat_model_cr3bp::multShoot_createOutput(const iterationData *it, const tpat_nodeset *nodes_in, bool findEvent) const{

    // Create a nodeset with the same system data as the input
    const tpat_sys_data_cr3bp *crSys = static_cast<const tpat_sys_data_cr3bp *>(it->sysData);
    tpat_nodeset_cr3bp *nodeset_out = new tpat_nodeset_cr3bp(crSys);

    int numNodes = (int)(it->origNodes.size());
    for(int i = 0; i < numNodes; i++){
        double state[6];
        std::copy(&(it->X[i*6]), &(it->X[i*6])+6, state);

        // Reverse scaling
        for(int i = 0; i < 6; i++){
            state[i] /= i < 3 ? it->freeVarScale[0] : it->freeVarScale[1];
        }

        tpat_node node(state, NAN);
        node.setVelCon(nodes_in->getNode(i).getVelCon());
        node.setConstraints(nodes_in->getNode(i).getConstraints());

        if(i + 1 < numNodes){
            // Get TOF, reverse variable scaling, save to node
            double tof;
            if(it->varTime){
                // Get data
                tof = it->equalArcTime ? it->X[6*it->numNodes]/(it->numNodes - 1) : it->X[6*it->numNodes+i];
                // Reverse scaling
                tof /= it->freeVarScale[2];    // Time scaling
            }else{
                tof = nodes_in->getTOF(i);
            }
            node.setTOF(tof);

            // Set Jacobi Constant
            node.setExtraParam(1, getJacobi(state, crSys->getMu()));
        }else{
            node.setTOF(NAN);
            // Set Jacobi Constant
            node.setExtraParam(1, getJacobi(state, crSys->getMu()));

            /* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
            state for the last node. We do this by appending the STM elements and dqdT elements to the
            end of the node array. This output nodeset should have two "nodes": the first 6 elements
            are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
            information*/
            if(findEvent){
                // Append the 36 STM elements to the node vector
                tpat_traj lastSeg = it->allSegs.back();
                MatrixXRd stm = lastSeg.getSTM(-1);
                std::vector<double> extraParam(stm.data(), stm.data()+36);
                
                node.setExtraParams(extraParam);
            }
        }

        nodeset_out->appendNode(node);
    }

    return nodeset_out;
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
int tpat_model_cr3bp::fullEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;

    // Extract mu from params
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const tpat_sys_data_cr3bp *sysData = static_cast<const tpat_sys_data_cr3bp *>(paramStruct->sysData);
    
    double mu = sysData->getMu();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4];

    // compute distance to primaries
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] =   2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3);
    sdot[4] =  -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3);
    sdot[5] =  -(1-mu)*z/pow(d,3) - mu*z/pow(r,3); 

    /*
     * Next step, compute STM
     */
    
    // Create A Matrix
    double ddots[6];
    getUDDots(mu, x, y, z, ddots);

    double a_data[] = { 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1,
                        ddots[0], ddots[3], ddots[4], 0, 2, 0,
                        ddots[3], ddots[1], ddots[5], -2, 0, 0,
                        ddots[4], ddots[5], ddots[2], 0, 0, 0};
    MatrixXRd A = Eigen::Map<MatrixXRd>(a_data, 6, 6);

    // Copy the STM states into a sub-array
    double stmElements[36];
    std::copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    MatrixXRd phi = Eigen::Map<MatrixXRd>(stmElements, 6, 6);
    
    // Compute derivative of STM
    MatrixXRd phiDot = A*phi;

    // Copy the elements of phiDot into the derivative array
    double *phiDotData = phiDot.data();
    std::copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 6-d state vector
 *  @param sdot the 6-d state derivative vector
 *  @param params points to an eomParamStruct object
 */
int tpat_model_cr3bp::simpleEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;
    
    // Extract mu from params
    // tpat_sys_data_cr3bp *sysData = static_cast<tpat_sys_data_cr3bp *>(params);
    eomParamStruct *paramStruct = static_cast<eomParamStruct *>(params);
    const tpat_sys_data_cr3bp *sysData = static_cast<const tpat_sys_data_cr3bp *>(paramStruct->sysData);
    double mu = sysData->getMu();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4];

    // compute distance to primaries
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] =   2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3);
    sdot[4] =  -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3);
    sdot[5] =  -(1-mu)*z/pow(d,3) - mu*z/pow(r,3);

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
 */
void tpat_model_cr3bp::getEquilibPt(const tpat_sys_data_cr3bp *sysData, int L, double tol, double pos[3]){
    if(L < 1 || L > 5){
        throw tpat_exception("Invalid Lagrange Point");
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
        printErr("Could not converge on L%d\n", L);
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
double tpat_model_cr3bp::getJacobi(const double s[], double mu){
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
void tpat_model_cr3bp::getUDDots(double mu, double x, double y, double z, double* ddots){
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

