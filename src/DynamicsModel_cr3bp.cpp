/**
 *  \file DynamicsModel_cr3bp.cpp
 *  \brief Derivative of DynamicsModel, specific to CR3BP
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

#include "DynamicsModel_cr3bp.hpp"

#include "Arcset_cr3bp.hpp"
#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "EigenDefs.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "MultShootEngine.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

#include <gsl/gsl_errno.h>

namespace astrohelion{
/**
 *  \brief Construct a CR3BP Dynamic DynamicsModel
 */
DynamicsModel_cr3bp::DynamicsModel_cr3bp() : DynamicsModel(DynamicsModel_tp::MODEL_CR3BP) {
    // Allow a few more constraints than the default
    allowedCons.push_back(Constraint_tp::JC);
    allowedCons.push_back(Constraint_tp::PSEUDOARC);
    allowedEvents.push_back(Event_tp::JC);
}//==============================================

/**
 *  \brief Copy Constructor
 *  \param m a model reference
 */
DynamicsModel_cr3bp::DynamicsModel_cr3bp(const DynamicsModel_cr3bp &m) : DynamicsModel(m) {}

/**
 *  \brief Assignment operator
 *  \param m a model reference
 */
DynamicsModel_cr3bp& DynamicsModel_cr3bp::operator =(const DynamicsModel_cr3bp &m){
	DynamicsModel::operator =(m);
	return *this;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for only the core states (i.e. simple)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp::getSimpleEOM_fcn() const{
	return &simpleEOMs;
}//==============================================

/**
 *  \brief Retrieve a pointer to the EOM function that computes derivatives
 *  for all states (i.e. full)
 */
DynamicsModel::eom_fcn DynamicsModel_cr3bp::getFullEOM_fcn() const{
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
std::vector<double> DynamicsModel_cr3bp::getPrimPos(double t, const SysData *pSysData) const{
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
void DynamicsModel_cr3bp::getPrimPos(double t, const SysData *pSysData, int pIx, double *pos) const{
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
std::vector<double> DynamicsModel_cr3bp::getPrimVel(double t, const SysData *pSysData) const{
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
void DynamicsModel_cr3bp::getPrimVel(double t, const SysData *pSysData, int pIx, double *vel) const{
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
std::vector<double> DynamicsModel_cr3bp::getStateDeriv(double t, std::vector<double> state, EOM_ParamStruct *params) const{
    if(state.size() != coreDim)
        throw Exception("DynamicsModel_cr3bp::getStateDeriv: State size does not match the core state size specified by the dynamical model");

    // Compute the acceleration
    std::vector<double> dsdt(coreDim,0);
    simpleEOMs(t, &(state[0]), &(dsdt[0]), params);
    
    return dsdt;
}//==================================================

//------------------------------------------------------------------------------------------------------
//      Simulation Engine Functions
//------------------------------------------------------------------------------------------------------

// Use the defaults

//------------------------------------------------------------------------------------------------------
//      Multiple Shooting Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Compute constraint function and partial derivative values for a constraint
 *  
 *  This function calls its relative in the DynamicsModel base class and appends additional
 *  instructions specific to the CR3BP
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param con the constraint being applied
 *  \param c the index of the constraint within the total constraint vector (which is, in
 *  turn, stored in the iteration data)
 */ 
void DynamicsModel_cr3bp::multShoot_applyConstraint(MultShootData *it, const Constraint& con, int c) const{

    // Let the base class do its thing first
    DynamicsModel::multShoot_applyConstraint(it, con, c);

    // Handle constraints specific to the CR3BP
    int row0 = it->conRows[c];

    switch(con.getType()){
        case Constraint_tp::JC:
            multShoot_targetJC(it, con, row0);
            break;
        case Constraint_tp::PSEUDOARC:
            multShoot_targetPseudoArc(it, con, row0);
        default: break;
    }
}//=========================================================

/**
 *  \brief Perform model-specific initializations on the MultShootData object
 *  \param it pointer to the object to be initialized
 */
void DynamicsModel_cr3bp::multShoot_initIterData(MultShootData *it) const{
    Arcset_cr3bp traj(static_cast<const SysData_cr3bp *>(it->nodesIn->getSysData()));
    it->propSegs.assign(it->nodesIn->getNumSegs(), traj);
}//====================================================

/**
 *  \brief Compute constraint function and partial derivative values for a Jacobi Constraint
 *
 *  \param it a pointer to the corrector's iteration data structure
 *  \param con the constraint being applied
 *  \param row0 the row this constraint begins on
 */
void DynamicsModel_cr3bp::multShoot_targetJC(MultShootData* it, const Constraint& con, int row0) const{
    std::vector<double> conData = con.getData();
    MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, con.getID());
    const SysData_cr3bp *crSys = static_cast<const SysData_cr3bp *> (it->nodesIn->getSysData());

    if(state_var.row0 == -1)
        throw Exception("DynamicsModel_cr3bp::multShoot_targetJC: Cannot constrain state that is not in the free variable vector.");

    // Compute the value of Jacobi at this node
    double mu = crSys->getMu();
    double nodeState[6];
    std::copy(&(it->X[state_var.row0]), &(it->X[state_var.row0])+6, nodeState); 

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

    it->DF_elements.push_back(Tripletd(row0, state_var.row0+0, (-2*(x + mu)*(1 - mu)/pow(d,3) - 2*(x + mu - 1)*mu/pow(r,3) + 2*x) ));
    it->DF_elements.push_back(Tripletd(row0, state_var.row0+1, (-2*y*(1 - mu)/pow(d,3) - 2*y*mu/pow(r,3) + 2*y) ));
    it->DF_elements.push_back(Tripletd(row0, state_var.row0+2, (-2*z*(1 - mu)/pow(d,3) - 2*z*mu/pow(r,3)) ));
    it->DF_elements.push_back(Tripletd(row0, state_var.row0+3, -2*vx));
    it->DF_elements.push_back(Tripletd(row0, state_var.row0+4, -2*vy));
    it->DF_elements.push_back(Tripletd(row0, state_var.row0+5, -2*vz));
}//=============================================

/**
 *  \brief Compute constraint function and partial derivative values for Pseudo Arc-Length
 *  
 *  \param it a pointer to the corrector's iteration data structure
 *  \param con the constraint being applied
 *  \param row0 the row this constraint begins on
 *  @throw Exception if the pseudo arclength constraint is not listed as the final constraint
 *  @throw Exception if the Jacobian matrix (w/o the PAL constraint) is nonsquare.
 */
void DynamicsModel_cr3bp::multShoot_targetPseudoArc(MultShootData *it, const Constraint& con, int row0) const{
    std::vector<double> conData = con.getData();

    if(row0 != it->totalCons-1)
        throw Exception("DynamicsModel_cr3bp::multShoot_targetPseudoArc: Pseudo Arc-Length constraint must be the final constraint; please re-create the nodeset accordingly");

    if(it->totalCons != it->totalFree)
        throw Exception("DynamicsModel_cr3bp::multShoot_targetPseudoArc: Jacobian matrix is not square; cannot apply pseudo arc-length");

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
        it->DF_elements.push_back(Tripletd(row0, i, N(i)));
    }
}//=============================================

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
 *  
 *  \return a pointer to a nodeset containing the corrected nodes
 */
void DynamicsModel_cr3bp::multShoot_createOutput(const MultShootData *it) const{

    // Create a nodeset with the same system data as the input
    std::vector<int> newNodeIDs;
    newNodeIDs.reserve(it->numNodes);
    
    for(int n = 0; n < it->numNodes; n++){
        MSVarMap_Obj state_var = it->getVarMap_obj(MSVar_tp::STATE, it->nodesIn->getNodeRefByIx_const(n).getID());
        std::vector<double> state;

        if(state_var.row0 == -1){
            state = it->nodesIn->getState(state_var.key.id);
        }else{
            state = std::vector<double>(it->X.begin()+state_var.row0, it->X.begin()+state_var.row0 + coreDim);
        }

        Node node(state, 0);
        node.setConstraints(it->nodesIn->getNodeRefByIx_const(n).getConstraints());

        // Add the node to the output nodeset and save the new ID
        newNodeIDs.push_back(it->nodesOut->addNode(node));
    }

    double tof;
    int newOrigID, newTermID;
    for(unsigned int s = 0; s < it->nodesIn->getNumSegs(); s++){
        const Segment &seg = it->nodesIn->getSegRefByIx_const(s);

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
        newSeg.setStateWidth(it->propSegs[s].getSegRef_const(0).getStateWidth());
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
                epoch = it->nodesOut->getEpoch(order[i].id);
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
                epoch += std::abs(it->nodesOut->getTOF(order[i].id));
            }
        }
    }

    // it->nodesOut->print();
    std::vector<Constraint> arcCons = it->nodesIn->getArcConstraints();
    for(unsigned int i = 0; i < arcCons.size(); i++){
        it->nodesOut->addConstraint(arcCons[i]);
    }
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Integrate the equations of motion for the CR3BP
 *  \param t the current time of the integration; not used for this system
 *  \param s the 42-d state vector
 *  \param sdot the 42-d state derivative vector
 *  \param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp::fullEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;

    // Extract mu from params
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp *sysData = static_cast<const SysData_cr3bp *>(paramStruct->pSysData);
    
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
 *  \brief Integrate the equations of motion for the CR3BP without the STM
 *  \param t time at integration step (unused)
 *  \param s the 6-d state vector
 *  \param sdot the 6-d state derivative vector
 *  \param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp::simpleEOMs(double t, const double s[], double sdot[], void *params){
    (void)t;
    
    // Extract mu from params
    // SysData_cr3bp *sysData = static_cast<SysData_cr3bp *>(params);
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp *sysData = static_cast<const SysData_cr3bp *>(paramStruct->pSysData);
    double mu = sysData->getMu();

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
 *  \brief Compute the location of a Lagrange point in the CR3BP
 *
 *  \param sysData pointer to an object describing the particular CR3BP
 *  \param L the Lagrange point number, 1 to 5
 *  \param tol the tolerance to use; if NAN is input, then a default value of 1e-14 will
 *  be used.
 *  \param pos a 3-element array to store the position of the Lagrange point
 *  \throws DivergeException if the Newton-Raphson process fails to converge on the 
 *  Lagrange point location
 */
void DynamicsModel_cr3bp::getEquilibPt(const SysData_cr3bp *sysData, int L, double tol, double pos[3]){
    if(L < 1 || L > 5){
        throw Exception("Invalid Lagrange Point");
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
        throw DivergeException("DynamicsModel_cr3bp::getEquilibPt: Could not converge on Lagrange point.");
    }
}//========================================

/**
 *  \brief Compute the Jacobi Constant for the CR3BP
 *
 *  \param s the state vector; only the position and velocity states are required
 *  \param mu the non-dimensional system mass ratio
 *
 *  \return the Jacobi Constant at this specific state and system
 */
double DynamicsModel_cr3bp::getJacobi(const double s[], double mu){
    double v_squared = s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);
    return 2*U - v_squared;
}//================================================

/**
 *  \brief Compute the second derivatives of the pseudo-potential function
 *
 *  \param mu the mass ratio of the system, non-dimensional
 *  \param x coordinate, non-dimensional units 
 *  \param y coordinate, non-dimensional units 
 *  \param z coordinate, non-dimensional units 
 *  \param ddots a pointer to a 6-element double array where the function will store 
 *  values for {Uxx, Uyy, Uzz, Uxy, Uxz, Uyz}. Note that Uyx = Uxy, etc.
 */
void DynamicsModel_cr3bp::getUDDots(double mu, double x, double y, double z, double* ddots){
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

}// END of Astrohelion namespace