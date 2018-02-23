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

#include "Arcset_cr3bp_lt.hpp"
#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
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
 *  @brief Construct a CR3BP Low-Thrust, Velocity Pointing Dynamic DynamicsModel
 */
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt() : DynamicsModel_cr3bp() {
    coreDim = 7;
    extraDim = 0;
    allowedEvents.push_back(Event_tp::MASS);
}//==============================================

/**
 *  @brief Copy Constructor
 *  @param m a model reference
 */ 
DynamicsModel_cr3bp_lt::DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt &m) : DynamicsModel_cr3bp(m) {}

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
 *  @brief Retrieve the state derivative
 *  @details Evaluate the equations of motion to compute the state time-derivative at 
 *  the specified time and state
 * 
 *  @param t time parameter
 *  @param state state vector
 *  @param params structure containing parameters relevant to the integration
 *  @return the time-derivative of the state vector
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
 *  @brief Create default events for a simulation run
 *  @details These events are intended to prevent numerical issues, e.g., to avoid singularities.
 * 
 *  @param pSys pointer to system data object
 *  @return A vector of events to use in the simulation
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
 *  @brief Perform model-specific initializations on the MultShootData object
 *  @param it a reference to the object to be initialized
 */
void DynamicsModel_cr3bp_lt::multShoot_initIterData(MultShootData& it) const{
    Arcset_cr3bp_lt traj(static_cast<const SysData_cr3bp_lt *>(it.pArcIn->getSysData()));
    it.propSegs.assign(it.pArcIn->getNumSegs(), traj);
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Static Calculation Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP
 *  @param t the current time of the integration
 *  @param s the state vector passed in from the SimEngine. This vector includes
 *  the core states, STM states, extra states, and control states, in that order.
 *  @param sdot the 43-d state derivative vector
 *  @param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::fullEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt *>(paramStruct->pSysData);
    const ControlLaw_cr3bp_lt *law = static_cast<const ControlLaw_cr3bp_lt *>(paramStruct->pCtrlLaw);
    const unsigned int coreDim = pSys->getDynamicsModel()->getCoreStateSize();

    double mu = pSys->getMu();           // nondimensional mass ratio

    // compute distance to primaries and velocity magnitude
    double r13 = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r23 = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Retrieve the control law acceleration values
    double control_accel[3] = {0};
    if(law){
        law->getLaw_Output(t, s, pSys, control_accel, 3);
    }

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(r13,3) - mu*(s[0]-1+mu)/pow(r23,3) + control_accel[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(r13,3) - mu*s[1]/pow(r23,3) + control_accel[1];
    sdot[5] = -(1-mu)*s[2]/pow(r13,3) - mu*s[2]/pow(r23,3) + control_accel[2];

    sdot[6] = law ? law->get_dmdt(t, s, pSys) : 0;

    // Save any time-derivatives of the control states
    const unsigned int ctrlDim = law ? law->getNumStates() : 0;
    const unsigned int ctrlOutDim = law ? law->getNumOutputs() : 0;
    if(law && ctrlDim > 0){
        std::vector<double> control_stateDeriv(ctrlDim, 0);
        law->getLaw_StateDeriv(t, s, pSys, &(control_stateDeriv.front()), ctrlDim);

        std::copy(control_stateDeriv.begin(), control_stateDeriv.begin() + ctrlDim, sdot + coreDim);
    }

    
    const unsigned int stmSide = (coreDim + ctrlDim);
    std::vector<double> A(stmSide*stmSide, 0);

    // Velocity relationships
    A[0*stmSide + 3] = 1;   // d/dvx (dx/dt)
    A[1*stmSide + 4] = 1;   // d/dvy (dy/dt)
    A[2*stmSide + 5] = 1;   // d/dvz (dz/dt)

    // Uxx = d/dx (dvx/dt)
    A[3*stmSide + 0] = 1 - (1-mu)/pow(r13,3) - mu/pow(r23,3) + 3*(1-mu)*pow((s[0] + mu),2)/pow(r13,5) + 
        3*mu*pow((s[0] + mu - 1), 2)/pow(r23,5);
    // Uxy = d/dy (dvx/dt)
    A[3*stmSide + 1] = 3*(1-mu)*(s[0] + mu)*s[1]/pow(r13,5) + 3*mu*(s[0] + mu - 1)*s[1]/pow(r23,5);
    // Uxz = d/dz (dvx/dt)
    A[3*stmSide + 2] = 3*(1-mu)*(s[0] + mu)*s[2]/pow(r13,5) + 3*mu*(s[0] + mu - 1)*s[2]/pow(r23,5);

    // Uyy = d/dy (dvy/dt)
    A[4*stmSide + 1] = 1 - (1-mu)/pow(r13,3) - mu/pow(r23,3) + 3*(1-mu)*s[1]*s[1]/pow(r13,5) + 3*mu*s[1]*s[1]/pow(r23,5);

    // Uyz = d/dz (dvy/dt)
    A[4*stmSide + 2] = 3*(1-mu)*s[1]*s[2]/pow(r13,5) + 3*mu*s[1]*s[2]/pow(r23,5);

    // Uzz = d/dz (dvz/dt)
    A[5*stmSide + 2] = -(1-mu)/pow(r13,3) - mu/pow(r23,3) + 3*(1-mu)*s[2]*s[2]/pow(r13,5) + 3*mu*s[2]*s[2]/pow(r23,5);

    // Symmetry
    A[4*stmSide + 0] = A[3*stmSide + 1];
    A[5*stmSide + 0] = A[3*stmSide + 2];
    A[5*stmSide + 1] = A[4*stmSide + 2];

    A[3*stmSide + 4] = 2;   // d/dvy (dvx/dt)
    A[4*stmSide + 3] = -2;   // d/dvx (dvy/dt)

    unsigned int r = 0, c = 0, k = 0;
    if(law){
        // Get partial derivatives of acceleration terms (which are part of EOMS 3, 4, 5) w.r.t. all state variables
        std::vector<double> law_accelPartials(ctrlOutDim*coreDim, 0);
        
        law->getLaw_OutputPartials(t, s, pSys, &(law_accelPartials.front()), law_accelPartials.size());

        // Add the control output partials to the existing partials (control outputs are part of core state EOMs)
        for(r = 3; r < 3+ctrlOutDim; r++){
            for(c = 0; c < coreDim; c++){
                A[r*stmSide + c] += law_accelPartials.at((r - 3)*coreDim + c);
            }
        }

        if(ctrlDim > 0){
            std::vector<double> law_eomPartials(coreDim*ctrlDim, 0), law_stateDerivPartials(ctrlDim*(coreDim+ctrlDim), 0);

            law->getLaw_StateDerivPartials(t, s, pSys, &(law_stateDerivPartials.front()), law_stateDerivPartials.size());
            law->getLaw_EOMPartials(t, s, pSys, &(law_eomPartials.front()), law_eomPartials.size());

            // Assign the partial derivatives of the time derivatives of the control states w.r.t. all core and control states
            for(r = coreDim; r < coreDim + ctrlDim; r++){
                for(c = 0; c < coreDim + ctrlDim; c++){
                    A[r*stmSide + c] = law_stateDerivPartials.at((r - coreDim)*(coreDim + ctrlDim) + c);
                }
            }

            // Assign the partial derivatives of the EOMs w.r.t. control states
            for(r = 0; r < coreDim; r++){
                for(c = coreDim; c < coreDim + ctrlDim; c++){
                    A[r*stmSide + c] = law_eomPartials.at(r*ctrlDim + c - coreDim);
                }
            }
        }
    }

    // Do the matrix multiplication Phi_dot = A*Phi to compute
    // the STM time-derivative.
    // * Phi element 0 is stored in s[coreDim + ctrlDim], size stmSide*stmSide
    // * Phi_dot element 0 is stored in the same location, but in the sdot array
    for(r = 0; r < stmSide; r++){
        for(c = 0; c < stmSide; c++){
            sdot[coreDim + ctrlDim + r*stmSide + c] = 0;
            for(k = 0; k < stmSide; k++){
                sdot[coreDim + ctrlDim + r*stmSide + c] += A[r*stmSide + k]*s[coreDim + ctrlDim + k*stmSide + c];
            }
        }
    }

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  @param t time at integration step (unused)
 *  @param s the state vector passed in from the SimEngine. This vector includes
 *  the core states and control states, in that order.
 *  @param sdot the state derivative vector
 *  @param params points to an EOM_ParamStruct object
 */
int DynamicsModel_cr3bp_lt::simpleEOMs(double t, const double s[], double sdot[], void *params){
    EOM_ParamStruct *paramStruct = static_cast<EOM_ParamStruct *>(params);
    const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt *>(paramStruct->pSysData);
    const ControlLaw_cr3bp_lt *law = static_cast<const ControlLaw_cr3bp_lt *>(paramStruct->pCtrlLaw);

    double mu = pSys->getMu();           // nondimensional mass ratio

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (s[0]+mu)*(s[0]+mu) + s[1]*s[1] + s[2]*s[2] );
    double r = sqrt( (s[0]-1+mu)*(s[0]-1+mu) + s[1]*s[1] + s[2]*s[2] );

    // Retrieve the control law acceleration values
    double control_accel[3] = {0};
    if(law)
        law->getLaw_Output(t, s, pSys, control_accel, 3);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*s[4] + s[0] - (1-mu)*(s[0]+mu)/pow(d,3) - mu*(s[0]-1+mu)/pow(r,3) + control_accel[0];
    sdot[4] = -2*s[3] + s[1] - (1-mu) * s[1]/pow(d,3) - mu*s[1]/pow(r,3) + control_accel[1];
    sdot[5] = -(1-mu)*s[2]/pow(d,3) - mu*s[2]/pow(r,3) + control_accel[2];
    
    sdot[6] = law ? law->get_dmdt(t, s, pSys) : 0;

    // Save any derivatives of the control states
    if(law){
        if(unsigned int ctrl_dim = law->getNumStates() > 0){
            std::vector<double> control_stateDeriv(ctrl_dim, 0);
            law->getLaw_StateDeriv(t, s, pSys, &(control_stateDeriv.front()), ctrl_dim);

            std::copy(control_stateDeriv.begin(), control_stateDeriv.begin() + ctrl_dim, sdot+7);
        }
    }

    return GSL_SUCCESS;
}//====================================================

/**
 * @brief      Compute the location of the equilibria given a planar low-thrust vector
 *
 * @param  pSys  Pointer to the system data object
 * @param  L     Lagrange point number [1 - 5]
 * @param  f     Nondimensional thrust magnitude
 * @param  tol   The tolerance
 * @param  zac   Pointer to a vector that stores the equilibria as rows of
 *               values: [alpha, x, y]
 * @param  verb  The verbosity with which to carry out the numerical processes
 */
void DynamicsModel_cr3bp_lt::getEquilibPt(const SysData_cr3bp_lt *pSys, int L, double f,
    double tol, std::vector<double> *zac, Verbosity_tp verb){

    unsigned int maxCount = 20;
    unsigned int maxOuterCount = 10000;
    double maxDiffAlpha = PI/360.0;
    double minStep_x = f*1e-6;
    double minStep_y = f*1e-6;
    double minStep_alpha = PI/500.0;
    double stepChangeFactor = 2;

    if(L < 1 || L > 5)
        throw Exception("DynamicsModel_cr3bp_lt::getEquilibPt: input L value is invalid");

    std::vector<double>& zacRef = *zac;
    if(f == 0){
        zac->assign(3, 0);
        double LPt[3];
        DynamicsModel_cr3bp::getEquilibPt(pSys, L, tol, LPt);
        zacRef[1] = LPt[1];
        zacRef[2] = LPt[2];
    }else{
        double mu = pSys->getMu();
        double initSol[3] = {0};    //[alpha, x, y]

        // First, compute location at alpha = 0 (or pi)
        if(L > 3){
            if(f < mu){
                double r13 = pow( (1-mu)/(1-mu+f), 1.0/3.0 );
                double r23 = pow( 1 - f/mu, -1.0/3.0);
                double theta = acos((1 + r13*r13 - r23*r23)/(2*r13));

                if(r13 - 1 < r23 && r13 + 1 > r23 && r13 + r23 > 1){
                    if(!std::isnan(theta)){
                        int sgn = L == 4 ? 1 : -1;
                        initSol[1] = r13*cos(sgn*theta) - mu;
                        initSol[2] = r13*sin(sgn*theta);
                    }else{
                        char msg[128];
                        sprintf(msg, "DynamicsModel_cr3bp_lt::getEquilibPt: Did not catch imaginary number in L%d at alpha = 0", L);
                        throw Exception(msg);
                    }
                }
            }else{
                // Start at alpha = pi instead of alpha = 0
                if(f < 1-mu){
                    double r13 = pow( (1.0 - mu)/(1.0 - mu - f), 1.0/3.0);
                    double r23 = pow( 1.0 + f/mu, -1.0/3.0);

                    if(r13 - 1 < r23 && r13 + 1 > r23 && r13 + r23 > 1){
                        double theta = acos( (1 + r13*r13 - r23*r23)/(2*r13) );
                        int sgn = L == 4 ? 1 : -1;

                        initSol[0] = PI;
                        initSol[1] = r13*cos(sgn*theta) - mu;
                        initSol[2] = r13*sin(sgn*theta);
                    }else{
                        char msg[128];
                        sprintf(msg, "DynamicsModel_cr3bp_lt::getEquilibPt: Could not compute initial value for L%d", L);
                        throw Exception(msg);
                    }
                }else{
                    // This shouldn't happen for "reasonable" values of f
                    char msg[128];
                    sprintf(msg, "DynamicsModel_cr3bp_lt::getEquilibPt: f = %e is too large for analytical solution for L4/5", f);
                    throw Exception(msg);
                }
            }

        }else{  // end of L > 3; thus L = 1, 2, or 3
            // No analytical solution is available for the locations of L1, L2, L3, so use
            // a Newton-Raphson algorithm to solve.
            //
            // Assume alpha0 = 0 for these

            double gamma = L > 2 ? 1 : 0.1;
            double gamma_prev = -999;
            unsigned int count = 0;

            while(std::abs(gamma - gamma_prev) > tol && count < maxCount){
                count++;
                gamma_prev = gamma;
                switch(L){
                    case 1:
                        gamma -= (mu/(gamma*gamma) - (1-mu)/((1-gamma)*(1-gamma)) - gamma - mu + 1 + f) / (-2*mu/pow(gamma, 3) - 2*(1-mu)/pow(1-gamma, 3) - 1);
                        break;
                    case 2:
                        gamma -= (-mu/(gamma*gamma) - (1-mu)/((1+gamma)*(1+gamma)) + gamma - mu + 1 + f) / (2*mu/pow(gamma, 3) + 2*(1-mu)/pow(1+gamma, 3) + 1);
                        break;
                    case 3:
                        gamma -= (mu/((1+gamma)*(1+gamma)) + (1-mu)/(gamma*gamma) - gamma - mu + f) / (-2*mu/pow(1+gamma, 3) - 2*(1-mu)/pow(gamma, 3) - 1);
                        break;
                    default:
                    {
                        char msg[128];
                        sprintf(msg, "DynamicsModel_cr3bp_lt::getEquilibPt: Invalid L = %d in Newton-Raphson process", L);
                        throw Exception(msg);
                    }
                }
            }

            if(std::abs(gamma - gamma_prev) > tol){
                char msg[128];
                sprintf(msg, "DynamicsModel_cr3bp_lt::getEquilibPt: L%d Newton process did not converge for alpha = 0", L);
                throw Exception(msg);
            }

            if(L < 3){
                int sgn = L == 1 ? -1 : 1;
                initSol[1] = 1-mu + sgn*gamma;
            }else{
                initSol[1] = -mu - gamma;
            }
        }

        // Ok - we have an initial solution from analytical or semi-analytical methods
        // Now, set up continuation process
        double slope = L < 4 ? 999 : 0;             // Force vertical step for collinear points, hoirzontal step for triangular points
        // double alphaSpan = L < 4 ? PI : 2*PI;       // Required angle range that must be computed (can use symmetry on collinear points)
        double step_x = f*1e-3, step_y = f*1e-3;    // Initial step sizes for continuation
        double step_alpha = (PI/100.0)*cos(initSol[0]); // Initial step size for continuation
        double diffAlpha = 0;                       // Change in alpha between solution iterations
        bool doAlphaStep = false;

        double alpha = initSol[0];
        double x = initSol[1];
        double y = initSol[2];

        double err = 999;
        unsigned int innerCount = 0;

        // Store initial solution
        zac->push_back(alpha);
        zac->push_back(x);
        zac->push_back(y);

        // Begin loop
        unsigned int outerCount = 0;
        bool reachedEnd = false;
        while( !reachedEnd && outerCount < maxOuterCount){

            if(L >= 4){
                // Old check: std::abs(alpha - initSol[0]) < alphaSpan - PI/200.0
                // Check for looped all the way around
                reachedEnd = (outerCount > 15) && std::abs(x - initSol[1]) < stepChangeFactor*std::abs(step_x) && 
                    std::abs(y - initSol[2]) < stepChangeFactor*std::abs(step_y) &&
                    std::abs(sin(alpha) - sin(initSol[0])) < stepChangeFactor*std::abs(step_alpha);
            }else{
                // Check for x-axis crossing
                reachedEnd = (outerCount > 15) && zacRef[zac->size() - 1]*zacRef[zac->size() - 4] < 0;
            }

            if(std::abs(slope) > 1.0){
                y += step_y;        // Take a step in y (this value is fixed)
                x += step_y/slope;  // Linear approximation

                // Newton-Raphson to solve for x and alpha given y
                err = 999;
                innerCount = 0;
                while(err >= tol && innerCount < maxCount){
                    double Uddots[6] = {0};
                    DynamicsModel_cr3bp::getUDDots(mu, x, y, 0, Uddots);
                    double r13 = sqrt((x+mu)*(x+mu) + y*y);
                    double r23 = sqrt((x-1+mu)*(x-1+mu) + y*y);

                    // 2D acceleration vector - we want to find the zeros of this guy
                    double F[] = {x - (1-mu)*(x+mu)/pow(r13, 3) - mu*(x-1+mu)/pow(r23,3) + f*cos(alpha),
                                    y - (1-mu)*y/pow(r13,3) - mu*y/pow(r23,3) + f*sin(alpha)};
                    // A is the Jacobian of F w.r.t. x and alpha; use analytical inverse; this is the determinant
                    double detA = Uddots[0]*f*cos(alpha) + f*sin(alpha)*Uddots[3];

                    // Update p -= inv(A)*F where p = [x; alpha]
                    x -= (f*cos(alpha)*F[0] + f*sin(alpha)*F[1])/detA;
                    alpha -= (-Uddots[3]*F[0] + Uddots[0]*F[1])/detA;

                    err = sqrt(F[0]*F[0] + F[1]*F[1]);
                    innerCount++;
                }

                // Check the change in angle
                diffAlpha = alpha - zacRef[zac->size()-3];

                // Decrease step size if Newton-Raphson does not converge OR if the change in angle is too large (makes pictures look better)
                if(err >= tol || std::abs(diffAlpha) > maxDiffAlpha){
                    if(std::abs(diffAlpha) > maxDiffAlpha){
                        alpha = zacRef[zac->size()-3]; // Reset to previous value
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  |diffAlpha = %f| > %f\n", diffAlpha, maxDiffAlpha);
                    }

                    if(std::abs(step_y/stepChangeFactor) > minStep_y){
                        y = zacRef[zac->size() - 1];    // Reset to previous value
                        x = zacRef[zac->size() - 2];
                        step_y /= stepChangeFactor;                  // Decrease step size
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  decreasing step_y to %e\n", step_y);
                    }else{
                        // Can't step any smaller in y, so try stepping in alpha
                        doAlphaStep = true;
                    }
                }else{
                    // No errors and angular separation is good!
                    // Save data, update slope
                    zac->push_back(alpha);
                    zac->push_back(x);
                    zac->push_back(y);

                    if(zac->size() >= 6){
                        // Adjust step_x as we move along y for the future when stepping changes to the horizontal direction
                        step_x = zacRef[zac->size() - 2] - zacRef[zac->size() - 5];
                        step_alpha = diffAlpha;
                        slope = (zacRef[zac->size() - 1] - zacRef[zac->size() - 4])/step_x;
                        if(std::isnan(slope)){ slope = 999; }

                        printVerb(verb >= Verbosity_tp::SOME_MSG, "STEP_Y: alpha = %f, [%f, %f]; step_x = %f, step_y = %f, step_alpha = %f, slope = %f\n",
                            alpha, x, y, step_x, step_y, step_alpha, slope);
                    }

                    if(innerCount < 5 && maxDiffAlpha/std::abs(diffAlpha) > stepChangeFactor && std::abs(slope) > 1.0){
                        step_y *= stepChangeFactor;    // Increase step size to expedite process
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  increasing step_y = %e\n", step_y);
                    }
                }
            }else{
                x += step_x;            // slope <= 1, so step in x
                y += step_x*slope;      // Linear approximation

                // Newton-Raphson to solve for y and alpha given x
                err = 999;
                innerCount = 0;
                while(err >= tol && innerCount < maxCount){
                    double Uddots[6] = {0};
                    DynamicsModel_cr3bp::getUDDots(mu, x, y, 0, Uddots);
                    double r13 = sqrt((x+mu)*(x+mu) + y*y);
                    double r23 = sqrt((x-1+mu)*(x-1+mu) + y*y);

                    // 2D acceleration vector - we want to find the zeros of this guy
                    double F[] = {x - (1-mu)*(x+mu)/pow(r13, 3) - mu*(x-1+mu)/pow(r23,3) + f*cos(alpha),
                                    y - (1-mu)*y/pow(r13,3) - mu*y/pow(r23,3) + f*sin(alpha)};

                    // A is the Jacobian of F w.r.t. y and alpha; use analytical inverse; this is the determinant
                    double detA = Uddots[3]*f*cos(alpha) + f*sin(alpha)*Uddots[1];

                    // Update p -= inv(A)*F where p = [y; alpha]
                    y -= (f*cos(alpha)*F[0] + f*sin(alpha)*F[1])/detA;
                    alpha -= (-Uddots[1]*F[0] + Uddots[3]*F[1])/detA;

                    err = sqrt(F[0]*F[0] + F[1]*F[1]);
                    innerCount++;
                }

                // Check the change in angle
                diffAlpha = alpha - zacRef[zac->size()-3];

                // Decrease step size if Newton-Raphson does not converge OR if the change in angle is too large (makes pictures look better)
                if(err >= tol || std::abs(diffAlpha) > maxDiffAlpha){
                    if(std::abs(diffAlpha) > maxDiffAlpha){
                        alpha = zacRef[zac->size()-3]; // Reset to previous value
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  |diffAlpha = %f| > %f\n", diffAlpha, maxDiffAlpha);
                    }

                    if(std::abs(step_x/stepChangeFactor) > minStep_x){
                        y = zacRef[zac->size() - 1];   // Reset to previous values
                        x = zacRef[zac->size() - 2];
                        step_x /= stepChangeFactor;              // Decrease step size
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  decreasing step_x to %e\n", step_x);
                    }else{
                        // Can't step any smaller in x, so try stepping in alpha
                        doAlphaStep = true;
                    }
                }else{
                    // No errors and angular separation is good!
                    // Save data, update slope
                    zac->push_back(alpha);
                    zac->push_back(x);
                    zac->push_back(y);

                    if(zac->size() >= 6){
                        // Adjust step_y as we move along x for the future when stepping changes in the vertical direction
                        step_y = zacRef[zac->size() - 1] - zacRef[zac->size() - 4];
                        step_alpha = diffAlpha;
                        slope = step_y/(zacRef[zac->size() - 2] - zacRef[zac->size() - 5]);
                        if(std::isnan(slope)){ slope = 999; }

                        printVerb(verb >= Verbosity_tp::SOME_MSG, "STEP_X: alpha = %f, [%f, %f]; step_x = %f, step_y = %f, step_alpha = %f, slope = %f\n",
                            alpha, x, y, step_x, step_y, step_alpha, slope);
                    }

                    if(innerCount < 5 && maxDiffAlpha/std::abs(diffAlpha) > stepChangeFactor && std::abs(slope) < 1.0){
                        step_x *= stepChangeFactor;    // Increase step size to expedite process
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  increasing step_x = %e\n", step_x);
                    }
                }
            }// End of if/else for abs(slope) > 1

            if(doAlphaStep){
                alpha = zacRef[zac->size() - 3] + step_alpha;   // Reset and take a step

                err = 999;
                innerCount = 0;
                while(err >= tol && innerCount < maxCount){
                    double Uddots[6] = {0};
                    DynamicsModel_cr3bp::getUDDots(mu, x, y, 0, Uddots);
                    double r13 = sqrt((x+mu)*(x+mu) + y*y);
                    double r23 = sqrt((x-1+mu)*(x-1+mu) + y*y);

                    // 2D acceleration vector - we want to find the zeros of this guy
                    double F[] = {x - (1-mu)*(x+mu)/pow(r13, 3) - mu*(x-1+mu)/pow(r23,3) + f*cos(alpha),
                                    y - (1-mu)*y/pow(r13,3) - mu*y/pow(r23,3) + f*sin(alpha)};

                    // A is the Jacobian of F w.r.t. x and y; use analytical inverse; this is the determinant
                    double detA = Uddots[0]*Uddots[1] - Uddots[3]*Uddots[3];

                    // Update p -= inv(A)*F where p = [x; y]
                    x -= (Uddots[1]*F[0] - Uddots[3]*F[1])/detA;
                    y -= (-Uddots[3]*F[0] + Uddots[0]*F[1])/detA;

                    err = sqrt(F[0]*F[0] + F[1]*F[1]);
                    innerCount++;
                }

                if(err >= tol){
                    if(std::abs(step_alpha/stepChangeFactor) > minStep_alpha){
                        step_alpha /= stepChangeFactor;
                        printVerb(verb >= Verbosity_tp::SOME_MSG, "  decreasing step_alpha to %e\n", step_alpha);
                    }else{
                        // Could not converge with alpha step... give up on this one
                        printWarn("Could not converge L%d with alpha step... ending\n", L);
                        break;
                    }
                }else{
                    // Save data, update slope
                    zac->push_back(alpha);
                    zac->push_back(x);
                    zac->push_back(y);

                    if(zac->size() >= 6){
                        step_x = zacRef[zac->size() - 2] - zacRef[zac->size() - 5];
                        step_y = zacRef[zac->size() - 1] - zacRef[zac->size() - 4];
                        slope = step_y/step_x;

                        if(std::isnan(slope)){ slope = 999; }

                        diffAlpha = zacRef[zac->size() - 6] - zacRef[zac->size() - 3];
                        doAlphaStep = false;

                        printVerb(verb >= Verbosity_tp::SOME_MSG, "STEP_A: alpha = %f, [%f, %f]; step_x = %f, step_y = %f, step_alpha = %f, slope = %f\n",
                            alpha, x, y, step_x, step_y, step_alpha, slope);
                    }
                }
            }

            outerCount++;
        }// End of loop

        if(reachedEnd)
            printVerb(verb >= Verbosity_tp::SOME_MSG, "L%d loop ended because reachedEnd = true\n", L);
        else if(outerCount < maxOuterCount)
            printVerb(verb >= Verbosity_tp::SOME_MSG, "L%d loop ended due to iteration max-out\n", L);

        
        if(L <= 3){
            std::vector<double> sym;
            for(unsigned int i = 0; i < zac->size()/3.0; i++){
                zacRef[3*i] = wrapToPi(zacRef[3*i]);    // Wrap angles to [-pi, pi]

                // Leverage Symmetry        
                sym.insert(sym.begin(), -1.0*zacRef[3*i+2]);
                sym.insert(sym.begin(),  zacRef[3*i+1]);
                sym.insert(sym.begin(), -1.0*zacRef[3*i+0]);
            }
            zac->insert(zac->begin(), sym.begin(), sym.end());
        }else{
            for(unsigned int i = 0; i < zac->size()/3.0; i++){
                zacRef[3*i] = wrapToPi(zacRef[3*i]);    // Wrap angles to [-pi, pi]
            }
        }
    }// End of if/else f == 0
}//====================================================

/**
 * @brief Compute the CR3BP-LT Hamiltonian at the specified state and time
 * 
 * @param t nondimensional time value associated with the state vector
 * @param q pointer to the full state (core + control)
 * @param pSys pointer to the system data object
 * @param pLaw pointer to the control law object
 * @return the CR3BP-LT Hamiltonian value at the specified state and time
 */
double DynamicsModel_cr3bp_lt::getHamiltonian(double t, const double *q, const SysData_cr3bp_lt *pSys, const ControlLaw_cr3bp_lt *pLaw){
    if(pSys ==  nullptr)
        throw Exception("DynamicsModel_cr3bp_lt::getHamiltonian: Cannot proceed with SysData = nullptr");

    if(pLaw == nullptr)
        throw Exception("DynamicsModel_cr3bp_lt::getHamiltonian: Cannot proceed with ControlLaw = nullptr");

    // Get the control law output (acceleration vector)
    Eigen::Vector3d a;
    pLaw->getLaw_Output(t, q, pSys, a.data(), 3);

    // Compute other useful parameters
    double mu = pSys->getMu();
    double r13 = sqrt((q[0] + mu)*(q[0] + mu) + q[1]*q[1] + q[2]*q[2]);
    double r23 = sqrt((q[0] - 1 + mu)*(q[0] -1 + mu) + q[1]*q[1] + q[2]*q[2]);

    // Compute the CR3BP-LT Hamiltonian
    return 0.5*(q[3]*q[3] + q[4]*q[4] + q[5]*q[5]) - 0.5*(q[0]*q[0] + q[1]*q[1]) -
        (1-mu)/r13 - mu/r23 - a(0)*q[0] - a(1)*q[1] - a(2)*q[2];
}//====================================================

/**
 *  @brief Construct a new control law and allocated it on the stack.
 *  @details Each dynamic model will return a pointer to the specific control
 *  law applicable to the system / model
 *  
 *  @param id the control law ID
 *  @param params the parameters that create the control law
 *  
 *  @return A pointer to a control law object. The object has been allocated
 *  on the stack so the delete() function must be employed to free the memory
 */
ControlLaw* DynamicsModel_cr3bp_lt::createControlLaw(unsigned int id,
    const std::vector<double> &params) const{

    return new ControlLaw_cr3bp_lt(id, params);
}//=====================================================

/**
 *  @brief Determine whether or not the model supports a specific control law
 *  @details By default, a DynamicsModel does not support any control laws.
 * 
 *  @param pLaw Pointer to a control law
 *  @return Whether or not the Dynamics model supports the specified control law
 */
bool DynamicsModel_cr3bp_lt::supportsControl(const ControlLaw *pLaw) const{
    // Does support control laws
    (void) pLaw;
    return true;
}//====================================================





}// END of Astrohelion namespace


