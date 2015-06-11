/**
 *  @file tpat_calculations.cpp
 *
 *   This library contains functions to numerically integrate the equations of 
 *   motion associated with various mathematical models. It also contains 
 *   miscellaneous functions to compute other quantities, like Jacobi Constant
 */

#include "tpat_calculations.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_bcr4bpr_traj.hpp"
#include "tpat_constants.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_cr3bp_traj.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_trajectory.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <iostream>


using namespace std;

//-----------------------------------------------------
//      Equations of Motion
//-----------------------------------------------------

/**
*   @brief Integrate the equations of motion for the CR3BP
*   @param t the current time of the integration
*   @param s the 42-d state vector
*   @param sdot the 42-d state derivative vector
*   @param *params pointer to extra parameters required for integration. For this
*   function, the pointer points to a cr3bp system data object
*/
int cr3bp_EOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params
    tpat_cr3bp_sys_data *sysData = static_cast<tpat_cr3bp_sys_data *>(params);
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
    cr3bp_getUDDots(mu, x, y, z, ddots);

    double a_data[] = { 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1,
                        ddots[0], ddots[3], ddots[4], 0, 2, 0,
                        ddots[3], ddots[1], ddots[5], -2, 0, 0,
                        ddots[4], ddots[5], ddots[2], 0, 0, 0};
    tpat_matrix A(6,6,a_data);

    // Copy the STM states into a sub-array
    double stmElements[36];
    copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    tpat_matrix phi(6,6, stmElements);
    
    // Compute derivative of STM
    tpat_matrix phiDot = A*phi;
    double *phiDotData = phiDot.getDataPtr();

    // Copy the elements of phiDot into the derivative array
    copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *   @brief Integrate the equations of motion for the CR3BP without the STM
 *   @param t time at integration step (unused)
 *   @param s the 6-d state vector
 *   @param sdot the 6-d state derivative vector
 *   @param params points to a cr3bp system data object
 */
int cr3bp_simple_EOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params
    tpat_cr3bp_sys_data *sysData = static_cast<tpat_cr3bp_sys_data *>(params);
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
 *   @brief Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   @param t time at integration step (unused)
 *   @param s the 48-d state vector
 *   @param sdot the 48-d state derivative vector
 *   @param params points to additional integration parameters wrapped in an 
 *  <tt>tpat_bcr4bpr_sys_data</tt> data object.
 */
int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    tpat_bcr4bpr_sys_data *sysData = static_cast<tpat_bcr4bpr_sys_data *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, *sysData, primPosData);
    tpat_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    copy(s, s+3, r_data);
    tpat_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    copy(s+3, s+6, v_data);
    tpat_matrix v(3,1,v_data);

    // Create relative position vectors between s/c and primaries
    tpat_matrix r_p1 = r - trans(primPos.getRow(0));
    tpat_matrix r_p2 = r - trans(primPos.getRow(1));
    tpat_matrix r_p3 = r - trans(primPos.getRow(2));
    double d1 = norm(r_p1);
    double d2 = norm(r_p2);
    double d3 = norm(r_p3);
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Create C-matrix
    double c[] = {0, 2*k, 0, -2*k, 0, 0, 0, 0, 0};
    tpat_matrix C(3,3,c);

    // truncated position vector used in EOMs
    double r_trunc_data[3] = {0};
    copy(s, s+2, r_trunc_data);
    tpat_matrix r_trunc(3,1,r_trunc_data);

    // Compute acceleration using matrix math
    tpat_matrix accel(3,1);
    accel = C*v + k*k*r_trunc - (1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2, 3) - 
            nu*r_p3/pow(d3, 3);

    // Compute psuedo-potential
    double dxdx = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(0),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) -
                3*pow(r_p3.at(0),2)/pow(d3,5));
    double dxdy = (1/k - mu)*3*r_p1.at(0)*r_p1.at(1)/pow(d1,5) +
            (mu - nu)*3*r_p2.at(0)*r_p2.at(1)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);
    double dxdz = (1/k - mu)*3*r_p1.at(0)*r_p1.at(2)/pow(d1,5) +
            (mu - nu)*3*r_p2.at(0)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);
    double dydy = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(1),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(1),2)/pow(d2,5)) - nu*(1/pow(d3,3) -
            3*pow(r_p3.at(1),2)/pow(d3,5));
    double dydz = (1/k - mu)*3*r_p1.at(1)*r_p1.at(2)/pow(d1,5) +
            (mu - nu)*3*r_p2.at(1)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);
    double dzdz = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(2),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) -
            3*pow(r_p3.at(2),2)/pow(d3,5));

    // Create A matrix for STM derivative
    double aData[] = {  0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1, 
                        dxdx, dxdy, dxdz, c[0], c[1], c[2],
                        dxdy, dydy, dydz, c[3], c[4], c[5],
                        dxdz, dydz, dzdz, c[6], c[7], c[8]};
    tpat_matrix A(6,6, aData);
    
    // Compute the STM derivative
    double phiData[36];
    copy(s+6, s+42, phiData);
    tpat_matrix Phi(6, 6, phiData);
    tpat_matrix PhiDot = A*Phi;

    // Compute partials of state w.r.t. primary positions; dont' compute partials
    // for P1 because its velocity is zero in the rotating frame
    double dfdr2[18] = {0};   double dfdr3[18] = {0};

    dfdr2[9] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);        //dxdx2
    dfdr2[10] = 3*r_p2.at(0)*r_p2.at(1)/pow(d2,5);                  //dxdy2
    dfdr2[11] = 3*r_p2.at(0)*r_p2.at(2)/pow(d2,5);                  //dxdz2
    dfdr2[13] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);       //dydy2
    dfdr2[14] = 3*r_p2.at(1)*r_p2.at(2)/pow(d2,5);                  //dydz2
    dfdr2[17] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);       //dzdz2

    dfdr2[12] = dfdr2[10];      // Fill in symmetric matrix
    dfdr2[15] = dfdr2[11];
    dfdr2[16] = dfdr2[14];

    dfdr3[9] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);        //dxdx3
    dfdr3[10] = 3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);                  //dxdy3
    dfdr3[11] = 3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);                  //dxdz3
    dfdr3[13] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);       //dydy3
    dfdr3[14] = 3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);                  //dydz3
    dfdr3[17] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);       //dzdz3

    dfdr3[12] = dfdr3[10];      // Fill in symmetric matrix
    dfdr3[15] = dfdr3[11];
    dfdr3[16] = dfdr3[14];

    tpat_matrix DfDr2(6,3, dfdr2);
    tpat_matrix DfDr3(6,3, dfdr3);

    // Scale by constants
    DfDr2 *= -1*(mu-nu);
    DfDr3 *= -1*nu;
    
    // Pull the state derivative w.r.t. Epoch time from the large state vector; create column vector
    double dqdT_data[6] = {0};
    copy(s+42,s+48, dqdT_data);
    tpat_matrix dqdT(6,1, dqdT_data);

    // Get the velocity of the primaries
    double primVelData[9] = {0};
    bcr4bpr_getPrimaryVel(t, *sysData, primVelData);
    tpat_matrix primVel(3,3, primVelData);

    // Compute derivative of dqdT
    tpat_matrix dot_dqdT = A*dqdT + DfDr2*trans(primVel.getRow(1)) + DfDr3*trans(primVel.getRow(2));

    // Save derivatives to output vector
    double *accelPtr = accel.getDataPtr();
    double *phiDotPtr = PhiDot.getDataPtr();
    double *dqdtDotPtr = dot_dqdT.getDataPtr();

    copy(s+3, s+6, sdot);
    copy(accelPtr, accelPtr+3, sdot+3);
    copy(phiDotPtr, phiDotPtr+36, sdot+6);
    copy(dqdtDotPtr, dqdtDotPtr+6, sdot+42);

    // printf("sdot: [%.4f %.4f %.4f %.4f %.4f %.4f]\n", sdot[0], sdot[1], sdot[2], sdot[3], sdot[4], sdot[5]);
    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================

/**
 *   @brief Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   @param t time at integration step (unused)
 *   @param s the 6-d state vector
 *   @param sdot the 6-d state derivative vector
 *   @param params points to additional integration parameters wrapped in an 
 *  <tt>tpat_bcr4bpr_sys_data</tt> data object.
 */
int bcr4bpr_simple_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    tpat_bcr4bpr_sys_data *sysData = static_cast<tpat_bcr4bpr_sys_data *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, *sysData, primPosData);
    tpat_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    copy(s, s+3, r_data);
    tpat_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    copy(s+3, s+6, v_data);
    tpat_matrix v(3,1,v_data);

    // Create relative position vectors between s/c and primaries
    tpat_matrix r_p1 = r - trans(primPos.getRow(0));
    tpat_matrix r_p2 = r - trans(primPos.getRow(1));
    tpat_matrix r_p3 = r - trans(primPos.getRow(2));
    double d1 = norm(r_p1);
    double d2 = norm(r_p2);
    double d3 = norm(r_p3);
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Create C-matrix
    double c[] = {0, 2*k, 0, -2*k, 0, 0, 0, 0, 0};
    tpat_matrix C(3,3,c);

    // truncated position vector used in EOMs
    double r_trunc_data[3] = {0};
    copy(s, s+2, r_trunc_data);
    tpat_matrix r_trunc(3,1,r_trunc_data);

    // Compute acceleration using matrix math
    tpat_matrix accel(3,1);
    accel = C*v + k*k*r_trunc - (1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2, 3) - 
            nu*r_p3/pow(d3, 3);

    // Save derivatives to output vector
    double *accelPtr = accel.getDataPtr();

    copy(s+3, s+6, sdot);
    copy(accelPtr, accelPtr+3, sdot+3);

    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================


//-----------------------------------------------------
//      CR3BP Utility Functions
//-----------------------------------------------------


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
void cr3bp_getUDDots(double mu, double x, double y, double z, double* ddots){
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

/**
 *  @brief Compute the Jacobi Constant for the CR3BP
 *
 *  @param s the state vector; only the position and velocity states are required
 *  @param mu the non-dimensional system mass ratio
 *
 *  @return the Jacobi Constant at this specific state and system
 */
double cr3bp_getJacobi(double s[], double mu){
    double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);
    return 2*U - v*v;
}//================================================


/**
 *  @brief Compute the location of a Lagrange point in the CR3BP
 *
 *  @param sysData object describing the particular CR3BP
 *  @param L the Lagrange point number, 1 to 5
 *  @param tol the tolerance to use; if NAN is input, then a default value of 1e-14 will
 *  be used.
 *  @param pos a 3-element array to store the position of the Lagrange point
 */
void cr3bp_getEquilibPt(tpat_cr3bp_sys_data sysData, int L, double tol, double pos[3]){
    if(L < 1 || L > 5){
        throw tpat_exception("Invalid Lagrange Point");
    }

    if(tol == NAN){
        tol = 1e-14;
    }

    double mu = sysData.getMu();
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
            while(abs(gamma - gamma_prev) > tol && count < maxCount){   // Newton-Raphson for L1
                gamma_prev = gamma;
                gamma = gamma - ( mu/(gamma*gamma) - (1-mu)/pow(1-gamma, 2) - gamma - mu + 1)/
                    ( -2*mu/pow(gamma,3) - 2*(1-mu)/pow(1-gamma,3) - 1 );
                count++;
            }
            pos[0] = 1 - mu - gamma;
            break;
        case 2:
            gamma = 0.1;    // Initial guess is 10% of orbital radius
            while(abs(gamma - gamma_prev) > tol && count < maxCount){
                gamma_prev = gamma;
                gamma = gamma - ( -1*mu/(gamma*gamma) - (1-mu)/pow(1+gamma, 2) - mu + 1 + gamma)/
                    ( 2*mu/pow(gamma, 3) + 2*(1-mu)/pow(1+gamma, 3) + 1 );
                count++;
            }
            pos[0] = 1 - mu + gamma;
            break;
        case 3:
            gamma = 1;  // Initial guess is 100% of orbital radius
            while(abs(gamma - gamma_prev) > tol && count < maxCount){
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

    if(L < 4 && abs(gamma - gamma_prev) > tol){
        printErr("Could not converge on L%d\n", L);
    }
}//========================================

/**
 *  @brief Transition a trajectory from the EM system to the SE system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, so if the start of the trajectory does not coincide
 *  with this initial time, be sure to shift the time coordinates to assure that the first
 *  value in the time vector reflects correct initial time.
 *
 *  @param EMTraj a CR3BP Earth-Moon trajectory
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_cr3bp_traj cr3bp_EM2SE(tpat_cr3bp_traj EMTraj, double thetaE0, double thetaM0, double gamma){
    // Create a trajectory in the Sun-Earth system
    tpat_cr3bp_sys_data SESys("sun", "earth");
    tpat_cr3bp_traj SETraj(SESys);
    vector<double>* state = SETraj.getState();

    double charTE = EMTraj.getSysData().getCharT();     // characteristic time in EM system
    double charLE = EMTraj.getSysData().getCharL();     // characteristic length in EM system
    double charTS = SESys.getCharT();                   // characteristic time in SE system
    double charLS = SESys.getCharL();                   // characteristic length in SE system

    for(int n = 0; n < ((int)EMTraj.getTime()->size()); n++){

        // Transform the state from EM coordinates to SE coordinates
        vector<double> state_SE = cr3bp_EM2SE_state(EMTraj.getState(n), EMTraj.getTime(n), thetaE0,
            thetaM0, gamma, charLE, charTE, charLS, charTS, SESys.getMu());

        // Recompute Jacobi
        SETraj.getJC()->push_back(cr3bp_getJacobi(&(state_SE[0]), SESys.getMu()));

        // Copy transformed state into new vector
        state->insert(state->end(), state_SE.begin(), state_SE.end());

        // Re-Scale Time
        SETraj.getTime()->push_back(EMTraj.getTime(n)*charTE/charTS);

        // Put in a bogus STM
        SETraj.getSTM()->push_back(tpat_matrix::I(6));
    }

    SETraj.setLength();
    return SETraj; 
}//=========================================================

/**
 *  @brief Transition a nodeset from the EM system to the SE system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, but you can use the <tt>t0</tt> argument to specify
 *  an initial time since CR3BP nodesets do not have an epoch for each node.
 *
 *  @param EMNodes a CR3BP Earth-Moon nodeset
 *  @param t0 epoch associated with the first node
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_cr3bp_nodeset cr3bp_EM2SE(tpat_cr3bp_nodeset EMNodes, double t0, double thetaE0, double thetaM0,
    double gamma){

    tpat_cr3bp_sys_data SESys("sun", "earth");
    tpat_cr3bp_nodeset SENodes(SESys);
    vector<double>* nodes = SENodes.getNodes();

    double charTE = EMNodes.getSysData()->getCharT();       // characteristic time in EM system
    double charLE = EMNodes.getSysData()->getCharL();       // characteristic length in EM system
    double charTS = SESys.getCharT();                       // characteristic time in SE system
    double charLS = SESys.getCharL();                       // characteristic length in SE system

    double epoch = t0;
    printColor(BLUE, "Converting EM to SE\nEM Sys:\n  %d Nodes\n  %d TOFs\n", EMNodes.getNumNodes(),
        ((int)EMNodes.getTOFs()->size()));
    for(int n = 0; n < EMNodes.getNumNodes(); n++){
        vector<double> node_SE = cr3bp_EM2SE_state(EMNodes.getNode(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys.getMu());

        // Copy first 6 elements of transformed state/node into new node vector
        nodes->insert(nodes->end(), node_SE.begin(), node_SE.begin()+6);

        if(n < EMNodes.getNumNodes()-1){
            // Re-Scale Time
            SENodes.getTOFs()->push_back(EMNodes.getTOF(n)*charTE/charTS);

            // Update epoch time
            epoch += EMNodes.getTOF(n);
        }
    }

    return SENodes;
}//=========================================================

/**
 *  @brief Transition a trajectory from the SE system to the EM system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, so if the start of the trajectory does not coincide
 *  with this initial time, be sure to shift the time coordinates to assure that the first
 *  value in the time vector reflects correct initial time.
 *
 *  @param SETraj a CR3BP Sun-Earth trajectory
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_cr3bp_traj cr3bp_SE2EM(tpat_cr3bp_traj SETraj, double thetaE0, double thetaM0, double gamma){
    // Create a trajectory in the Earth-Moon system
    tpat_cr3bp_sys_data EMSys("earth", "moon");
    tpat_cr3bp_traj EMTraj(EMSys);
    vector<double>* state = EMTraj.getState();

    // Shift coordinates to EM barcyenter from SE barycenter
    tpat_matrix posShift = tpat_matrix::e_j(3, 1)*(1 - SETraj.getSysData().getMu());

    double charTE = EMSys.getCharT();               // characteristic time in EM system
    double charLE = EMSys.getCharL();               // characteristic length in EM system
    double charTS = SETraj.getSysData().getCharT(); // characteristic time in SE system
    double charLS = SETraj.getSysData().getCharL(); // characteristic length in SE system

    for(int n = 0; n < ((int)SETraj.getTime()->size()); n++){
        
        // Transform the state from SE coordinates to EM coordinates
        vector<double> state_EM = cr3bp_SE2EM_state(SETraj.getState(n), SETraj.getTime(n), thetaE0,
            thetaM0, gamma, charLE, charTE, charLS, charTS, SETraj.getSysData().getMu());

        // Recompute Jacobi
        EMTraj.getJC()->push_back(cr3bp_getJacobi(&(state_EM[0]), EMSys.getMu()));

        // Copy transformed state into new vector
        state->insert(state->end(), state_EM.begin(), state_EM.end());

        // Re-Scale Time
        EMTraj.getTime()->push_back(SETraj.getTime(n)*charTS/charTE);

        // Put in a bogus STM
        EMTraj.getSTM()->push_back(tpat_matrix::I(6));
    }

    EMTraj.setLength();
    return EMTraj; 
}//=========================================================

/**
 *  @brief Transition a nodeset from the SE system to the EM system
 *  
 *  The relative orientation between the two systems is described by the three angles
 *  <tt>thetaE0</tt>, <tt>thetaM0</tt>, and <tt>gamma</tt>. These angles describe the
 *  system orientation at time t = 0, but you can use the <tt>t0</tt> argument to specify
 *  an initial time since CR3BP nodesets do not have an epoch for each node.
 *
 *  @param SENodes a CR3BP Sun-Earth nodeset
 *  @param t0 epoch associated with the first node
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_cr3bp_nodeset cr3bp_SE2EM(tpat_cr3bp_nodeset SENodes, double t0, double thetaE0, double thetaM0,
    double gamma){

    tpat_cr3bp_sys_data EMSys("earth", "moon");
    tpat_cr3bp_sys_data SESys("sun", "earth");
    tpat_cr3bp_nodeset EMNodes(EMSys);
    vector<double>* nodes = EMNodes.getNodes();

    double charTE = EMSys.getCharT();                   // characteristic time in EM system
    double charLE = EMSys.getCharL();                   // characteristic length in EM system
    double charTS = SENodes.getSysData()->getCharT();    // characteristic time in SE system
    double charLS = SENodes.getSysData()->getCharL();    // characteristic length in SE system

    double epoch = t0;
    printColor(BLUE, "Converting SE to EM\nSE Sys:\n  %d Nodes\n  %d TOFs\n", SENodes.getNumNodes(),
        ((int)SENodes.getTOFs()->size()));
    for(int n = 0; n < SENodes.getNumNodes(); n++){
        // Transform a single node
        vector<double> node_EM = cr3bp_SE2EM_state(SENodes.getNode(n), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys.getMu());

        // Copy first 6 elements of transformed state/node into new node vector
        nodes->insert(nodes->end(), node_EM.begin(), node_EM.begin()+6);

        if(n < SENodes.getNumNodes()-1){
            // Re-Scale Time
            EMNodes.getTOFs()->push_back(SENodes.getTOF(n)*charTS/charTE);

            // Update epoch time
            epoch += SENodes.getTOF(n);
        }
    }

    return EMNodes;
}//=========================================================

/**
 *  @brief Transform a single state from EM coordinates to SE coordinates
 *
 *  @param state_EM a 6- or 9-element state vector
 *  @param t non-dimensional time associated with the state
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  @param charLE EM characteristic length
 *  @param charTE EM characterstic time
 *  @param charLS SE characteristic length
 *  @param charTS SE characteristic time
 *  @param mu_SE SE mass ratio
 *
 *  @return a 9-element state vector in EM coordinates
 */
vector<double> cr3bp_EM2SE_state(vector<double> state_EM, double t, double thetaE0, double thetaM0,
    double gamma, double charLE, double charTE, double charLS, double charTS, double mu_SE){

    // Shift coordinates to SE barcyenter from EM barycenter
    tpat_matrix posShift = tpat_matrix::e_j(3, 1)*(1 - mu_SE);

    // Compute Earth's position in inertial frame at this time
    double thetaE_k = thetaE0 + t*charTE/charTS;
    double thetaM_k = thetaM0 + t;

    // Compute DCMs
    double inertToLunar[] = {cos(gamma), 0, sin(gamma), 0, 1, 0, -sin(gamma), 0, cos(gamma)};
    double inertToSE[] = {cos(thetaE_k), -sin(thetaE_k), 0, sin(thetaE_k), cos(thetaE_k), 0,
                                0, 0, 1};
    double lunarToEM[] = {cos(thetaM_k), -sin(thetaM_k), 0, sin(thetaM_k), cos(thetaM_k), 0,
                                0, 0, 1};
    tpat_matrix DCM_I2L(3,3, inertToLunar);
    tpat_matrix DCM_I2S(3,3, inertToSE);
    tpat_matrix DCM_L2E(3,3, lunarToEM);

    tpat_matrix posEM(1,3, &(state_EM[0]));
    tpat_matrix velEM(1,3, &(state_EM[0])+3);

    // Rotate the position into SE frameand shift basepoint to SE barycenter (EM working frame)
    tpat_matrix posSE = posEM*trans(DCM_L2E)*trans(DCM_I2L)*DCM_I2S + posShift*charLS/charLE;

    // Angular velocity of SE frame in EM frame (working frame = EM)
    tpat_matrix omega = -1*charTE/charTS*tpat_matrix::e_j(3,3) *
        trans(DCM_I2S)*DCM_I2L*DCM_L2E + charTE/charTE*tpat_matrix::e_j(3,3);

    // Apply BKE to get velocity with SE observer, still in EM working frame
    tpat_matrix velSE = velEM + cross(omega, posEM);

    // Rotate the velocity into the SE working frame
    velSE *= trans(DCM_L2E)*trans(DCM_I2L)*DCM_I2S;

    // Units are still in non-dim EM, so change to SE non-dim
    posSE *= charLE/charLS;
    velSE *= (charLE/charTE)/(charLS/charTS);

    // Put new data into state vector
    vector<double> state_SE(3,0);   // three zeros for acceleration; insert pos and vel in front
    state_SE.insert(state_SE.begin(), posSE.getDataPtr(), posSE.getDataPtr()+3);
    state_SE.insert(state_SE.begin()+3, velSE.getDataPtr(), velSE.getDataPtr()+3);

    return state_SE;
}

/**
 *  @brief Transform a single state from SE coordinates to EM coordinates
 *
 *  @param state_SE a 6- or 9-element state vector
 *  @param t non-dimensional time associated with the state
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 *  @param charLE EM characteristic length
 *  @param charTE EM characterstic time
 *  @param charLS SE characteristic length
 *  @param charTS SE characteristic time
 *  @param mu_SE SE mass ratio
 *
 *  @return a 9-element state vector in SE coordinates
 */
vector<double> cr3bp_SE2EM_state(vector<double> state_SE, double t, double thetaE0, double thetaM0,
    double gamma, double charLE, double charTE, double charLS, double charTS, double mu_SE){
    
    tpat_matrix posShift = tpat_matrix::e_j(3, 1)*(1 - mu_SE);

    // Compute Earth's position in inertial frame at this time
    double thetaE_k = thetaE0 + t;
    double thetaM_k = thetaM0 + t*charTS/charTE;

    // Compute DCMs
    double inertToLunar[] = {cos(gamma), 0, sin(gamma), 0, 1, 0, -sin(gamma), 0, cos(gamma)};
    double inertToSE[] = {cos(thetaE_k), -sin(thetaE_k), 0, sin(thetaE_k), cos(thetaE_k), 0,
                                0, 0, 1};
    double lunarToEM[] = {cos(thetaM_k), -sin(thetaM_k), 0, sin(thetaM_k), cos(thetaM_k), 0,
                                0, 0, 1};
    tpat_matrix DCM_I2L(3,3, inertToLunar);
    tpat_matrix DCM_I2S(3,3, inertToSE);
    tpat_matrix DCM_L2E(3,3, lunarToEM);

    tpat_matrix posSE(1,3, &(state_SE[0]));
    tpat_matrix velSE(1,3, &(state_SE[0])+3);

    // Rotate the position into SE frame, coordinates are EM ND
    tpat_matrix posEM = (posSE - posShift)*trans(DCM_I2S)*DCM_I2L*DCM_L2E;

    // Angular velocity of EM frame in SE frame (working frame = SE)
    tpat_matrix omega = charTS/charTS*tpat_matrix::e_j(3,3) - charTS/charTE*tpat_matrix::e_j(3,3)*
        trans(DCM_I2L)*DCM_I2S;

    // Compute velocity in EM frame (working frame = SE)
    tpat_matrix velEM = velSE + cross(omega, (posSE-posShift));

    // Rotate the velocity into the EM working frame
    velEM *= trans(DCM_I2S)*DCM_I2L*DCM_L2E;

    // Units are still in non-dim SE, so change to EM non-dim
    posEM *= charLS/charLE;
    velEM *= (charLS/charTS)/(charLE/charTE);

    // Put new data into state vector
    vector<double> state_EM(3, 0);  // three zeros for acceleration; insert pos and vel in front
    state_EM.insert(state_EM.begin(), posEM.getDataPtr(), posEM.getDataPtr()+3);
    state_EM.insert(state_EM.begin()+3, velEM.getDataPtr(), velEM.getDataPtr()+3);

    return state_EM;
}

//-----------------------------------------------------
//      BCR4BP Utility Functions
//-----------------------------------------------------


/**
 *  @brief Compute the location of the three primaries in the BCR4BP (rotating coord.)
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param primPos a pointer to a 1x9 double array that will hold the positions of the three primaries in 
 *  row-major order. The first three elements are the position of P1, etc.
 */
void bcr4bpr_getPrimaryPos(double t, tpat_bcr4bpr_sys_data sysData, double *primPos){
    double k = sysData.getK();
    double mu = sysData.getMu();
    double nu = sysData.getNu();
    double theta0 = sysData.getTheta0();
    double phi0 = sysData.getPhi0();
    double gamma = sysData.getGamma();

    // Compute the angles for the system at the specified time
    double theta = theta0 + k*t;
    double phi = phi0 + sqrt(mu/pow(sysData.getCharLRatio(), 3)) * t;

    // P1 position
    primPos[0] = -mu;
    primPos[1] = 0;
    primPos[2] = 0;

    // P2 position
    primPos[3] = 1/k - mu - nu/mu*sysData.getCharLRatio() * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[4] = -nu/mu*sysData.getCharLRatio() * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[5] = nu/mu*sysData.getCharLRatio() * cos(phi) * sin(gamma);

    // P3 position
    primPos[6] = 1/k - mu + (1 - nu/mu)*sysData.getCharLRatio() * 
        (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[7] = (1 - nu/mu)*sysData.getCharLRatio() * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[8] = (nu/mu - 1)*sysData.getCharLRatio() * cos(phi)*sin(gamma);
}//================================================

/**
 *  @brief Compute the velocity of the three primaries in the BCR4BP, rotating coordinates.
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param primVel a pointer to a 3x3 double array that will hold the velocities of the three primaries in
 *  row-major order. The first three elements hold the velocity of P1, etc.
 */
void bcr4bpr_getPrimaryVel(double t, tpat_bcr4bpr_sys_data sysData, double *primVel){

    double k = sysData.getK();
    double mu = sysData.getMu();
    double nu = sysData.getNu();
    double theta0 = sysData.getTheta0();
    double phi0 = sysData.getPhi0();
    double gamma = sysData.getGamma();

    double thetaDot = k;
    double phiDot = sqrt(mu/pow(sysData.getCharLRatio(), 3));

    double theta = theta0 + thetaDot*t;
    double phi = phi0 + phiDot * t;    

    // P1 is stationary in this coordinate system
    primVel[0] = 0;  primVel[1] = 0;  primVel[2] = 0;

    // Angular velocity of P2-P3 line
    double v_P2P3Line[3] = {0};
    v_P2P3Line[0] = thetaDot*(sin(phi)*cos(theta) - cos(phi)*sin(theta)*cos(gamma)) + 
        phiDot*(cos(phi)*sin(theta) - sin(phi)*cos(theta)*cos(gamma));
    v_P2P3Line[1] = (phiDot - thetaDot)*cos(phi - theta);
    v_P2P3Line[2] = phiDot*sin(phi)*sin(gamma);

    // Multiply by radii of P2 and P3 to get their velocities
    primVel[3] = v_P2P3Line[0] * (-nu/mu)*sysData.getCharLRatio();
    primVel[4] = v_P2P3Line[1] * (-nu/mu)*sysData.getCharLRatio();
    primVel[5] = v_P2P3Line[2] * (-nu/mu)*sysData.getCharLRatio();

    primVel[6] = v_P2P3Line[0] * (1-nu/mu)*sysData.getCharLRatio();
    primVel[7] = v_P2P3Line[1] * (1-nu/mu)*sysData.getCharLRatio();
    primVel[8] = v_P2P3Line[2] * (1-nu/mu)*sysData.getCharLRatio();
}//===================================================================

/**
 *  @brief Convert a CR3BP Sun-Earth trajectory into a BCR4BPR Sun-Earth-Moon trajectory
 *
 *  @param crTraj a CR3BP Sun-Earth trajectory
 *  @param bcSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  @param t0 the epoch at the initial point on the CR3BP Trajectory
 *
 *  @return a BCR4BPR Trajectory object
 */
tpat_bcr4bpr_traj bcr4bpr_SE2SEM(tpat_cr3bp_traj crTraj, tpat_bcr4bpr_sys_data bcSys, double t0){
    string sun("Sun");      string earth("Earth");      string moon("Moon");

    if(crTraj.getSysData().getPrimID(0) != 10 || crTraj.getSysData().getPrimID(1) != 399){
        throw tpat_exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(bcSys.getPrimID(0) != 10 || bcSys.getPrimID(1) != 399 || bcSys.getPrimID(2) != 301){
        throw tpat_exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Trajectory
    tpat_bcr4bpr_traj bcTraj(bcSys);

    vector<double>* state = bcTraj.getState();
    vector<double>* times = bcTraj.getTime();
    vector<double>* dqdT = bcTraj.get_dqdT();
    vector<tpat_matrix>* STMs = bcTraj.getSTM();

    double charL2 = crTraj.getSysData().getCharL();
    double charT2 = crTraj.getSysData().getCharT();
    double charL3 = bcSys.getCharL();
    double charT3 = bcSys.getCharT();

    vector<double> blank(6, NAN);

    for(int n = 0; n < crTraj.getLength(); n++){
        vector<double> crState = crTraj.getState(n);
        printColor(BLUE, "crState[%d]: %d/%d\n", ((int)crState.size()), n, crTraj.getLength());
        for(int r = 0; r < ((int)crState.size()); r++){
            if(r < 3)   // Convert position
                state->push_back(crState[r]*charL2/charL3);
            else if(r < 6)  // Convert velocity
                state->push_back(crState[r]*(charL2/charL3)*(charT3/charT2));
            else    // Convert acceleration
                state->push_back(crState[r]*(charL2/charL3)*(charT3/charT2)*(charT3/charT2));
        }
        
        times->push_back(crTraj.getTime(n)*charT2/charT3 + t0);

        // Put bogus values into dqdT and STMs
        dqdT->insert(dqdT->end(), blank.begin(), blank.end());
        STMs->push_back(tpat_matrix::I(6));
    }

    bcTraj.setLength();
    return bcTraj;
}//==================================================



