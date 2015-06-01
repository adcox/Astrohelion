/**
 *  @file adtk_calculations.cpp
 *
 *   This library contains functions to numerically integrate the equations of 
 *   motion associated with various mathematical models. It also contains 
 *   miscellaneous functions to compute other quantities, like Jacobi Constant
 */

#include "adtk_calculations.hpp"

#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_matrix.hpp"

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
    adtk_cr3bp_sys_data *sysData = static_cast<adtk_cr3bp_sys_data *>(params);
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
    adtk_matrix A(6,6,a_data);

    // Copy the STM states into a sub-array
    double stmElements[36];
    copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    adtk_matrix phi(6,6, stmElements);
    
    // Compute derivative of STM
    adtk_matrix phiDot = A*phi;
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
    adtk_cr3bp_sys_data *sysData = static_cast<adtk_cr3bp_sys_data *>(params);
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
 *  <tt>adtk_bcr4bpr_sys_data</tt> data object.
 */
int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    adtk_bcr4bpr_sys_data *sysData = static_cast<adtk_bcr4bpr_sys_data *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, *sysData, primPosData);
    adtk_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    copy(s, s+3, r_data);
    adtk_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    copy(s+3, s+6, v_data);
    adtk_matrix v(3,1,v_data);

    // Create relative position vectors between s/c and primaries
    adtk_matrix r_p1 = r - primPos.getRow(0).trans();
    adtk_matrix r_p2 = r - primPos.getRow(1).trans();
    adtk_matrix r_p3 = r - primPos.getRow(2).trans();
    double d1 = r_p1.norm();
    double d2 = r_p2.norm();
    double d3 = r_p3.norm();
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Create C-matrix
    double c[] = {0, 2*k, 0, -2*k, 0, 0, 0, 0, 0};
    adtk_matrix C(3,3,c);

    // truncated position vector used in EOMs
    double r_trunc_data[3] = {0};
    copy(s, s+2, r_trunc_data);
    adtk_matrix r_trunc(3,1,r_trunc_data);

    // Compute acceleration using matrix math
    adtk_matrix accel(3,1);
    accel = C*v + k*k*r_trunc - (1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2, 3) - 
            nu*r_p3/pow(d3, 3);

    // Compute psuedo-potential
    double dxdx = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(0),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 3*pow(r_p3.at(0),2)/pow(d3,5));
    double dxdy = (1/k - mu)*3*r_p1.at(0)*r_p1.at(1)/pow(d1,5) + (mu - nu)*3*r_p2.at(0)*r_p2.at(1)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);
    double dxdz = (1/k - mu)*3*r_p1.at(0)*r_p1.at(2)/pow(d1,5) + (mu - nu)*3*r_p2.at(0)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);
    double dydy = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(1),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(1),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 3*pow(r_p3.at(1),2)/pow(d3,5));
    double dydz = (1/k - mu)*3*r_p1.at(1)*r_p1.at(2)/pow(d1,5) + (mu - nu)*3*r_p2.at(1)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);
    double dzdz = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(2),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 3*pow(r_p3.at(2),2)/pow(d3,5));

    // Create A matrix for STM derivative
    double aData[] = {  0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1, 
                        dxdx, dxdy, dxdz, c[0], c[1], c[2],
                        dxdy, dydy, dydz, c[3], c[4], c[5],
                        dxdz, dydz, dzdz, c[6], c[7], c[8]};
    adtk_matrix A(6,6, aData);
    
    // Compute the STM derivative
    double phiData[36];
    copy(s+6, s+42, phiData);
    adtk_matrix Phi(6, 6, phiData);
    adtk_matrix PhiDot = A*Phi;

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

    adtk_matrix DfDr2(6,3, dfdr2);
    adtk_matrix DfDr3(6,3, dfdr3);

    // Pull the state derivative w.r.t. Epoch time from the large state vector; create column vector
    double dqdT_data[6] = {0};
    copy(s+42,s+48, dqdT_data);
    adtk_matrix dqdT(6,1, dqdT_data);

    // Get the velocity of the primaries
    double primVelData[9] = {0};
    bcr4bpr_getPrimaryVel(t, *sysData, primVelData);
    adtk_matrix primVel(3,3, primVelData);

    // Compute derivative of dqdT
    adtk_matrix dot_dqdT = A*dqdT + DfDr2*(primVel.getRow(1).trans()) + DfDr3*(primVel.getRow(2).trans());

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
 *  <tt>adtk_bcr4bpr_sys_data</tt> data object.
 */
int bcr4bpr_simple_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    adtk_bcr4bpr_sys_data *sysData = static_cast<adtk_bcr4bpr_sys_data *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, *sysData, primPosData);
    adtk_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    copy(s, s+3, r_data);
    adtk_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    copy(s+3, s+6, v_data);
    adtk_matrix v(3,1,v_data);

    // Create relative position vectors between s/c and primaries
    adtk_matrix r_p1 = r - primPos.getRow(0).trans();
    adtk_matrix r_p2 = r - primPos.getRow(1).trans();
    adtk_matrix r_p3 = r - primPos.getRow(2).trans();
    double d1 = r_p1.norm();
    double d2 = r_p2.norm();
    double d3 = r_p3.norm();
    
    // Save constants to short variables for readability
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();

    // Create C-matrix
    double c[] = {0, 2*k, 0, -2*k, 0, 0, 0, 0, 0};
    adtk_matrix C(3,3,c);

    // truncated position vector used in EOMs
    double r_trunc_data[3] = {0};
    copy(s, s+2, r_trunc_data);
    adtk_matrix r_trunc(3,1,r_trunc_data);

    // Compute acceleration using matrix math
    adtk_matrix accel(3,1);
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


//-----------------------------------------------------
//      BCR4BP Utility Functions
//-----------------------------------------------------


/**
 *  @brief Compute the location of the three primaries in the BCR4BP (rotating coord.)
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi9
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param primPos a pointer to a 1x9 double array that will hold the positions of the three primaries in 
 *  row-major order. The first three elements are the position of P1, etc.
 */
void bcr4bpr_getPrimaryPos(double t, adtk_bcr4bpr_sys_data sysData, double *primPos){
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
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi9
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param primVel a pointer to a 3x3 double array that will hold the velocities of the three primaries in
 *  row-major order. The first three elements hold the velocity of P1, etc.
 */
void bcr4bpr_getPrimaryVel(double t, adtk_bcr4bpr_sys_data sysData, double *primVel){

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
}