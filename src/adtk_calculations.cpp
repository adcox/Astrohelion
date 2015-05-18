/**
 *  @file adtk_calculations.cpp
 *
 *   This library contains functions to numerically integrate the equations of 
 *   motion associated with various mathematical models. It also contains 
 *   miscellaneous functions to compute other quantities, like Jacobi Constant
 */

#include "adtk_calculations.hpp"

#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_matrix.hpp"

#include <math.h>


using namespace std;

//-----------------------------------------------------
//      Equations of Motion
//-----------------------------------------------------

/**
*   Integrate the equations of motion for the CR3BP
*   @param t the current time of the integration
*   @param s the 42-d state vector
*   @param sdot the 42-d state derivative vector
*   @param *params pointer to extra parameters required for integration. For this
*   function, the only extra parameter is mu
*/
int cr3bp_EOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params by casting it to a dobule
    double mu = *(double *)params;

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
    double *phiDotData = phiDot.getGSLMat()->data;

    // Copy the elements of phiDot into the derivative array
    copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *   Integrate the equations of motion for the CR3BP without the STM
 *   @param t time at integration step (unused)
 *   @param s the 6-d state vector
 *   @param sdot the 6-d state derivative vector
 *   @param params points to a double that contains the value of mu
 */
int cr3bp_simple_EOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params by casting it to a dobule
    double mu = *(double *)params;

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
 *   Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   @param t time at integration step (unused)
 *   @param s the 48-d state vector
 *   @param sdot the 48-d state derivative vector
 *   @param params points to additional integration parameters wrapped in an 
 *  <tt>adtk_bcr4bpr_eomData</tt> data object.
 */
int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    adtk_bcr4bpr_eomData eomData = *(adtk_bcr4bpr_eomData *)params;

    // Pull out the system data object
    adtk_bcr4bpr_sys_data sysData = eomData.sysData;

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, sysData, eomData.theta0, eomData.phi0, eomData.gamma, primPosData);
    adtk_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double tmp[3] = {0};
    copy(s, s+3, tmp);
    adtk_matrix r(3,1,tmp);

    // Put velocity states into a 3-element column vector
    copy(s+3, s+6, tmp);
    adtk_matrix v(3,1,tmp);

    // Create relative position vectors between s/c and primaries
    adtk_matrix r_p1 = r - primPos.getRow(0).trans();
    adtk_matrix r_p2 = r - primPos.getRow(1).trans();
    adtk_matrix r_p3 = r - primPos.getRow(2).trans();

    
    // Save constants to short variables for readability
    double k = sysData.getK();
    double mu = sysData.getMu();
    double nu = sysData.getNu();

    // Create C-matrix
    double c[] = {0, 2*k, 0, -2*k, 0, 0, 0, 0, 0};
    adtk_matrix C(3,3,c);

    // truncated position vector used in EOMs
    tmp[0] = s[0]; tmp[1] = s[1]; tmp[2] = 0;
    adtk_matrix r_trunc(3,1,tmp);

    adtk_matrix accel(3,1);
    accel = C*v + k*k*r_trunc - (1/k - mu)*r_p1/pow(r_p1.norm(), 3) - (mu - nu)*r_p2/pow(r_p2.norm(), 3) - 
            nu*r_p3/pow(r_p3.norm(), 3);

    // TODO Put in computations for these!
    double dxdx = 0;
    double dxdy = 0;
    double dxdz = 0;
    double dydy = 0;
    double dydz = 0;
    double dzdz = 0;

    double aData[] = {  0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 1, 
                        dxdx, dxdy, dxdz, c[0], c[1], c[2],
                        dxdy, dydy, dydz, c[3], c[4], c[5],
                        dxdz, dydz, dzdz, c[6], c[7], c[8]};

    adtk_matrix A(6,6, aData);
    
    double phiData[36];
    copy(s+6, s+42, phiData);
    adtk_matrix Phi(6,6, phiData);

    adtk_matrix PhiDot = A*Phi;

    // Save derivatives to output vector
    double *accelPtr = accel.getGSLMat()->data;
    double *phiDotPtr = PhiDot.getGSLMat()->data;

    copy(s+3, s+6, sdot);
    copy(accelPtr, accelPtr+3, sdot+3);
    copy(phiDotPtr, phiDotPtr+36, sdot+6);

    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================


//-----------------------------------------------------
//      CR3BP Utility Functions
//-----------------------------------------------------


/**
 *  Compute the second derivatives of the pseudo-potential function
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
 *  Compute the Jacobi Constant for the CR3BP
 *
 *  @param s the state vector; only the position and velocity states are required
 *  @para mu the non-dimensional system mass ratio
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
 *  Compute the location of the three primaries in the BCR4BP (rotating coord.)
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi9
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param theta0 the angle between the P1/P2 system and the inertial x-axis, in radians
 *  @param phi0 the angle between the P2/P3 system and inertial x-axis (projected into XY plane), in radians
 *  @param gamma inclination of the P2/P3 system relative to the P1/P2 plane, in radians. This inclination
 *  will remain constant.
 *  @param primPos a pointer to a 1x9 double array that will hold the positions of the three primaries in 
 *  row-major order. The first three elements are the position of P1, etc.
 */
void bcr4bpr_getPrimaryPos(double t, adtk_bcr4bpr_sys_data sysData, double theta0, double phi0, double gamma, 
        double *primPos){
    double k = sysData.getK();
    double mu = sysData.getMu();
    double nu = sysData.getNu();

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
 *  Compute the velocity of the three primaries in the BCR4BP, rotating coordinates.
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi9
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param theta0 the angle between the P1/P2 system and the inertial x-axis, in radians
 *  @param phi0 the angle between the P2/P3 system and inertial x-axis (projected into XY plane), in radians
 *  @param gamma inclination of the P2/P3 system relative to the P1/P2 plane, in radians. This inclination
 *  will remain constant.
 *  @param primVel a pointer to a 3x3 double array that will hold the velocities of the three primaries in
 *  row-major order. The first three elements hold the velocity of P1, etc.
 */
void bcr4bpr_getPrimaryVel(double t, adtk_bcr4bpr_sys_data sysData, double theta0, double phi0, double gamma,
    double *primVel){

    double k = sysData.getK();
    double mu = sysData.getMu();
    double nu = sysData.getNu();

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