/**
 *  @file tpat_calculations.cpp
 *
 *  @brief  Contains non-member calculation functions
 *
 *   This library contains functions to numerically integrate the equations of 
 *   motion associated with various mathematical models. It also contains 
 *   miscellaneous functions to compute other quantities, like Jacobi Constant
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
#include "tpat.hpp"

#include "tpat_calculations.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_constants.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj.hpp"
#include "tpat_traj_bcr4bpr.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_traj_step.hpp"
#include "tpat_utilities.hpp"

#include "cspice/SpiceUsr.h"
#include "gsl/gsl_linalg.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>

//-----------------------------------------------------
//      Equations of Motion
//-----------------------------------------------------

/**
 *  @brief Integrate the equations of motion for the CR3BP
 *  @param t the current time of the integration; not used for this system
 *  @param s the 42-d state vector
 *  @param sdot the 42-d state derivative vector
 *  @param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to a cr3bp system data object
 */
int cr3bp_EOMs(double t, const double s[], double sdot[], void *params){
    (void)t;

    // Extract mu from params
    tpat_sys_data_cr3bp *sysData = static_cast<tpat_sys_data_cr3bp *>(params);
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
    std::copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    tpat_matrix phi(6,6, stmElements);
    
    // Compute derivative of STM
    tpat_matrix phiDot = A*phi;
    double *phiDotData = phiDot.getDataPtr();

    // Copy the elements of phiDot into the derivative array
    std::copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 6-d state vector
 *  @param sdot the 6-d state derivative vector
 *  @param params points to a cr3bp system data object
 */
int cr3bp_simple_EOMs(double t, const double s[], double sdot[], void *params){
    (void)t;
    
    // Extract mu from params
    tpat_sys_data_cr3bp *sysData = static_cast<tpat_sys_data_cr3bp *>(params);
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
 *  @brief Integrate the equations of motion for the CR3BP LTVP
 *  @param t the current time of the integration
 *  @param s the 43-d state vector. The first 6 elements are position and velocity,
 *  the 7th is mass, and the final 36 are STM elements
 *  @param sdot the 43-d state derivative vector
 *  @param *params pointer to extra parameters required for integration. For this
 *  function, the pointer points to a cr3bp system data object
 */
int cr3bp_ltvp_EOMs(double t, const double s[], double sdot[], void *params){
    // Extract mu from params
    tpat_sys_data_cr3bp_ltvp *sysData = static_cast<tpat_sys_data_cr3bp_ltvp *>(params);
    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    
    double g0_nonDim = (G_GRAV_0/charL)*charT*charT;
    double m = sysData->getM0() - T/(Isp*g0_nonDim)*t;    // assumes t began at 0

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );
    double v = sqrt( (xdot - y)*(xdot - y) + (ydot + x)*(ydot + x) + zdot*zdot);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3) - T/(m*v)*(xdot - y);
    sdot[4] = -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3) - T/(m*v)*(ydot + x);
    sdot[5] = -(1-mu)*z/pow(d,3) - mu*z/pow(r,3) - T/(m*v)*zdot; 

    /*
     * Next step, compute STM
     */
    
    // Partials of x_ddot w.r.t. state variables
    double dxdx = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*pow((x + mu),2)/pow(d,5) + 
        3*mu*pow((x + mu - 1), 2)/pow(r,5) + (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    double dxdy = 3*(1-mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) +
        (T/m)*(v*v + (xdot - y)*(xdot - y))/pow(v,3);
    double dxdz = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5);
    double dxdxdot = -1*(T/m)*(v*v - (xdot - y)*(xdot - y))/pow(v,3);
    double dxdydot = (T/m)*(xdot - y)*(ydot + x)/pow(v,3);
    double dxdzdot = (T/m)*(xdot - y)*zdot/pow(v,3);

    // Partials of y_ddot w.r.t. state variables
    double dydx = 3*(1 - mu)*(x + mu)*y/pow(d,5) + 3*mu*(x + mu - 1)*y/pow(r,5) - 
        (T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    double dydy = 1 - (1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*y*y/pow(d,5) + 3*mu*y*y/pow(r,5) +
        +(T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    double dydz = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5);
    double dydxdot = -2 + (T/m)*(ydot + x)*(xdot - y)/pow(v,3);
    double dydydot = -1*(T/m)*(v*v - (ydot + x)*(ydot + x))/pow(v,3);
    double dydzdot = (T/m)*(ydot + x)*zdot/pow(v,3);

    // Partials of z_ddot w.r.t. state variables
    double dzdx = 3*(1-mu)*(x + mu)*z/pow(d,5) + 3*mu*(x + mu - 1)*z/pow(r,5) + (T/m)*(ydot + x)*zdot/pow(v,3);
    double dzdy = 3*(1-mu)*y*z/pow(d,5) + 3*mu*y*z/pow(r,5) - (T/m)*(xdot - y)*zdot/pow(v,3);
    double dzdz = -(1-mu)/pow(d,3) - mu/pow(r,3) + 3*(1-mu)*z*z/pow(d,5) + 3*mu*z*z/pow(r,5);
    double dzdxdot = (T/m)*zdot*(xdot - y)/pow(v,3);
    double dzdydot = (T/m)*zdot*(ydot + x)/pow(v,3);
    double dzdzdot = -1*(T/m)*(v*v - zdot*zdot)/pow(v,3);

    // Create A Matrix
    double a_data[] = { 0,    0,    0,    1,       0,       0,
                        0,    0,    0,    0,       1,       0,
                        0,    0,    0,    0,       0,       1,
                        dxdx, dxdy, dxdz, dxdxdot, dxdydot, dxdzdot,
                        dydx, dydy, dydz, dydxdot, dydydot, dydzdot,
                        dzdx, dzdy, dzdz, dzdxdot, dzdydot, dzdzdot};
    tpat_matrix A(6,6,a_data);


    // Copy the STM states into a sub-array
    double stmElements[36];
    std::copy(s+6, s+42, stmElements);

    // Turn sub-array into matrix object for math stuffs
    tpat_matrix phi(6,6, stmElements);
    
    // Compute derivative of STM
    tpat_matrix phiDot = A*phi;
    double *phiDotData = phiDot.getDataPtr();

    // Copy the elements of phiDot into the derivative array
    std::copy(phiDotData, phiDotData+36, sdot+6);

    return GSL_SUCCESS;
}//===============================================================

/**
 *  @brief Integrate the equations of motion for the CR3BP LTVP without the STM
 *  @param t time at integration step (unused)
 *  @param s the 7-d state vector
 *  @param sdot the 7-d state derivative vector
 *  @param params points to a cr3bp system data object
 */
int cr3bp_ltvp_simple_EOMs(double t, const double s[], double sdot[], void *params){
    tpat_sys_data_cr3bp_ltvp *sysData = static_cast<tpat_sys_data_cr3bp_ltvp *>(params);
    double mu = sysData->getMu();
    double T = sysData->getThrust();
    double Isp = sysData->getIsp();
    double charT = sysData->getCharT();
    double charL = sysData->getCharL();

    double x = s[0];    double y = s[1];    double z = s[2];
    double xdot = s[3]; double ydot = s[4]; double zdot = s[5];
    
    double g0_nonDim = G_GRAV_0/charL*charT*charT;
    double m = sysData->getM0() - T/(Isp*g0_nonDim)*t;    // assumes t began at 0

    // compute distance to primaries and velocity magnitude
    double d = sqrt( (x+mu)*(x+mu) + y*y + z*z );
    double r = sqrt( (x-1+mu)*(x-1+mu) + y*y + z*z );
    double v = sqrt( (xdot - y)*(xdot - y) + (ydot + x)*(ydot + x) + zdot*zdot);

    sdot[0] = s[3];
    sdot[1] = s[4];
    sdot[2] = s[5];

    sdot[3] = 2*ydot + x - (1-mu)*(x+mu)/pow(d,3) - mu*(x-1+mu)/pow(r,3) - T/(m*v)*(xdot - y);
    sdot[4] = -2*xdot + y - (1-mu) * y/pow(d,3) - mu*y/pow(r,3) - T/(m*v)*(ydot + x);
    sdot[5] = -(1-mu)*z/pow(d,3) - mu*z/pow(r,3) - T/(m*v)*zdot; 

    return GSL_SUCCESS;
}//=====================================================

/**
 *   @brief Integrate the equations of motion for the BCR4BP, rotating coordinates.
 *
 *   @param t time at integration step (unused)
 *   @param s the 48-d state vector
 *   @param sdot the 48-d state derivative vector
 *   @param params points to additional integration parameters wrapped in an 
 *  <tt>tpat_sys_data_bcr4bpr</tt> data object.
 */
int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    tpat_sys_data_bcr4bpr *sysData = static_cast<tpat_sys_data_bcr4bpr *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, sysData, primPosData);
    tpat_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    std::copy(s, s+3, r_data);
    tpat_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    std::copy(s+3, s+6, v_data);
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
    std::copy(s, s+2, r_trunc_data);
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
    std::copy(s+6, s+42, phiData);
    tpat_matrix Phi(6, 6, phiData);
    tpat_matrix PhiDot = A*Phi;

    // Compute partials of state w.r.t. primary positions; dont' compute partials
    // for P1 because its velocity is zero in the rotating frame
    double dfdr2[18] = {0};   double dfdr3[18] = {0};

    dfdr2[9] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);        //dxdx2
    dfdr2[10] = 3*r_p2.at(0)*r_p2.at(1)/pow(d2,5);                  //dxdy2
    dfdr2[11] = 3*r_p2.at(0)*r_p2.at(2)/pow(d2,5);                  //dxdz2
    dfdr2[13] = -1/pow(d2,3) + 3*pow(r_p2.at(1),2)/pow(d2,5);       //dydy2
    dfdr2[14] = 3*r_p2.at(1)*r_p2.at(2)/pow(d2,5);                  //dydz2
    dfdr2[17] = -1/pow(d2,3) + 3*pow(r_p2.at(2),2)/pow(d2,5);       //dzdz2

    dfdr2[12] = dfdr2[10];      // Fill in symmetric matrix
    dfdr2[15] = dfdr2[11];
    dfdr2[16] = dfdr2[14];

    dfdr3[9] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);        //dxdx3
    dfdr3[10] = 3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);                  //dxdy3
    dfdr3[11] = 3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);                  //dxdz3
    dfdr3[13] = -1/pow(d3,3) + 3*pow(r_p3.at(1),2)/pow(d3,5);       //dydy3
    dfdr3[14] = 3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);                  //dydz3
    dfdr3[17] = -1/pow(d3,3) + 3*pow(r_p3.at(2),2)/pow(d3,5);       //dzdz3

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
    std::copy(s+42,s+48, dqdT_data);
    tpat_matrix dqdT(6,1, dqdT_data);

    // Get the velocity of the primaries
    double primVelData[9] = {0};
    bcr4bpr_getPrimaryVel(t, sysData, primVelData);
    tpat_matrix primVel(3,3, primVelData);

    // Compute derivative of dqdT
    tpat_matrix dot_dqdT = A*dqdT + DfDr2*trans(primVel.getRow(1)) + DfDr3*trans(primVel.getRow(2));

    // Save derivatives to output vector
    double *accelPtr = accel.getDataPtr();
    double *phiDotPtr = PhiDot.getDataPtr();
    double *dqdtDotPtr = dot_dqdT.getDataPtr();

    std::copy(s+3, s+6, sdot);
    std::copy(accelPtr, accelPtr+3, sdot+3);
    std::copy(phiDotPtr, phiDotPtr+36, sdot+6);
    std::copy(dqdtDotPtr, dqdtDotPtr+6, sdot+42);

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
 *  <tt>tpat_sys_data_bcr4bpr</tt> data object.
 */
int bcr4bpr_simple_EOMs(double t, const double s[], double sdot[], void *params){
    // Dereference the eom data object
    tpat_sys_data_bcr4bpr *sysData = static_cast<tpat_sys_data_bcr4bpr *>(params);

    // Put the positions of the three primaries in a 3x3 matrix
    double primPosData[9] = {0};
    bcr4bpr_getPrimaryPos(t, sysData, primPosData);
    tpat_matrix primPos(3, 3, primPosData);

    // Put the position states into a 3-element column vector
    double r_data[3] = {0};
    std::copy(s, s+3, r_data);
    tpat_matrix r(3,1,r_data);

    // Put velocity states into a 3-element column vector
    double v_data[3] = {0};
    std::copy(s+3, s+6, v_data);
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
    std::copy(s, s+2, r_trunc_data);
    tpat_matrix r_trunc(3,1,r_trunc_data);

    // Compute acceleration using matrix math
    tpat_matrix accel(3,1);
    accel = C*v + k*k*r_trunc - (1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2, 3) - 
            nu*r_p3/pow(d3, 3);

    // Save derivatives to output vector
    double *accelPtr = accel.getDataPtr();

    std::copy(s+3, s+6, sdot);
    std::copy(accelPtr, accelPtr+3, sdot+3);

    return GSL_SUCCESS;
}//============== END OF BCR4BPR EOMs ======================

//-----------------------------------------------------
//      General Utility Functions
//-----------------------------------------------------

/**
 *  @brief convert a date to epoch time
 *  
 *  @param date a string representing the date. The string can be formatted in one
 *  of two ways. First, a Gregorian-style date: 'YYYY/MM/DD HH:II:SS' (UTC, 24-hour clock);
 *  The time 'HH:II:SS' can be ommited, and the function will assume the time is 0:0:00
 *  Second, a Julian date (UTC) can be input with the format 'jd #' where '#' represents
 *  the Julian date.
 *
 *  @return the J2000 epoch time, or number of seconds after Jan 1, 2000 at 0:0:00 UTC.
 *  
 */
double dateToEpochTime(const char *date){
    // Load time kernel
    char timeKernel[] = "/Users/andrew/Documents/Purdue/Astrodynamics_Research/Code/C++/libTPAT/share/naif0010.tls.pc";
    // char timeKernel[] = "/home/andrew/Projects/Astrodynamics_Research/Code/C++/libTPAT/share/naif0010.tls.pc";
    furnsh_c(timeKernel);

    if(failed_c()){
        char errMsg[26];
        getmsg_c("short", 25, errMsg);
        printErr("Spice Error: %s\n", errMsg);
        reset_c();  // reset error status
        throw tpat_exception("Could not load SPICE time kernel");
    }

    // Convert the date to epoch time
    double et = 0;
    str2et_c(date, &et);

    // Unload the kernel
    unload_c(timeKernel);
    if(failed_c()){
        char errMsg[26];
        getmsg_c("short", 25, errMsg);
        printErr("Spice Error: %s\n", errMsg);
        reset_c();  // reset error status
        throw tpat_exception("Could not unload SPICE time kernel");
    }

    return et;
}//==========================================

/**
 *  @brief use least squares to predict new values of variables in a continuation process
 *
 *  This function uses a 2nd-order polynomial fit to predict a set of
 *  dependent variables using past relationships between an independent
 *  variable and the dependent variables. The number of points considered
 *  in those past relationships can be adjusted to vary the "stiffness"
 *  of the fit.
 *
 *  Ocasionally the 2nd-order fit does not work well and we encounter a 
 *  singular (or very near singular) matrix. In this case, the algorithm
 *  will apply linear regression, which generally solves the problem and
 *  results in a safe inversion.
 *
 *  For all these inputs, the "state vector" must take the following form:
 *      [x, y, z, xdot, ydot, zdot, Period, Jacobi Constant]
 *
 *  @param indVarIx the index of the state variable to use as the indpendent variable
 *  @param nextInd The next value of the independent variable
 *  @param depVars a vector specifying the indices of the states that will be
 *  dependent variables. The algorithm will predict fugure values for these
 *  variables based on how they have changed with the independent variable.
 *  @param varHistory a vector representing an n x 8 matrix which contains
 *  information about previous states, period, and JC. n should be at least 
 *  3. If it is larger than MemSize, only the first set of MemSize rows will
 *  be used.
 *
 *  @return an 8-element vector with predictions for the dependent variables.
 *  If a particular variable has not been predicted, its place will be kept with
 *  a NAN value.
 *
 *  An example may make things more clear:
 *
 *  Say I am continuing a family and am using x as the natural parameter in the
 *  continuation. indVarIx would be 0 to represent x. We input the value of x
 *  for the next orbit in the family (nextInd) and specify which variables (from
 *  the 8-element "state") we would like to have predicted by least-squares.
 */
std::vector<double> familyCont_LS(int indVarIx, double nextInd, std::vector<int> depVars, std::vector<double> varHistory){
    const int STATE_SIZE = 8;
    const double EPS = 1e-14;

    // Form A and B matrices
    std::vector<double> A_data;
    std::vector<double> B_data;

    for(size_t n = 0; n < varHistory.size()/STATE_SIZE; n++){
        double d = varHistory[n*STATE_SIZE + indVarIx];
        A_data.push_back(d*d);  // ind. var^2
        A_data.push_back(d);    // ind. var^1
        A_data.push_back(1);    // ind. var^0

        for(size_t p = 0; p < depVars.size(); p++){
            // vector of dependent variables
            B_data.push_back(varHistory[n*STATE_SIZE + depVars[p]]);
        }
    }

    tpat_matrix A(varHistory.size()/STATE_SIZE, 3, A_data);
    tpat_matrix B(varHistory.size()/STATE_SIZE, depVars.size(), B_data);

    // Generate coefficient matrix; these are coefficients for second-order
    // polynomials in the new independent variable
    tpat_matrix G = trans(A)*A;    // also a Gramm matrix
    tpat_matrix Gcopy = G;
    gsl_matrix *V = gsl_matrix_alloc(Gcopy.getRows(), Gcopy.getCols());
    gsl_vector *S = gsl_vector_alloc(Gcopy.getRows());
    gsl_vector *work = gsl_vector_alloc(Gcopy.getRows());

    int status = gsl_linalg_SV_decomp (Gcopy.getGSLMat(), V, S, work);

    if(status){
        printErr("tpat_calculations::familyCont_LS: GSL ERR: %s\n", gsl_strerror(status));
        throw tpat_linalg_err("Unable to take singular value decomposition");
    }


    double *smallestVal = std::min_element(S->data, S->data + S->size);
    tpat_matrix P(1,depVars.size());

    if(*smallestVal > EPS){
        // Use 2nd-order polynomial fit
        tpat_matrix C = solveAX_eq_B(G, trans(A)*B);

        double indMatData[] = {nextInd*nextInd, nextInd, 1};
        tpat_matrix indMat(1,3, indMatData);

        P = indMat*C;
    }else{
        // User 1st-order polynomial fit
        std::vector<double> A_lin_data;
        for(size_t n = 0; n < varHistory.size()/STATE_SIZE; n++){
            A_lin_data.push_back(varHistory[n*STATE_SIZE + indVarIx]);
            A_lin_data.push_back(1);
        }

        tpat_matrix A_lin(varHistory.size()/STATE_SIZE, 2, A_lin_data);
        tpat_matrix G = trans(A_lin)*A_lin;
        tpat_matrix C = solveAX_eq_B(G, trans(A_lin)*B);

        double indMatData[] = {nextInd, 1};
        tpat_matrix indMat(1,2, indMatData);

        P = indMat*C;
    }

    // Insert NAN for states that have not been predicted
    std::vector<double> predicted;
    predicted.assign(STATE_SIZE, NAN);
    for(size_t i = 0; i < depVars.size(); i++)
        predicted[depVars[i]] = P.at(i);

    return predicted;
}//====================================================

/**
 *  @brief construct a matrix to mirror a 6-d state over the specified plane or axis
 *  @param mirrorType describes how to mirror a 6-d state
 *  @return a 6x6 matrix that will mirror a 6-d state over the specified plane or axis
 */
tpat_matrix getMirrorMat(mirror_t mirrorType){
    switch(mirrorType){
        case MIRROR_XZ:
        {
            double data[] = {1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0,
                            0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, -1};
            return tpat_matrix(6,6,data);
        }
        case MIRROR_X_AX_H:
        case MIRROR_X_AX_V:
        {
            double data[] = {1, 0, 0, 0, 0, 0,
                            0, -1, 0, 0, 0, 0,
                            0, 0, -1, 0, 0, 0,
                            0, 0, 0, -1, 0, 0,
                            0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 1};
            return tpat_matrix(6,6,data);
        }
        default:
            printErr("tpat_calculations::getMirrorMat: Mirror type is not implemented; returning identiy\n");
            return tpat_matrix::I(6);
    }
}//===================================================

/**
 *  @brief Solve a matrix equation by iteratively applying LU factorization
 *
 *  Consider the matrix equation AX = B where A is an n x n square matrix,
 *  X is an n x m matrix or vector, and B is an n x m matrix or vector. We 
 *  solve this system by factoring A and solving individually for each column
 *  of X.
 *
 *  @param A an m x n matrix
 *  @param B an m x n matrix
 *  @return an n x n matrix that solves the systsem AX = B
 */
tpat_matrix solveAX_eq_B(tpat_matrix A, tpat_matrix B){
    if(A.getCols() != A.getRows()){
        printErr("tpat_calculations::solveAX_eq_B: A must be square");
        throw tpat_sizeMismatch("Cannot solve system unless A is square; cannot invert non-square matrix");
    }

    // printf("A = \n"); A.print("%12.8f");
    // printf("B = \n"); B.print("%12.8f");

    gsl_permutation *perm = gsl_permutation_alloc(A.getCols());
    int permSign = 0;
    int status = gsl_linalg_LU_decomp(A.getGSLMat(), perm, &permSign);
    if(status){
        printErr("tpat_calculations::solveAX_eq_B: GSL ERR: %s\n", gsl_strerror(status));
        throw tpat_linalg_err("Cannot form LU Decomp");
    }

    std::vector<double> solvedCols;
    gsl_vector *x = gsl_vector_alloc(B.getRows());

    for(int c = 0; c < B.getCols(); c++){
        tpat_matrix col = B.getCol(c);
        // printf("Applying inverse to column %d:\n",c);
        // col.print();
        gsl_vector_view gsl_col = gsl_vector_view_array(col.getDataPtr(), col.getRows());

        status = gsl_linalg_LU_solve(A.getGSLMat(), perm, &(gsl_col.vector), x);
        if(status){
            printErr("tpat_calculations::solveAX_eq_B: GSL ERR: %s\n", gsl_strerror(status));
            throw tpat_linalg_err("Could not solve for collumn of X");
        }
        tpat_matrix newCol(x, false);
        // printf("New column:\n");
        // newCol.print();
        solvedCols.insert(solvedCols.end(), x->data, x->data + x->size);
    }

    tpat_matrix allSolvedCols(B.getCols(), A.getCols(), solvedCols);
    return trans(allSolvedCols);
}//====================================================

/**
 *  @brief Sort a list of eigenvalues
 *
 *  This algorithm is intended to sort eigenvalues from the monodromy matrices of a family of
 *  periodic orbits. Because of this, the code attempts to locate at least one pair of 
 *  eigenvalues that are equal to 1.0. Other eigenvalues come in either real, reciprocal pairs
 *  or in complex conjugate pairs. The first row is sorted such that the pairs of eigenvalues
 *  are located next to each other in the row. All subsequent rows are sorted so that the order
 *  of the eigenvalues remains the same
 *
 *
 *  @param eigVals vector of eigenvalues in row-major order. We assume that 6 eigenvalues are 
 *  included in each row.
 *  @param sortedIxs a pointer to an integer vector that will store the indices of the original
 *  eigenvalues in their new order (after sorting)
 *  @return the sorted eigenvalues, again in row-major order in a vector
 */
std::vector<cdouble> sortEig(std::vector<cdouble> eigVals, std::vector<int> *sortedIxs){
    if(eigVals.size() == 0){
        printErr("tpat_calculations::sortEig: Cannot sort eigenvalues: there are no family members\n");
        return eigVals;
    }

    if(eigVals.size() % 6 != 0){
        printErr("tpat_calculations::sortEig: Must have 6n eigenvalues\n");
        return eigVals;
    }

    const double MAX_ONES_ERR = 1e-5;
    std::vector<cdouble> sortedEigs;

    // Figure out which eigenvalues are closest to one, save their indices
    std::vector<double> onesErr;
    std::vector<cdouble> e1(eigVals.begin(), eigVals.begin()+6);
    for(int i = 0; i < 6; i++){
        onesErr.push_back(std::abs(std::real(e1[i])-1) + std::abs(std::imag(e1[i])));
    }
    std::vector<double>::iterator smallestErr = std::min_element(onesErr.begin(), onesErr.end());
    int onesIx[] = {0,0};   // Indices of the eigenvalues that (nearly) exactly 1.0
    if(*smallestErr < MAX_ONES_ERR){
        int smallIx = smallestErr - onesErr.begin();
        if(smallIx % 2 == 0 && smallIx != 0){
            onesIx[0] = smallIx;
            onesIx[1] = smallIx+1;
        }else{
            onesIx[0] = smallIx-1;
            onesIx[1] = smallIx;
        }
    }else{
        printWarn("tpat_family_cr3bp::sortEigs: did not find eigenvalues at 1.0\n");
    }

    // Generate all permutations of the indices 0 through 5
    std::vector<int> vals {0,1,2,3,4,5};
    std::vector<int> ixPerms = tpat_util::generatePerms<int>(vals);
    cdouble predict[6];
    std::copy(e1.begin(), e1.end(), predict);

    for(size_t m = 0; m < eigVals.size()/6; m++){
        if(m == 1){
            if(*smallestErr < MAX_ONES_ERR){
                // Update the indices of the ones in case they got moved in the first sorting
                onesIx[0] = sortedIxs->at(onesIx[0]);
                onesIx[1] = sortedIxs->at(onesIx[1]);
            }
            std::copy(sortedEigs.begin(), sortedEigs.begin()+6, predict);
        }else if(m > 1){
            // Use linear extrapolation of m-1 and m-2 to predict m
            for(int i = 0; i < 6; i++){
                cdouble deltaEig = sortedEigs[(m-1)*6 + i] - sortedEigs[(m-2)*6 + i];
                predict[i] = sortedEigs[(m-1)*6 + i] + deltaEig;
            }
        }

        std::vector<double> cost;
        cost.assign(ixPerms.size()/6, 0);
        cdouble one(1,0);
        for(size_t p = 0; p < ixPerms.size()/6; p++){

            // rearrange the splitOrig vector using the current permutation
            cdouble swappedOrig[6];
            for(int i = 0; i < 6; i++){
                swappedOrig[i] = eigVals[m*6 + ixPerms[6*p + i]];
            }

            // Compute difference between this permutation and prediction
            for(int i = 0; i < 6; i++){
                cost[p] += std::abs(swappedOrig[i] - predict[i]);
            }
            
            // Add infinite cost if the ones eigenvalues change spots
            if(*smallestErr < MAX_ONES_ERR){
                double distToOne[6];
                for(int i = 0; i < 6; i++){
                    distToOne[i] = std::abs(swappedOrig[i]-one);
                }
                double *closestToOne = std::min_element(distToOne, distToOne+6);
                int closestIx = closestToOne - distToOne;

                if(closestIx != onesIx[0] && closestIx != onesIx[1]){
                    cost[p] += 1e20;
                    continue;
                }
            }

            // Add infinite cost if each consecutive pair are not inverses
            // or complex conjugates
            for(int i = 0; i < 6; i+=2){
                /* don't consider the eigenvalues that are closest to one
                    because they may have small innaccuracies that trigger
                    these costs, screwing up the algorithm */
                if(i != onesIx[0]){
                    // Complex parts don't sum to zero (conjugates will, two
                    // real eigenvalues will)
                    if(std::abs(std::imag(swappedOrig[i]) + std::imag(swappedOrig[i+1]))){
                        cost[p] += 1e20;
                        break;
                    }

                    // We have established that the ones are in the same place; Compute
                    // the error in a 1.0 eigenvalue; this error is acceptable in the
                    // reciprocal calculation (error tends to get bad at the larger orbits)
                    double okErr = std::abs(swappedOrig[onesIx[0]] - one);

                    // Both are real but are not reciprocals
                    if(std::imag(swappedOrig[i]) == 0 && std::imag(swappedOrig[i+1]) == 0){
                        double recipErr = 1.0 - std::real(swappedOrig[i])*std::real(swappedOrig[i+1]);
                        if(recipErr > okErr){
                            cost[p] += 1e20;
                            break;
                        }
                    }
                }
            }// end of checking consecutive pairs
        }// end of loop through permutations

        // Find the minimum cost
        std::vector<double>::iterator minCost = std::min_element(cost.begin(), cost.end());
        if(*minCost > 1e10){
            printErr("Minimum cost on set #%d is %f - probably a bug with the sorter!\n", m, *minCost);
        }
        int ix = minCost - cost.begin();
        for(int i = 0; i < 6; i++){
            sortedEigs.push_back(eigVals[m*6 + ixPerms[ix*6 + i]]);
        }

        sortedIxs->insert(sortedIxs->end(), ixPerms.begin() + ix*6, ixPerms.begin() + (ix+1)*6);
    }// end of loop through all members/eigenvalue sets

    return sortedEigs;
}//=====================================================

/**
 *
 *  Notes: No support for epoch time (yet)
 */
std::vector<tpat_traj_cr3bp> getManifolds(manifold_t type, tpat_traj_cr3bp *perOrbit, int numMans, double tof){
    // Get eigenvalues of monodromy matrix
    tpat_matrix mono = perOrbit->getSTM(-1);
    std::vector< std::vector<cdouble> > eigData = eig(mono);

    // Sort eigenvalues to put them in a "propper" order, get indices
    // to sort the eigenvectors to match
    std::vector<int> sortedIx;
    std::vector<cdouble> sortedEig = sortEig(eigData[0], &sortedIx);

    // Figure out which eigenvalues are the stable and unstable ones,
    // delete the rest
    std::vector<cdouble> nonCenterVals;
    std::vector<cdouble> nonCenterVecs;
    for(size_t c = 0; c < sortedEig.size(); c++){
        double realErr = std::real(sortedEig[c]) - 1.0;
        double imagErr = std::imag(sortedEig[c]);

        if(std::abs(realErr) > 1e-5 && std::abs(imagErr) < 1e-5){
            // Keep this eigenvalue/eigenvector pair
            nonCenterVals.push_back(sortedEig[c]);
            int vecIx = sortedIx[c];
            nonCenterVecs.insert(nonCenterVecs.end(), 
                eigData.at(1).begin()+vecIx*6,
                eigData.at(1).begin()+(vecIx+1)*6);
        }
    }

    std::vector<tpat_traj_cr3bp> allManifolds;

    if(nonCenterVals.size() == 0){
        printWarn("tpat_calculations::getManifolds: No stable/unstable eigenvalues were found\n");
        return allManifolds;
    }

    if(nonCenterVals.size() == 1){
        throw tpat_exception("tpat_calculations::getManifolds: Only found one stable/unstable eigenvalue");
    }

    if(nonCenterVals.size() > 2){
        printWarn("tpat_calculations::getManifolds: Stable/Unstable subspace is larger than 2D. Only the first pair will be considered\n");
    }

    /** TODO: If and when I can flexibily integrate generic trajectories and nodesets,
     *  I would like to replace this with a nodeset to space the points evenly
     *  in time and/or arclength
     */
    // Get a bunch of points to use as starting guesses for the manifolds
    if(numMans > perOrbit->getLength()){
        printWarn("tpat_calculations::getManifolds: Requested too many manifolds... will return fewer\n");
        numMans = perOrbit->getLength();
    }

    double stepSize = ((double)perOrbit->getLength())/((double)numMans);
    std::vector<int> pointIx(numMans, NAN);
    for(int i = 0; i < numMans; i++){
        pointIx[i] = floor(i*stepSize+0.5);
    }

    tpat_matrix temp(nonCenterVecs.size()/6, 6, tpat_util::real(nonCenterVecs));
    tpat_matrix vecs = trans(temp); // Transpose so eigenvectors are columns

    // NOW, copute the manifolds!
    tpat_simulation_engine sim(perOrbit->getSysData());
    double stepDist = 200;
    double charL = perOrbit->getSysData()->getCharL();
    for(int n = 0; n < numMans; n++){
        // Transform the eigenvectors to this updated time
        tpat_matrix newVecs = perOrbit->getSTM(pointIx[n])*vecs;

        // Pick the direction from one of the transformed eigenvectors
        tpat_matrix direction(6,1);
        for(size_t v = 0; v < 2; v++){
            tpat_matrix eigVec = newVecs.getCol(v);
            double mag = sqrt(eigVec.at(0)*eigVec.at(0) + 
                eigVec.at(1)*eigVec.at(1) + eigVec.at(2)*eigVec.at(2));
            if(type == MAN_U_P || type == MAN_U_M){
                if(std::abs(nonCenterVals[v]) > 1){
                    direction = eigVec/mag;
                    sim.setRevTime(false);
                    break;
                }
            }else{
                if(std::abs(nonCenterVals[v]) < 1){
                    direction = eigVec/mag;
                    sim.setRevTime(true);
                    break;
                }
            }
        }

        // Make sure it is pointing in +x direction
        direction *= tpat_util::sign(direction.at(0));
        
        // Orient according to specified type
        if(type == MAN_U_M || type == MAN_S_M)
            direction *= -1;

        // Step away from the point on the arc in the direction of the eigenvector
        tpat_matrix q0(6, 1, perOrbit->getState(pointIx[n]));
        q0 += stepDist/charL * direction;

        // Simulate for some time to generate a manifold arc
        sim.runSim(q0.getDataPtr(), tof);
        allManifolds.push_back(sim.getCR3BP_Traj());
    }

    return allManifolds;
}//====================================================

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
double cr3bp_getJacobi(const double s[], double mu){
    double v_squared = s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double d = sqrt((s[0] + mu)*(s[0] + mu) + s[1]*s[1] + s[2]*s[2]);
    double r = sqrt((s[0] - 1 + mu)*(s[0] - 1 + mu) + s[1]*s[1] + s[2]*s[2]);
    double U = (1 - mu)/d + mu/r + 0.5*(s[0]*s[0] + s[1]*s[1]);
    return 2*U - v_squared;
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
void cr3bp_getEquilibPt(tpat_sys_data_cr3bp sysData, int L, double tol, double pos[3]){
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
 *  @brief Compute a periodic orbit in the CR3BP system
 *
 *  The initial and final states are constrained based on the mirror type, but
 *  all other states are allowed to vary during the corrections process. If you
 *  wish to constrain a specific states, see 
 *  cr3bp_getPeriodic(sys, IC, period, numNodes, mirrorType, fixedStates)
 *  This function also uses only two nodes; to specify more, see the function above.
 *
 *  This function also assumes the order of the periodic orbit is 1
 *
 *  @param sys the dynamical system
 *  @param IC non-dimensional initial state vector
 *  @param period non-dimensional period for the orbit
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
tpat_traj_cr3bp cr3bp_getPeriodic(tpat_sys_data_cr3bp *sys, std::vector<double> IC,
    double period, mirror_t mirrorType){
    
    std::vector<int> fixedStates;   // Initialize an empty vector
    return cr3bp_getPeriodic(sys, IC, period, 2, 1, mirrorType, fixedStates);
}//========================================

/**
 *  @brief Compute a periodic orbit in the CR3BP system
 *
 *  @param sys the dynamical system
 *  @param IC non-dimensional initial state vector
 *  @param period non-dimensional period for the orbit
 *  @param numNodes the number of nodes to use; more nodes may result in a more robust correction
 *  @param order the number of revolutions about the system/primary this orbit completes before
 *  it repeats periodically. Think of a period-3 DRO (order = 3) or a butterfly (order = 2)
 *  @param mirrorType how this periodic orbit mirrors in the CR3BP
 *  @param fixedStates a vector containing the indices of which initial states
 *  we would like to fix. Not all states are possible for each mirror condition.
 *  See the enum definition for specific details.
 *  
 *  @return A periodic orbit. Note that this algorithm only enforces the mirror
 *  condition at the initial state and halfway point. To increase the accuracy
 *  of the periodic orbit, run it through a corrector to force the final state 
 *  to equal the first
 */
tpat_traj_cr3bp cr3bp_getPeriodic(tpat_sys_data_cr3bp *sys, std::vector<double> IC,
    double period, int numNodes, int order, mirror_t mirrorType, std::vector<int> fixedStates){

    tpat_simulation_engine sim(sys);    // Engine to perform simulation
    std::vector<int> zeroStates;            // Which states must be zero to ensure a perpendicular crossing

    tpat_event mirrorEvt;
    // Determine which states must be zero for mirroring
    switch(mirrorType){
        case MIRROR_XZ:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt = tpat_event(sys, tpat_event::XZ_PLANE, 0, true);  // Tell the sim to quit once it reaches the XZ plane
            break;
        case MIRROR_YZ:
            zeroStates.push_back(0);    // x
            zeroStates.push_back(4);    // y-dot
            zeroStates.push_back(5);    // z-dot
            mirrorEvt = tpat_event(sys, tpat_event::YZ_PLANE, 0, true);
            break;
        case MIRROR_XY:
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            zeroStates.push_back(4);    // y-dot
            mirrorEvt = tpat_event(sys, tpat_event::YZ_PLANE, 0, true);
            break;
        case MIRROR_X_AX_H:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt = tpat_event(sys, tpat_event::XZ_PLANE, 0, true);
            break;
        case MIRROR_X_AX_V:
            zeroStates.push_back(1);    // y
            zeroStates.push_back(2);    // z
            zeroStates.push_back(3);    // x-dot
            mirrorEvt = tpat_event(sys, tpat_event::XY_PLANE, 0, true);
            break;
        default:
            throw tpat_exception("Mirror type either not defined or not implemented");
    }

    mirrorEvt.setStopCount(order);
    sim.addEvent(mirrorEvt);

    // Create a constraint to enforce mirror condition at the beginning and end of arc
    double mirrorCon0[] = {NAN,NAN,NAN,NAN,NAN,NAN};
    double mirrorCon1[] = {NAN,NAN,NAN,NAN,NAN,NAN};

    // Set the zero states to zero in both constraints
    for(size_t i = 0; i < zeroStates.size(); i++){
        mirrorCon0[zeroStates[i]] = 0;
        mirrorCon1[zeroStates[i]] = 0;
    }

    // Fix states at the initial point
    for(size_t i = 0; i < fixedStates.size(); i++){
        bool okToFix = true;
        for(size_t n = 0; n < zeroStates.size(); n++){
            if(fixedStates[i] == zeroStates[n]){
                printWarn("Cannot fix state %d; it must be zero for this mirror condition; ignoring\n", fixedStates[i]);
                okToFix = false;
                break;
            }
        }
        if(okToFix){
            mirrorCon0[fixedStates[i]] = IC[fixedStates[i]];
        }
    }

    tpat_constraint initStateCon(tpat_constraint::STATE, 0, mirrorCon0, 6);
    tpat_constraint finalStateCon(tpat_constraint::STATE, numNodes-1, mirrorCon1, 6);

    // Run the sim until the event is triggered
    sim.runSim(IC, period);
    tpat_traj_cr3bp halfOrbArc = sim.getCR3BP_Traj();

    double halfOrbTOF = halfOrbArc.getTime(-1);
    double tofErr = 100*std::abs(halfOrbTOF-period/2.0)/(period/2.0);

    if(tofErr > 10)
        printWarn("tpat_calculations::cr3bp_getPeriodic: Half-Period arc TOF varies from input half-period by more than 10%%\n");

    // Create a nodeset from arc
    tpat_nodeset_cr3bp halfOrbNodes(halfOrbArc, numNodes, tpat_nodeset::DISTRO_TIME);
    halfOrbNodes.addConstraint(initStateCon);
    halfOrbNodes.addConstraint(finalStateCon);

    halfOrbArc.saveToMat("HalfOrbArc.mat");
    halfOrbNodes.saveToMat("HalfOrbNodes.mat");

    // Use differential corrections to enforce the mirror conditions
    tpat_correction_engine corrector;
    try{
        corrector.correct(&halfOrbNodes);
        tpat_nodeset_cr3bp correctedHalfPer = corrector.getCR3BP_Output();

        // Grab the last node, change its TOF and re-append it to the set
        tpat_node lastNode = correctedHalfPer.getNode(-1);
        lastNode.setTOF(correctedHalfPer.getNode(-2).getTOF());
        correctedHalfPer.deleteNode(-1);
        correctedHalfPer.appendNode(lastNode);

        // Now add more nodes, use MS to get an accurate orbit
        tpat_matrix mirrorMat = getMirrorMat(mirrorType);
        for(int i = correctedHalfPer.getNumNodes()-2; i >= 0; i--){
            // Use mirroring to populate second half of the orbit via nodes
            tpat_matrix stateVec(1,6, correctedHalfPer.getNode(i).getPosVelState());
            tpat_matrix newStateVec = stateVec*mirrorMat;

            double tof = NAN;
            if(i > 0)
                tof = correctedHalfPer.getNode(i-1).getTOF();

            tpat_node newNode(newStateVec.getDataPtr(), tof);
            correctedHalfPer.appendNode(newNode);
        }

        // Remove intermediate constraint, constrain both end points
        correctedHalfPer.clearConstraints();
        correctedHalfPer.addConstraint(initStateCon);
        tpat_constraint finalCon(tpat_constraint::STATE, correctedHalfPer.getNumNodes()-1, mirrorCon1, 6);
        correctedHalfPer.addConstraint(finalCon);

        correctedHalfPer.saveToMat("FullOrbNodes.mat");

        // Reconverge the solution
        corrector.correct(&correctedHalfPer);

        // Return the corrected solution in trajectory form
        tpat_nodeset_cr3bp finalSet = corrector.getCR3BP_Output();
        
        finalSet.saveToMat("FinalSetNodes.mat");
        waitForUser();
        
        return tpat_traj_cr3bp::fromNodeset(finalSet);
    }catch(tpat_diverge &e){
        throw tpat_diverge("tpat_calculations::cr3bp_getPeriodic: Could not converge half-period arc with mirroring condition");
    }
}//================================================

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
 *  @param SESys a Sun-Earth CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_traj_cr3bp cr3bp_EM2SE(tpat_traj_cr3bp EMTraj, tpat_sys_data_cr3bp *SESys, double thetaE0, double thetaM0, double gamma){
    // Create a trajectory in the Sun-Earth system
    tpat_traj_cr3bp SETraj(SESys);

    double charTE = EMTraj.getSysData()->getCharT();     // characteristic time in EM system
    double charLE = EMTraj.getSysData()->getCharL();     // characteristic length in EM system
    double charTS = SESys->getCharT();                   // characteristic time in SE system
    double charLS = SESys->getCharL();                   // characteristic length in SE system

    for(int n = 0; n < EMTraj.getLength(); n++){

        // Transform the state from EM coordinates to SE coordinates
        std::vector<double> state_SE = cr3bp_EM2SE_state(EMTraj.getState(n), EMTraj.getTime(n), thetaE0,
            thetaM0, gamma, charLE, charTE, charLS, charTS, SESys->getMu());

        // Adjust time to correct non-dim units
        double t = EMTraj.getTime(n)*charTE/charTS;

        // Create bogus accel and STM
        double accel[] = {NAN, NAN, NAN};
        tpat_matrix stm = tpat_matrix::I(6);

        // Create a new step
        tpat_traj_step step(&(state_SE[0]), t, accel, stm.getDataPtr());

        SETraj.appendStep(step);
        
        // Recompute Jacobi
        SETraj.setJacobi(-1, cr3bp_getJacobi(&(state_SE[0]), SESys->getMu()));
    }

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
 *  @param SESys a Sun-Earth CR3BP system data object
 *  @param t0 epoch associated with the first node
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_nodeset_cr3bp cr3bp_EM2SE(tpat_nodeset_cr3bp EMNodes, tpat_sys_data_cr3bp *SESys, double t0, double thetaE0, double thetaM0,
    double gamma){

    tpat_nodeset_cr3bp SENodes(SESys);

    double charTE = EMNodes.getSysData()->getCharT();       // characteristic time in EM system
    double charLE = EMNodes.getSysData()->getCharL();       // characteristic length in EM system
    double charTS = SESys->getCharT();                       // characteristic time in SE system
    double charLS = SESys->getCharL();                       // characteristic length in SE system

    double epoch = t0;
    printColor(BLUE, "Converting EM to SE\nEM Sys:\n  %d Nodes\n  %d TOFs\n", EMNodes.getNumNodes(),
        EMNodes.getNumNodes()-1);
    for(int n = 0; n < EMNodes.getNumNodes(); n++){
        std::vector<double> node_SE = cr3bp_EM2SE_state(EMNodes.getNode(n).getPosVelState(), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys->getMu());

        double tof = n + 1 < EMNodes.getNumNodes() ? EMNodes.getTOF(n)*charTE/charTS : NAN;
        tpat_node node(node_SE, tof);
        SENodes.appendNode(node);


        if(n < EMNodes.getNumNodes()-1){
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
 *  @param EMSys an Earth-Moon CR3BP system data object
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_traj_cr3bp cr3bp_SE2EM(tpat_traj_cr3bp SETraj, tpat_sys_data_cr3bp *EMSys, double thetaE0, double thetaM0, double gamma){
    // Create a trajectory in the Earth-Moon system
    tpat_traj_cr3bp EMTraj(EMSys);

    // Shift coordinates to EM barcyenter from SE barycenter
    tpat_sys_data_cr3bp *seSys = static_cast<tpat_sys_data_cr3bp*>(SETraj.getSysData());
    tpat_matrix posShift = tpat_matrix::e_j(3, 1)*(1 - seSys->getMu());

    double charTE = EMSys->getCharT();               // characteristic time in EM system
    double charLE = EMSys->getCharL();               // characteristic length in EM system
    double charTS = SETraj.getSysData()->getCharT(); // characteristic time in SE system
    double charLS = SETraj.getSysData()->getCharL(); // characteristic length in SE system

    for(int n = 0; n < SETraj.getLength(); n++){
        
        // Transform the state from SE coordinates to EM coordinates
        std::vector<double> state_EM = cr3bp_SE2EM_state(SETraj.getState(n), SETraj.getTime(n), thetaE0,
            thetaM0, gamma, charLE, charTE, charLS, charTS, seSys->getMu());
        
        // Rescale Time
        double t = SETraj.getTime(n)*charTS/charTE;

        // Bogus values for accel and STM
        double accel[] = {NAN, NAN, NAN};
        tpat_matrix stm = tpat_matrix::I(6);

        tpat_traj_step step(&(state_EM[0]), t, accel, stm.getDataPtr());
        EMTraj.appendStep(step);

        // Recompute Jacobi
        EMTraj.setJacobi(-1, cr3bp_getJacobi(&(state_EM[0]), EMSys->getMu()));
    }

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
 *  @param EMSys an Earth-Moon CR3BP system data object
 *  @param t0 epoch associated with the first node
 *  @param thetaE0 the angle (radians) between the Sun-Earth line and the 
 *  inertial x-axis at time t = 0.
 *  @param thetaM0 the angle (radians) between the Earth-Moon line and 
 *  lunar "periapse" at time t = 0.
 *  @param gamma the inclination (radians) of the lunar orbital plane relative 
 *  to the ecliptic; this value is held constant.
 */
tpat_nodeset_cr3bp cr3bp_SE2EM(tpat_nodeset_cr3bp SENodes, tpat_sys_data_cr3bp *EMSys, double t0, double thetaE0, double thetaM0,
    double gamma){

    tpat_nodeset_cr3bp EMNodes(EMSys);

    double charTE = EMSys->getCharT();                   // characteristic time in EM system
    double charLE = EMSys->getCharL();                   // characteristic length in EM system
    double charTS = SENodes.getSysData()->getCharT();    // characteristic time in SE system
    double charLS = SENodes.getSysData()->getCharL();    // characteristic length in SE system

    tpat_sys_data_cr3bp *SESys = static_cast<tpat_sys_data_cr3bp*>(SENodes.getSysData());

    double epoch = t0;
    printColor(BLUE, "Converting SE to EM\nSE Sys:\n  %d Nodes\n  %d TOFs\n", SENodes.getNumNodes(),
        SENodes.getNumNodes()-1);
    for(int n = 0; n < SENodes.getNumNodes(); n++){
        // Transform a single node
        std::vector<double> node_EM = cr3bp_SE2EM_state(SENodes.getNode(n).getPosVelState(), epoch, thetaE0, thetaM0,
            gamma, charLE, charTE, charLS, charTS, SESys->getMu());

        double tof = n + 1 < SENodes.getNumNodes() ? SENodes.getTOF(n)*charTS/charTE : NAN;
        tpat_node node(node_EM, tof);
        EMNodes.appendNode(node);


        if(n < SENodes.getNumNodes()-1){
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
 *  @return a 6-element state vector in EM coordinates
 */
std::vector<double> cr3bp_EM2SE_state(std::vector<double> state_EM, double t, double thetaE0, double thetaM0,
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
    std::vector<double> state_SE;
    state_SE.insert(state_SE.begin(), posSE.getDataPtr(), posSE.getDataPtr()+3);
    state_SE.insert(state_SE.begin()+3, velSE.getDataPtr(), velSE.getDataPtr()+3);

    return state_SE;
}//====================================================

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
 *  @return a 6-element state vector in SE coordinates
 */
std::vector<double> cr3bp_SE2EM_state(std::vector<double> state_SE, double t, double thetaE0, double thetaM0,
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
    std::vector<double> state_EM;
    state_EM.insert(state_EM.begin(), posEM.getDataPtr(), posEM.getDataPtr()+3);
    state_EM.insert(state_EM.begin()+3, velEM.getDataPtr(), velEM.getDataPtr()+3);

    return state_EM;
}//====================================================

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
void bcr4bpr_getPrimaryPos(double t, tpat_sys_data_bcr4bpr *sysData, double *primPos){
    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();
    double theta0 = sysData->getTheta0();
    double phi0 = sysData->getPhi0();
    double gamma = sysData->getGamma();

    // Compute the angles for the system at the specified time
    double theta = theta0 + k*t;
    double phi = phi0 + sqrt(mu/pow(sysData->getCharLRatio(), 3)) * t;

    // P1 position
    primPos[0] = -mu;
    primPos[1] = 0;
    primPos[2] = 0;

    // P2 position
    primPos[3] = 1/k - mu - nu/mu*sysData->getCharLRatio() * (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[4] = -nu/mu*sysData->getCharLRatio() * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[5] = nu/mu*sysData->getCharLRatio() * cos(phi) * sin(gamma);

    // P3 position
    primPos[6] = 1/k - mu + (1 - nu/mu)*sysData->getCharLRatio() * 
        (cos(phi)*cos(gamma)*cos(theta) + sin(phi)*sin(theta));
    primPos[7] = (1 - nu/mu)*sysData->getCharLRatio() * (sin(phi)*cos(theta) - cos(phi)*sin(theta));
    primPos[8] = (nu/mu - 1)*sysData->getCharLRatio() * cos(phi)*sin(gamma);
}//================================================

/**
 *  @brief Compute the velocity of the three primaries in the BCR4BP, rotating coordinates.
 *
 *  @param t non-dimensional time since t0, where t0 coincides with the positions specified by theta0 and phi0
 *  @param sysData a system data object containing information about the BCR4BP primaries
 *  @param primVel a pointer to a 3x3 double array that will hold the velocities of the three primaries in
 *  row-major order. The first three elements hold the velocity of P1, etc.
 */
void bcr4bpr_getPrimaryVel(double t, tpat_sys_data_bcr4bpr *sysData, double *primVel){

    double k = sysData->getK();
    double mu = sysData->getMu();
    double nu = sysData->getNu();
    double theta0 = sysData->getTheta0();
    double phi0 = sysData->getPhi0();
    double gamma = sysData->getGamma();

    double thetaDot = k;
    double phiDot = sqrt(mu/pow(sysData->getCharLRatio(), 3));

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
    primVel[3] = v_P2P3Line[0] * (-nu/mu)*sysData->getCharLRatio();
    primVel[4] = v_P2P3Line[1] * (-nu/mu)*sysData->getCharLRatio();
    primVel[5] = v_P2P3Line[2] * (-nu/mu)*sysData->getCharLRatio();

    primVel[6] = v_P2P3Line[0] * (1-nu/mu)*sysData->getCharLRatio();
    primVel[7] = v_P2P3Line[1] * (1-nu/mu)*sysData->getCharLRatio();
    primVel[8] = v_P2P3Line[2] * (1-nu/mu)*sysData->getCharLRatio();
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
tpat_traj_bcr4bpr bcr4bpr_SE2SEM(tpat_traj_cr3bp crTraj, tpat_sys_data_bcr4bpr *bcSys, double t0){
    if(crTraj.getSysData()->getPrimID(0) != 10 || crTraj.getSysData()->getPrimID(1) != 399){
        throw tpat_exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(bcSys->getPrimID(0) != 10 || bcSys->getPrimID(1) != 399 || bcSys->getPrimID(2) != 301){
        throw tpat_exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Trajectory
    tpat_traj_bcr4bpr bcTraj(bcSys);

    double charL2 = crTraj.getSysData()->getCharL();
    double charT2 = crTraj.getSysData()->getCharT();
    double charL3 = bcSys->getCharL();
    double charT3 = bcSys->getCharT();

    std::vector<double> blank(6, NAN);

    for(int n = 0; n < crTraj.getLength(); n++){
        std::vector<double> crState = crTraj.getState(n);

        double bcState[6];
        for(int r = 0; r < 6; r++){
            if(r < 3)   // Convert position
                bcState[r] = crState[r]*charL2/charL3;
            else if(r < 6)  // Convert velocity
                bcState[r] = crState[r]*(charL2/charL3)*(charT3/charT2);
        }

        // Convert acceleration
        std::vector<double> crAccel = crTraj.getAccel(n);
        double bcAccel[3];
        for(int r = 0; r < 3; r++){
            bcAccel[r] = crAccel[r]*(charL2/charL3)*(charT3/charT2)*(charT3/charT2);
        }
        
        double t = crTraj.getTime(n)*charT2/charT3 + t0;

        // Bogus values for dqdT and STM
        double dqdT[] = {NAN, NAN, NAN, NAN, NAN, NAN};
        tpat_matrix stm = tpat_matrix::I(6);

        // Create step
        tpat_traj_step step(bcState, t, bcAccel, stm.getDataPtr());
        bcTraj.appendStep(step);
        bcTraj.set_dqdT(-1, dqdT);
    }

    return bcTraj;
}//==================================================

/**
 *  @brief Convert a Sun-Earth CR3BP Nodeset to a Sun-Earth-Moon BCR4BPR Nodeset
 *
 *  @param crNodes a CR3BP Sun-Earth nodeset
 *  @param bcSys a BCR4BPR Sun-Earth-Moon system data object; contains information about system
 *  scaling and orientation at time t = 0
 *  @param t0 the epoch at the first node on the CR3BP nodeset
 *
 *  @return a BCR4BPR nodeset
 */
tpat_nodeset_bcr4bpr bcr4bpr_SE2SEM(tpat_nodeset_cr3bp crNodes, tpat_sys_data_bcr4bpr *bcSys, double t0){
    if(crNodes.getSysData()->getPrimID(0) != 10 || crNodes.getSysData()->getPrimID(1) != 399){
        throw tpat_exception("CR3BP trajectory is not in the Sun-Earth System");
    }

    if(bcSys->getPrimID(0) != 10 || bcSys->getPrimID(1) != 399 || bcSys->getPrimID(2) != 301){
        throw tpat_exception("BCR4BPR system is not Sun-Earth-Moon");
    }

    // Create a BCR4BPR Trajectory
    tpat_nodeset_bcr4bpr bcNodes(bcSys);

    double charL2 = crNodes.getSysData()->getCharL();
    double charT2 = crNodes.getSysData()->getCharT();
    double charL3 = bcSys->getCharL();
    double charT3 = bcSys->getCharT();

    double ellapsed = t0;

    for(int n = 0; n < crNodes.getNumNodes(); n++){
        std::vector<double> bcNodeState;
        std::vector<double> crNode = crNodes.getNode(n).getPosVelState();
        for(int r = 0; r < ((int)crNode.size()); r++){
            if(r < 3)   // Convert position
                bcNodeState.push_back(crNode[r]*charL2/charL3);
            else  // Convert velocity
                bcNodeState.push_back(crNode[r]*(charL2/charL3)*(charT3/charT2));
        }
        
        double tof = crNodes.getTOF(n)*charT2/charT3;
        tpat_node bcNode(bcNodeState, tof);
        bcNode.setExtraParam(1, ellapsed);  // set epoch

        bcNodes.appendNode(bcNode);
        if(n < crNodes.getNumNodes()-1){
            ellapsed += tof;
        }
    }

    return bcNodes;
}//==================================================

/**
 *  @brief Orient a BCR4BPR system so that t = 0 corresponds to the specified epoch time
 *  
 *  @param et epoch time (seconds, J2000, UTC)
 *  @param sysData pointer to system data object; A new theta0 and phi0 will be stored
 *  in this data object
 */
void bcr4bpr_orientAtEpoch(double et, tpat_sys_data_bcr4bpr *sysData){
    double time_nonDim = (et - tpat_sys_data_bcr4bpr::REF_EPOCH)/sysData->getCharT();
    
    // Compute theta and phi
    double theta = sysData->getK()*time_nonDim;
    double phi = sqrt(sysData->getMu()/pow(sysData->getCharLRatio(), 3))*time_nonDim;

    // Adjust theta and phi to be between 0 and 2*PI
    theta -= floor(theta/(2*PI))*2*PI;
    phi -= floor(phi/(2*PI))*2*PI;

    sysData->setTheta0(theta);
    sysData->setPhi0(phi);
}//===========================================





