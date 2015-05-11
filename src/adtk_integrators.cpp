/**
*   This library contains functions to numerically integrate the equations of 
*   motion associated with various mathematical models.
*/

#include "adtk_integrators.hpp"
#include "adtk_matrix.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>


using namespace std;

/**
*   Integrate the equations of motion for the CR3BP
*   @param t the current time of the integration
*   @param s the 6-d state vector
*   @param sdot the 6-d state derivative vector
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
 *  Integrate the 6 state EOMs and 36 STM EOMs for the CR3BP
 *  This function uses GSL's explicity embedded Runge-Kutta Dormand-Prince (8,9) method
 *
 *  @param ic a 6-element initial state for the trajectory
 *  @param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  @param mu the non-dimensional mass ratio of the system begin integrated
 *  @param t_dim the dimension of t
 *
 *  @return a vector containing all 42 states with a column of time valeus appended to the end
 */
vector<double> adtk_cr3bp_integrate(double ic[], double t[], double mu, int t_dim){

    // Construct the full 42-element IC from the state ICs plus the STM ICs
    const int ic_dim = 42;
    double fullIC[ic_dim] = {0};
    copy(ic, ic+6, fullIC);
    fullIC[6] = 1;          // Make the identity matrix
    fullIC[13] = 1;
    fullIC[20] = 1;
    fullIC[27] = 1;
    fullIC[34] = 1;
    fullIC[41] = 1;

    int steps = 0;  // count number of integration steps
    vector<double> state;
    double *y = new double[ic_dim];
    double *intTimes = new double[t_dim];

    // Define the step type (or, which integrator are we using?)
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk8pd;

    // Pointer to time-stepping function
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(stepType, ic_dim);

    // Define a control that will keep the error in the state y within the specified tolerances
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new (absTol, relTol);

    // Allocate space for the integrated solution to evolve in
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(ic_dim);

    // Create a system to integrate; we don't include a Jacobian (NULL)
    gsl_odeiv2_system sys = {cr3bp_EOMs, NULL, ic_dim, &mu};

    // Set y equal to the initial state
    for(int i = 0; i < ic_dim; i++){
        y[i] = fullIC[i];
    }

    // copy t into intTimes
    for(int i = 0; i < t_dim; i++){
        intTimes[i] = t[i];
    }

    // put initial state vector [ic, t(0)] into the state vector
    for(int i = 0; i < ic_dim+1; i++){
        state.push_back(i < ic_dim ? y[i] : t[0]);
    }
    steps++;

    if(t_dim == 2){
        double t0 = t[0], tf = t[1];            // start and finish times for integration
        double dt = tf > t0 ? dtGuess : -dtGuess;   // step size (initial guess)
        int sgn = tf > t0 ? 1 : -1;

        while (sgn*t0 < sgn*tf){
            // apply integrator for one step; new state and time are stored in y and t0
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);

            if(status != GSL_SUCCESS){
                printf("GSL did not succeed; Ending integration\n");
                break;
            }

            // Put newly integrated state and time into state vector
            for (int i = 0; i < ic_dim + 1; i++){
                state.push_back(i < ic_dim ? y[i] : t0);
            }

            steps++;
        }
    }else{
        // Integrate each segment between the input times
        for (int j = 0; j < t_dim - 1; j++){
            // define start and end times; t0 will be updated by integrator
            double t0 = t[j], tf = t[j+1];
            double dt = tf > t0 ? dtGuess : -dtGuess;
            int sgn = tf > t0 ? 1 : -1;

            while(sgn*t0 < sgn*tf){
                int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
                if(status != GSL_SUCCESS){ break; }
            }

            // Add the newly integrated state and current time fo the state vector
            for (int i = 0; i < ic_dim + 1; i++){
                state.push_back(i < ic_dim ? y[i] : t0);
            }

            steps++;
        }
    }

    // Clean Up
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    delete(y);
    delete(intTimes);

    return state;   // return the array of states
}//END of cr3bp_integrate