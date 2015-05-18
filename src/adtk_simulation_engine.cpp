/**
 *  @file adtk_simulation_engine.cpp
 *
 * 	Simulation Engine
 *
 *	This object handles all simulation tasks.
 *
 */
#include "adtk_simulation_engine.hpp"

#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_bcr4bpr_traj.hpp"
#include "adtk_calculations.hpp"
#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_cr3bp_traj.hpp"
#include "adtk_matrix.hpp"

#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_simulation_engine::adtk_simulation_engine(){
	revTime = false;
	verbose = false;
    simpleIntegration = false;
	absTol = 1e-12;
	relTol = 1e-14;
	dtGuess = absTol;
}//===========================================

adtk_simulation_engine::~adtk_simulation_engine(){
    try{
        delete(traj);
    }catch(...){}
}//===========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return whether or not the simulation will run time in reverse
 */
bool adtk_simulation_engine::usesRevTime(){return revTime;}

/**
 *	@return whether or not the engine will be verbose in its outputs
 */
bool adtk_simulation_engine::isVerbose(){return verbose;}

/**
 *	@return the absolute tolerance for the engine, non-dimensional units
 */
double adtk_simulation_engine::getAbsTol(){return absTol;}

/**
 *	@return the relative tolerance for the engine, non-dimensional units
 */
double adtk_simulation_engine::getRelTol(){return relTol;}

/**
 *  Retrieve the trajectory. To avoid static casts in driver programs,
 *  we create several different getTraj() type functions that will perform
 *  the static cast and return the specific type of trajectory object rather
 *  than a generic one.
 *
 *  @return a CR3BP Trajectory object
 */
adtk_cr3bp_traj adtk_simulation_engine::getCR3BPTraj(){
    if(sysData->getType() == adtk_sys_data::CR3BP_SYS)
        /* Use a static cast to convert the general trajectory pointer
         * into a specific CR3BP trajectory pointer, then dereference
         * and return a COPY of the trajectory
         */
        return *( static_cast<adtk_cr3bp_traj*>(traj) );
    else{
        cout << "This trajectory was not propagated in CR3BP" << endl;
        throw;
    }
}//==============================================

/**
 *  Retrieve the trajectory. To avoid static casts in driver programs,
 *  we create several different getTraj() type functions that will perform
 *  the static cast and return the specific type of trajectory object rather
 *  than a generic one.
 *
 *  @return a BCR4BP, Rotating Coordinate Trajectory object
 */
adtk_bcr4bpr_traj adtk_simulation_engine::getBCR4BPRTraj(){
    if(sysData->getType() == adtk_sys_data::BCR4BPR_SYS)
        return *( static_cast<adtk_bcr4bpr_traj *>(traj) );
    else{
        cout << "This trajectory was not propagated in BCR4BP, Rot. Coord." << endl;
        throw;
    }
}//=====================================

/**
 *	Specify the system the engine will be using for integration
 *	@param d a pointer to the system data object (use &sys)
 */
void adtk_simulation_engine::setSysData(adtk_sys_data *d){ sysData = d; }

/**
 *	Specify whether or not the engine should run in reverse time
 *	@param b whether or not the engine should run in reverse time
 */
void adtk_simulation_engine::setRevTime(bool b){ revTime = b; }

/**
 *	Specify the verbosity of the engine
 *	@param b whether or not the engine should output verbose statements
 */
void adtk_simulation_engine::setVerbose(bool b){ verbose = b; }

/**
 *	Specify the absolute integration tolerance, non-dimensional units.
 *	The default value is 1e-12
 *	@param t the tolerance
 */
void adtk_simulation_engine::setAbsTol(double t){ absTol = t; }

/**
 *	Specify the absolute integration tolerance, non-dimensional units
 *	The default value is 1e-14
 *	@param t the tolerance
 */
void adtk_simulation_engine::setRelTol(double t){ relTol = t; }

//-----------------------------------------------------
// 		Simulation Functions
//-----------------------------------------------------

/**
 *	Run a simulation given a set of initial conditions and run time. It is
 *	assumed that t0 = 0
 *	@param ic a 6-element array containting the non-dimensional initial state
 *	@param tof the total integration time, or time-of-flight
 */
void adtk_simulation_engine::runSim(double *ic, double tof){
	this->runSim(ic, 0, tof);
}//=======================================================

/**
 *	Run a simulation in the specified system starting with a set of initial conditions,
 *  at a specified initial time, and integrating for a specified time-of-flight
 *	@param ic a 6-element array of non-dimensional initial states
 *	@param t0 the time at the start of the integration, non-dimensional units
 *	@param tof time-of-flight, non-dimensional time units
 */
void adtk_simulation_engine::runSim(double *ic, double t0, double tof){
	// Create time span using revTime to adjust limits
	double t_span[2] = {t0, 0};
	t_span[1] = revTime ? t0 - tof : t0 + tof;

    cout << "Running sim..." << endl;
	switch(sysData->getType()){
		case adtk_sys_data::CR3BP_SYS:
		{
            // Initialize trajectory (will only have one set of values)
            traj = new adtk_cr3bp_traj();
			break;
        }
		case adtk_sys_data::BCR4BPR_SYS:
            traj = new adtk_bcr4bpr_traj();
			break;
		default:
			cout << "Cannnot simulate with system type: " << sysData->getTypeStr() << endl;
			return;
	}

    // Run the simulation
    integrate(ic, t_span, 2);
}//============================================

//-----------------------------------------------------
// 		Numerical Integration
//-----------------------------------------------------

/**
 *  Integrate the 6 state EOMs and 36 STM EOMs with additional integration as required by 
 *  specific systems.
 *
 *  This function uses GSL's explicity embedded Runge-Kutta Dormand-Prince (8,9) method
 *
 *  @param ic a 6-element initial state for the trajectory
 *  @param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  @param t_dim the dimension of t
 */
void adtk_simulation_engine::integrate(double ic[], double t[], int t_dim){
    // Default IC dimension
    const int maxICDim = 48;
    int ic_dim = 42;

    if(simpleIntegration){
        ic_dim = 6;
    }else{
        if(sysData->getType() == adtk_sys_data::BCR4BPR_SYS){
            ic_dim = 48;    // requires 6 extra states for numerically integrated epoch dependencies
        }
    }

    // Construct the full 42-element IC from the state ICs plus the STM ICs
    double fullIC[maxICDim] = {0};  // store extra spots (only takes up a few bytes anyway)
    copy(ic, ic+6, fullIC);         // Copy in the 6 position & velocity states

    if(ic_dim > 6){
        fullIC[6] = 1;          // Make the identity matrix
        fullIC[13] = 1;
        fullIC[20] = 1;
        fullIC[27] = 1;
        fullIC[34] = 1;
        fullIC[41] = 1;
    }

    int steps = 0;                                  // count number of integration steps
    double *y = new double[maxICDim];                 // hold ONE integrated state; resets every step 

    // Define the step type (or, which integrator are we using?)
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_msadams;

    int (*eomFcn)(double, const double[], double[], void*) = 0;     // Pointer for the EOM function
    void *params = 0;   // Pointer to additional parameters
    // Choose EOM function based on system type and simplicity
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
        {
            eomFcn = simpleIntegration ? &cr3bp_simple_EOMs : &cr3bp_EOMs;
            
            adtk_cr3bp_sys_data *cr3bpData = static_cast<adtk_cr3bp_sys_data *>(sysData);
            double mu = cr3bpData->getMu();
            params = &mu;
            break;
        }
        case adtk_sys_data::BCR4BPR_SYS:
        {
            eomFcn = &bcr4bpr_EOMs;

            // Cast sys data to the specific system type
            adtk_bcr4bpr_sys_data thisSysData = *(static_cast<adtk_bcr4bpr_sys_data *>(sysData));

            double theta0 = 0;  // TODO make initializing these more intelligent
            double phi0 = 0;
            double gamma = 0;
            adtk_bcr4bpr_eomData eomData(thisSysData, theta0, phi0, gamma);
            
            params = &eomData;
            break;
        }
        default:
            throw;
    }

    // Create a system to integrate; we don't include a Jacobian (NULL)
    gsl_odeiv2_system sys = {eomFcn, NULL, static_cast<size_t>(ic_dim), params};
    
    // Create the driver
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, stepType, dtGuess, absTol, relTol);

    // Pointer to time-stepping function
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(stepType, ic_dim);
    gsl_odeiv2_step_set_driver(s, d);

    // Define a control that will keep the error in the state y within the specified tolerances
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new (absTol, relTol);
    gsl_odeiv2_control_set_driver(c, d);

    // Allocate space for the integrated solution to evolve in
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(ic_dim);
    gsl_odeiv2_evolve_set_driver(e, d);

    // Set y equal to the initial state
    copy(fullIC, fullIC+ic_dim, y);

    // Save the initial state, time, and STM
    saveIntegratedData(y, t[0], true);

    steps++;    // We've put in the IC, next step will be a new state

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
            saveIntegratedData(y, t0, false);
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
            saveIntegratedData(y, t0, false);
            steps++;
        }
    }

    // Clean Up
    gsl_odeiv2_driver_free(d);
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    delete(y);
    
    // Check lengths of vectors and set the numPoints value in traj
    traj->setLength();
}//================================END of cr3bp_integrate

/**
 *  Take an integrated state and the time, and save those variables into the 
 *  appropriate vectors in the trajectory data object
 *
 *  @param y a pointer to the 42-element state array computed by the integrator
 *  @param t the time at which the integrated state occurs
 *  @param first a flag telling the function if this is the first data point or not.
 *  A value of true indicates it IS the first, false indicates it is NOT
 */
void adtk_simulation_engine::saveIntegratedData(double *y, double t, bool first){

    // Grab pointers to the trajectory object's vectors
    vector<double>* state = traj->getState();       // hold the entire integrated state
    vector<double>* times = traj->getTime();        // hold all times along trajectory
    vector<adtk_matrix>* allSTM = traj->getSTM();   // hold all STM along trajectory

    // Save the position and velocity states
    for(int i = 0; i < 6; i++){
        if(first)
            state->at(i) = y[i];
        else
            state->push_back(y[i]);
    }

    printf("t=%5.2f :: %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);

    // Save time
    if(first)
        times->at(0) = t;
    else
        times->push_back(t);

    // Save STM
    double stmElm[36];
    copy(y+6, y+42, stmElm);
    if(first)
        allSTM->at(0) = adtk_matrix(6,6,stmElm);
    else
        allSTM->push_back(adtk_matrix(6,6,stmElm));

    // Compute acceleration
    double s[6] = {0};
    copy(y, y+6, s);
    double dsdt[6] = {0};
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
        {
            // Use the simple EOMs to compute the velocity/acceleration. 
            adtk_cr3bp_sys_data *cr3bpData = static_cast<adtk_cr3bp_sys_data *>(sysData);
            double mu = cr3bpData->getMu();
            cr3bp_simple_EOMs(0, y, dsdt, &mu);

            // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
            adtk_cr3bp_traj *cr3bpTraj = static_cast<adtk_cr3bp_traj*>(traj);
            
            if(first)
                cr3bpTraj->getJC()->at(0) = cr3bp_getJacobi(y, mu);
            else
                cr3bpTraj->getJC()->push_back(cr3bp_getJacobi(y, mu));
            break;
        }
        case adtk_sys_data::BCR4BPR_SYS:
            // Compute acceleration for BCR4BP
            break;
        default:
            cout << "Unknown sim type " << sysData->getTypeStr() << endl;
    }

    // Save the accelerations
    if(first){
        printf("1st Accel: %12.4f %12.4f %12.4f\n", dsdt[3], dsdt[4], dsdt[5]);
        state->at(6) = dsdt[3];
        state->at(7) = dsdt[4];
        state->at(8) = dsdt[5];
    }else{
        state->push_back(dsdt[3]);
        state->push_back(dsdt[4]);
        state->push_back(dsdt[5]);
    }
}//=========================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------