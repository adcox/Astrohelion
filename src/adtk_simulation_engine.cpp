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
#include "adtk_constants.hpp"
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

/**
 *  Construct a new simulation engine. Most variables will be intiailized, but
 *  you MUST set the system data via <tt>setSysData()</tt>
 */
adtk_simulation_engine::adtk_simulation_engine(){
    printMessage("Created Simulation Engine\n");
}//===========================================

/**
 *  Construct a simulation engine for a specific dynamical system
 *  @param data a pointer to a system data object
 */
adtk_simulation_engine::adtk_simulation_engine(adtk_sys_data *data){
    sysData = data;
    printMessage("Created Simulation Engine for %s system\n", data->getTypeStr().c_str());
}//===========================================

adtk_simulation_engine::~adtk_simulation_engine(){
    printMessage("Destroying simulation engine...\n");
    reset();    // Function handles deallocation and resetting of data
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
 *  @return whether or not the engine uses variable step size
 */
bool adtk_simulation_engine::usesVarStepSize(){ return varStepSize; }

/**
 *	@return the absolute tolerance for the engine, non-dimensional units
 */
double adtk_simulation_engine::getAbsTol(){return absTol;}

/**
 *	@return the relative tolerance for the engine, non-dimensional units
 */
double adtk_simulation_engine::getRelTol(){return relTol;}

/**
 *  @return the number of steps the integrator will be forced to take.
 *  The integrator may take intermediate steps between those enforced
 *  by the algorithm, but only <tt>numSteps</tt> data points will be output.
 */
int adtk_simulation_engine::getNumSteps(){ return numSteps; }

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
        printf("Wrong system type: %s\n", sysData->getTypeStr().c_str());
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
        printf("Wrong system type: %s\n", sysData->getTypeStr().c_str());
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
 *  Specify whether or not the engine should use variable step size.
 *  @param b whether or not the engine should use variable step size
 */
void adtk_simulation_engine::setVarStepSize(bool b){ varStepSize = b; }

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

/**
 *  Specify the number of steps the integrator must take during the 
 *  the integration. Only these points will be output to the 
 *  trajectory object, although the GSL driver may take steps in between
 *  those specified to maintain numerical accuracy.
 */
void adtk_simulation_engine::setNumSteps(int n){ numSteps = n; }

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
    printMessage("Running simulation...\n");
    if(!isClean){
        cleanEngine();
    }

    vector<double> t_span;
    // Compute the final time based on whether or not we're using reverse time integration
    double tf = revTime ? t0 - tof : t0 + tof;
    printMessage("  time will span from %.4f to %.4f\n", t0, tf);

    if(varStepSize){
        t_span.reserve(2);
        t_span.push_back(t0);
        t_span.push_back(tf);
    }else{
        t_span.reserve(numSteps);
        double dt = (tf - t0)/(numSteps - 1);

        for(int n = 0; n < numSteps; n++){
            t_span.push_back(t0 + dt*n);
        }
    }

	switch(sysData->getType()){
		case adtk_sys_data::CR3BP_SYS:
		{
            // Initialize trajectory (will only have one set of values)
            printMessage("  initializing CR3BP trajectory\n");
            adtk_cr3bp_sys_data *data = static_cast<adtk_cr3bp_sys_data *>(sysData);
            traj = new adtk_cr3bp_traj(*data);
			break;
        }
		case adtk_sys_data::BCR4BPR_SYS:
        {
            printMessage("  initializing BCR4BPR trajectory\n");
            adtk_bcr4bpr_sys_data *data = static_cast<adtk_bcr4bpr_sys_data *>(sysData);
            traj = new adtk_bcr4bpr_traj(*data);
			break;
        }
		default:
			cout << "Cannnot simulate with system type: " << sysData->getTypeStr() << endl;
			return;
	}

    // Run the simulation
    integrate(ic, &(t_span.front()), varStepSize ? 2 : numSteps);

    isClean = false;
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
    int ic_dim = 42;

    if(simpleIntegration){
        ic_dim = 6;
    }else{
        if(sysData->getType() == adtk_sys_data::BCR4BPR_SYS){
            ic_dim = 48;    // requires 6 extra states for numerically integrated epoch dependencies
        }
    }
    printMessage("  IC has %d initial states\n", ic_dim);

    // Construct the full IC from the state ICs plus the STM ICs and any other ICs for more complex systems
    vector<double> fullIC(ic_dim, 0);
    copy(ic, ic+6, &(fullIC.front()));

    if(ic_dim > 6){
        fullIC.at(6) = 1;       // STM initial condition: 6x6 identity matrix
        fullIC.at(13) = 1;
        fullIC.at(20) = 1;
        fullIC.at(27) = 1;
        fullIC.at(34) = 1;
        fullIC.at(41) = 1;      // Elements 42-48 for BCR4BPR ICs are all zero
    }

    int steps = 0;                      // count number of integration steps
    double *y = &(fullIC.front());      // double "array" that is passed to the integrator

    // set the parameters that will be given to the EOM function
    setEOMParams();

    // Choose EOM function based on system type and simplicity
    printMessage("  using %s integration\n", simpleIntegration ? "simple (no STM)" : "full (+ STM)");
    int (*eomFcn)(double, const double[], double[], void*) = 0;     // Pointer for the EOM function
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
        {
            eomFcn = simpleIntegration ? &cr3bp_simple_EOMs : &cr3bp_EOMs;
            break;
        }
        case adtk_sys_data::BCR4BPR_SYS:
        {
            eomFcn = simpleIntegration ? &bcr4bpr_simple_EOMs : &bcr4bpr_EOMs;
            break;
        }
        default:
            cout << "Unknown sim type " << sysData->getTypeStr() << endl;
            throw;
    }

    // Create a system to integrate; we don't include a Jacobian (NULL)
    gsl_odeiv2_system sys = {eomFcn, NULL, static_cast<size_t>(ic_dim), eomParams};
    
    // Define ODE objects, define them conditionaly based on varStepSize
    gsl_odeiv2_step *s;
    gsl_odeiv2_control *c;
    gsl_odeiv2_evolve *e;
    gsl_odeiv2_driver *d;

    if(varStepSize){
        printMessage("  variable step size, using Runge-Kutta Cash-Karp 4-5 method\n");
        // Allocate space for the stepping object; use the rkck algorithm (doesn't require driver)
        s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, ic_dim);
        // Define a control that will keep the error in the state y within the specified tolerances
        c = gsl_odeiv2_control_y_new (absTol, relTol);
        // Allocate space for the integrated solution to evolve in
        e = gsl_odeiv2_evolve_alloc(ic_dim);
    }else{
        printMessage("  fixed step size, using Adams-Bashforth, Adams-Moulton method\n");
        // Allocate space for a driver; the msadams algorithm requires access to the driver
        d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msadams, dtGuess, absTol, relTol);
        // Allocate space for the stepping object
        s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_msadams, ic_dim);
        gsl_odeiv2_step_set_driver(s, d);
    }

    // Save the initial state, time, and STM
    saveIntegratedData(y, t[0], true);

    steps++;    // We've put in the IC, next step will be a new state
    int status; // integrator status
    if(t_dim == 2){
        double t0 = t[0], tf = t[1];            // start and finish times for integration
        double dt = tf > t0 ? dtGuess : -dtGuess;   // step size (initial guess)
        int sgn = tf > t0 ? 1 : -1;

        while (sgn*t0 < sgn*tf){
            // apply integrator for one step; new state and time are stored in y and t0
            if(varStepSize){
                status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
            }else{
                status = gsl_odeiv2_driver_apply(d, &t0, tf, y);
            }

            if(status != GSL_SUCCESS){
                printf("Integration did not succeed:\n  GSL error: %s\n", gsl_strerror(status));
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

                if(varStepSize){
                    status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
                }else{
                    status = gsl_odeiv2_driver_apply(d, &t0, tf, y);
                }

                if(status != GSL_SUCCESS){
                    printf("Integration did not succeed:\n  GSL error: %s\n", gsl_strerror(status));
                    break;
                }
            }

            // Add the newly integrated state and current time fo the state vector
            saveIntegratedData(y, t0, false);
            steps++;
        }
    }

    // Clean Up - some have been allocated, some have not
    if(varStepSize){
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c);
    }else{
        gsl_odeiv2_driver_free(d);
    }
    gsl_odeiv2_step_free(s);
    
    // Check lengths of vectors and set the numPoints value in traj
    printMessage("  **Integration complete**\n  Total: %d data points\n", steps);
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

    // printf("t=%5.2f :: %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);

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
    double dsdt[6] = {0};
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
        {
            // Use the simple EOMs to compute the velocity/acceleration. 
            cr3bp_simple_EOMs(0, y, dsdt, eomParams);

            // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
            adtk_cr3bp_traj *cr3bpTraj = static_cast<adtk_cr3bp_traj*>(traj);
            adtk_cr3bp_sys_data *cr3bpData = static_cast<adtk_cr3bp_sys_data *>(sysData);
            
            if(first)
                cr3bpTraj->getJC()->at(0) = cr3bp_getJacobi(y, cr3bpData->getMu());
            else
                cr3bpTraj->getJC()->push_back(cr3bp_getJacobi(y, cr3bpData->getMu()));
            break;
        }
        case adtk_sys_data::BCR4BPR_SYS:
        {
            bcr4bpr_simple_EOMs(t, y, dsdt, eomParams);

            // Cast the trajectory to a BCR4BPR guy and save the dqdT data
            adtk_bcr4bpr_traj *bcrTraj = static_cast<adtk_bcr4bpr_traj*>(traj);
            if(first){
                vector<double>* vecPtr = bcrTraj->get_dqdT();
                copy(&(vecPtr->front()), &(vecPtr->front())+6, y+42);
            }
            else{
                bcrTraj->get_dqdT()->push_back(y[42]);
                bcrTraj->get_dqdT()->push_back(y[43]);
                bcrTraj->get_dqdT()->push_back(y[44]);
                bcrTraj->get_dqdT()->push_back(y[45]);
                bcrTraj->get_dqdT()->push_back(y[46]);
                bcrTraj->get_dqdT()->push_back(y[47]);
            }
            break;
        }
        default:
            cout << "Unknown sim type " << sysData->getTypeStr() << endl;
    }

    // Save the accelerations
    if(first){
        // printf("1st Accel: %14.8f %14.8f %14.8f\n", dsdt[3], dsdt[4], dsdt[5]);
        copy(dsdt+3, dsdt+6, &(state->front())+6);
    }else{
        state->push_back(dsdt[3]);
        state->push_back(dsdt[4]);
        state->push_back(dsdt[5]);
    }
}//=========================================

/**
 *  Set the pointer for EOM Parameters for each type of system
 */
void adtk_simulation_engine::setEOMParams(){
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
        {
            adtk_cr3bp_sys_data *cr3bpData = static_cast<adtk_cr3bp_sys_data *>(sysData);
            double *mu = new double(cr3bpData->getMu());
            eomParams = mu;
            break;
        }
        case adtk_sys_data::BCR4BPR_SYS:
        {
            // Cast sys data to the specific system type
            adtk_bcr4bpr_sys_data thisSysData = *(static_cast<adtk_bcr4bpr_sys_data *>(sysData));

            double theta0 = 0;  // TODO make initializing these more intelligent
            double phi0 = 0;
            double gamma = 5.14*PI/180;
            // Create a NEW object so that it isn't deleted once we exit this scope
            adtk_bcr4bpr_eomData *eomData = new adtk_bcr4bpr_eomData(thisSysData, theta0, phi0, gamma);
            
            eomParams = eomData;
            break;
        }
        default:
            cout << "Unknown sim type " << sysData->getTypeStr() << endl;
    }
}//==============================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  Clean out the trajectory storage variable so a new simulation can be run and store its data
 */
void adtk_simulation_engine::cleanEngine(){
    printMessage("Cleaning the engine...\n");
    delete traj;   // de-allocate the memory
    traj = 0;       // set pointer to 0 (null pointer)

    // Free eomParams; MUST cast to correct type before freeing so that the appropriate
    // destructor is called
    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
            delete (double *)eomParams;
            break;
        case adtk_sys_data::BCR4BPR_SYS:
            delete (adtk_bcr4bpr_eomData *)eomParams;
            break;
        default:
            cout << "Unrecongized system type; cannot free eomParams!!" << endl;
    }
    eomParams = 0;

    isClean = true;
}//====================================================

/**
 *  A wrapper function to print a message
 */
void adtk_simulation_engine::printMessage(const char * format, ...){
    if(verbose){
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
    }
}

/**
 *  Completely resets the simulation engine, reverting all variables (including ones the user
 *  has modified with set() functions) to their default values.
 */
void adtk_simulation_engine::reset(){
    if(!isClean)
        cleanEngine();

    revTime = false;
    verbose = false;
    varStepSize = true;
    absTol = 1e-12;
    relTol = 1e-14;
    dtGuess = 1e-6;
    numSteps = 1000;
}