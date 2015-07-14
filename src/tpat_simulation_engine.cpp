/**
 *  @file tpat_simulation_engine.cpp
 *
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
 
#include "tpat_simulation_engine.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_bcr4bpr_nodeset.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_bcr4bpr_traj.hpp"
#include "tpat_body_data.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_cr3bp_traj.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Construct a new simulation engine. 
 *
 *  Most variables will be intiailized, but
 *  you MUST set the system data via <tt>setSysData()</tt>
 */
tpat_simulation_engine::tpat_simulation_engine(){
    events.clear();
    eventOccurs.clear();
    printVerb(verbose, "Created Simulation Engine\n");
}//===========================================

/**
 *  @brief Construct a simulation engine for a specific dynamical system
 *
 *  This constructor will also create the default crash event detectors
 *
 *  @param data a pointer to a system data object
 */
tpat_simulation_engine::tpat_simulation_engine(tpat_sys_data data){
    events.clear();
    eventOccurs.clear();
    sysData = data;
    createCrashEvents();
    printVerb(verbose, "Created Simulation Engine for %s system\n", data.getTypeStr().c_str());
}//===========================================

/**
 *  @brief Copy constructor
 *  @param s a simulation engine 
 */
tpat_simulation_engine::tpat_simulation_engine(const tpat_simulation_engine& s){
    copyEngine(s);
}//=====================================

/**
 *  @brief Free memory and clean up
 */
tpat_simulation_engine::~tpat_simulation_engine(){
    printVerb(verbose, "Destroying simulation engine...\n");
    reset();    // Function handles deallocation and resetting of data
}//===========================================

/**
 *  Copy data from an input engine to this one
 *  @param s an input simulation engine
 */
void tpat_simulation_engine::copyEngine(const tpat_simulation_engine &s){
    sysData = s.sysData;

    // Copy the trajectory object using the correct casting
    switch(sysData.getType()){
        case tpat_sys_data::CR3BP_SYS:
            traj = new tpat_cr3bp_traj(*(static_cast<tpat_cr3bp_traj *>(s.traj)));
            break;
        case tpat_sys_data::BCR4BPR_SYS:
            traj = new tpat_bcr4bpr_traj(*(static_cast<tpat_bcr4bpr_traj *>(s.traj)));
            break;
        default:
            traj = new tpat_trajectory(*(s.traj));
    }

    eomParams = 0;                              // void*, will get set again by the runSim() method
    revTime = s.revTime;
    verbose = s.verbose;
    varStepSize = s.varStepSize;
    simpleIntegration = s.simpleIntegration;
    isClean = s.isClean;
    absTol = s.absTol;
    relTol = s.relTol;
    dtGuess = s.dtGuess;
    numSteps = s.numSteps;
    events = s.events;
    eventOccurs = s.eventOccurs;
}//=====================================

/**
 *  Create default crash events for the system
 */
void tpat_simulation_engine::createCrashEvents(){
    if(!madeCrashEvents){
        for(int p = 0; p < sysData.getNumPrimaries(); p++){
            // Put primary index # into an array, create event
            double Pix = (double)p;
            tpat_event crashEvt(&sysData, tpat_event::CRASH, -1, true, &Pix);
            // Add event to list by default
            addEvent(crashEvt);
        }
        madeCrashEvents = true;
    }else{
        printWarn("Crash events have already been created!\n");
    }
}//=====================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *  @brief Assignment operator; make this engine equal another by copying its data
 *  @param s another simulation engine
 */
tpat_simulation_engine& tpat_simulation_engine::operator =(const tpat_simulation_engine& s){
    copyEngine(s);
    return *this;
}//=====================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return whether or not the simulation will run time in reverse
 */
bool tpat_simulation_engine::usesRevTime() const {return revTime;}

/**
 *	@return whether or not the engine will be verbose in its outputs
 */
bool tpat_simulation_engine::isVerbose() const {return verbose;}

/**
 *  @return whether or not the engine uses variable step size
 */
bool tpat_simulation_engine::usesVarStepSize() const { return varStepSize; }

/**
 *	@return the absolute tolerance for the engine, non-dimensional units
 */
double tpat_simulation_engine::getAbsTol() const {return absTol;}

/**
 *	@return the relative tolerance for the engine, non-dimensional units
 */
double tpat_simulation_engine::getRelTol() const {return relTol;}

/**
 *  @return the number of steps the integrator will be forced to take.
 *  The integrator may take intermediate steps between those enforced
 *  by the algorithm, but only <tt>numSteps</tt> data points will be output.
 */
int tpat_simulation_engine::getNumSteps() const { return numSteps; }

/**
 *  @brief Retrieve the trajectory as a generic trajectory object
 *  @return a trajectory object
 */
tpat_trajectory tpat_simulation_engine::getTraj() const {
    // Make a copy and return that
    tpat_trajectory temp (*traj);
    return temp;
}

/**
 *  @brief Retrieve the CR3BP trajectory. 
 *
 *  To avoid static casts in driver programs,
 *  we create several different getTraj() type functions that will perform
 *  the static cast and return the specific type of trajectory object rather
 *  than a generic one.
 *
 *  @return a CR3BP Trajectory object
 */
tpat_cr3bp_traj tpat_simulation_engine::getCR3BPTraj() const{
    if(sysData.getType() == tpat_sys_data::CR3BP_SYS){
        /* Use a static cast to convert the general trajectory pointer
         * into a specific CR3BP trajectory pointer, then dereference
         * and return a COPY of the trajectory
         */
        tpat_cr3bp_traj temp( *(static_cast<tpat_cr3bp_traj *>(traj) ) );
        return temp;
    }else{
        throw tpat_exception("Wrong system type");
    }
}//==============================================

/**
 *  @brief Retrieve the BCR4BPR trajectory. 
 *
 *  To avoid static casts in driver programs,
 *  we create several different getTraj() type functions that will perform
 *  the static cast and return the specific type of trajectory object rather
 *  than a generic one.
 *
 *  @return a BCR4BP, Rotating Coordinate Trajectory object
 */
tpat_bcr4bpr_traj tpat_simulation_engine::getBCR4BPRTraj() const{
    if(sysData.getType() == tpat_sys_data::BCR4BPR_SYS){
        // Make a copy and return it
        tpat_bcr4bpr_traj temp( *( static_cast<tpat_bcr4bpr_traj *>(traj) ) );
        return temp;
    }
    else{
        throw tpat_exception("Wrong system type");
    }
}//=====================================

/**
 *  @brief Add an event for this integration
 *  @param type the event type
 *  @param dir -1 for negative direction, +1 for positive, 0 for both
 *  @param stop whether or not to stop integration when this event occurs
 */
void tpat_simulation_engine::addEvent(tpat_event::event_t type, int dir, bool stop){
    tpat_event temp(&sysData, type, dir, stop);
    events.push_back(temp);
}//======================================

/**
 *  @brief Add an event for this integration
 *  @param evt an event
 */
void tpat_simulation_engine::addEvent(tpat_event evt){
    events.push_back(evt);
}//======================================

/**
 *	@brief Specify the system the engine will be using for integration
 *	@param d a system data object
 */
void tpat_simulation_engine::setSysData(tpat_sys_data d){ sysData = d; }

/**
 *	@brief Specify whether or not the engine should run in reverse time
 *	@param b whether or not the engine should run in reverse time
 */
void tpat_simulation_engine::setRevTime(bool b){ revTime = b; }

/**
 *	@brief Specify the verbosity of the engine
 *	@param b whether or not the engine should output verbose statements
 */
void tpat_simulation_engine::setVerbose(bool b){ verbose = b; }

/**
 *  @brief Specify whether or not the engine should use variable step size.
 *  @param b whether or not the engine should use variable step size
 */
void tpat_simulation_engine::setVarStepSize(bool b){ varStepSize = b; }

/**
 *	@brief Specify the absolute integration tolerance, non-dimensional units.
 *	The default value is 1e-12
 *	@param t the tolerance
 */
void tpat_simulation_engine::setAbsTol(double t){ absTol = t; }

/**
 *	@brief Specify the absolute integration tolerance, non-dimensional units
 *	The default value is 1e-14
 *	@param t the tolerance
 */
void tpat_simulation_engine::setRelTol(double t){ relTol = t; }

/**
 *  @brief Specify the number of steps the integrator must take during the 
 *  the integration. 
 *
 *  Only these points will be output to the 
 *  trajectory object, although the GSL driver may take steps in between
 *  those specified to maintain numerical accuracy.
 */
void tpat_simulation_engine::setNumSteps(int n){ numSteps = n; }

//-----------------------------------------------------
// 		Simulation Functions
//-----------------------------------------------------

/**
 *	@brief Run a simulation given a set of initial conditions and run time. 
 *
 *  It is assumed that t0 = 0
 *	@param ic a 6-element array containting the non-dimensional initial state
 *	@param tof the total integration time, or time-of-flight (non-dim units)
 */
void tpat_simulation_engine::runSim(double *ic, double tof){
	runSim(ic, 0, tof);
}//=======================================================

/**
 *  @brief Run a simulation given a set of initial conditions and run time
 *
 *  It is assumed that t0 = 0
 *  @param ic a vector containing the IC (six non-dim states)
 *  @param tof the total integration time, or time-of-flight (non-dim units)
 */
void tpat_simulation_engine::runSim(std::vector<double> ic, double tof){
    runSim(ic, 0, tof);
}//=======================================================

/**
 *  @brief Run a simulation given a set of initial conditions and run time
 *
 *  It is assumed that t0 = 0
 *  @param ic a vector containing the IC (six non-dim states)
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units)
 */
void tpat_simulation_engine::runSim(std::vector<double> ic, double t0, double tof){
    if(ic.size() >= 6){
        std::vector<double> tempIC = ic;
        runSim(&(tempIC[0]), t0, tof);
    }else{
        throw tpat_exception("IC must have at least six elements");
    }
}//=======================================================

/**
 *	@brief Run a simulation in the specified system starting with a set of initial conditions,
 *  at a specified initial time, and integrating for a specified time-of-flight
 *	@param ic a 6-element array of non-dimensional initial states
 *	@param t0 the time at the start of the integration, non-dimensional units
 *	@param tof time-of-flight, non-dimensional time units
 */
void tpat_simulation_engine::runSim(double *ic, double t0, double tof){
    printVerbColor(verbose, GREEN, "Running simulation...\n");
    if(!isClean){
        cleanEngine();
    }

    std::vector<double> t_span;
    // Compute the final time based on whether or not we're using reverse time integration
    double tf = revTime ? t0 - std::abs(tof) : t0 + std::abs(tof);
    printVerb(verbose, "  time will span from %.3e to %.3e\n", t0, tf);
    printVerb(verbose, "  (Reverse Time is %s)\n", revTime ? "ON" : "OFF");

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

	switch(sysData.getType()){
		case tpat_sys_data::CR3BP_SYS:
		{
            // Initialize trajectory (will only have one set of values)
            printVerb(verbose, "  initializing CR3BP trajectory\n");
            tpat_cr3bp_sys_data data = *static_cast<tpat_cr3bp_sys_data *>(&sysData);
            traj = new tpat_cr3bp_traj(data);
			break;
        }
		case tpat_sys_data::BCR4BPR_SYS:
        {
            printVerb(verbose, "  initializing BCR4BPR trajectory\n");
            tpat_bcr4bpr_sys_data data = *static_cast<tpat_bcr4bpr_sys_data *>(&sysData);
            traj = new tpat_bcr4bpr_traj(data);
			break;
        }
		default:
			throw tpat_exception("Cannot simulate for this system type");
	}

    // Run the simulation
    integrate(ic, &(t_span.front()), varStepSize ? 2 : numSteps);

    isClean = false;
}//============================================

//-----------------------------------------------------
// 		Numerical Integration
//-----------------------------------------------------

/**
 *  @brief Integrate the 6 state EOMs and 36 STM EOMs with additional integration as required by 
 *  specific systems.
 *
 *  This function uses values stored in object-wide variables to determine the direction time flows,
 *  whether or not to use simple integration, and whether or not to use variable step sizes. Dad
 *  is saved to object-wide storage vectors via the <tt>saveIntegratedData()</tt> function.
 *
 *  @param ic a 6-element initial state for the trajectory
 *  @param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  @param t_dim the dimension of t
 */
void tpat_simulation_engine::integrate(double ic[], double t[], int t_dim){
    // Default IC dimension
    int ic_dim = 42;

    if(simpleIntegration){
        ic_dim = 6;
    }else{
        if(sysData.getType() == tpat_sys_data::BCR4BPR_SYS){
            ic_dim = 48;    // requires 6 extra states for numerically integrated epoch dependencies
        }
    }
    printVerb(verbose, "  IC has %d initial states\n", ic_dim);

    // Construct the full IC from the state ICs plus the STM ICs and any other ICs for more complex systems
    std::vector<double> fullIC(ic_dim, 0);
    std::copy(ic, ic+6, &(fullIC.front()));

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
    printVerb(verbose, "  using %s integration\n", simpleIntegration ? "simple (no STM)" : "full (+ STM)");
    int (*eomFcn)(double, const double[], double[], void*) = 0;     // Pointer for the EOM function
    switch(sysData.getType()){
        case tpat_sys_data::CR3BP_SYS:
            eomFcn = simpleIntegration ? &cr3bp_simple_EOMs : &cr3bp_EOMs;
            break;
        case tpat_sys_data::BCR4BPR_SYS:
            eomFcn = simpleIntegration ? &bcr4bpr_simple_EOMs : &bcr4bpr_EOMs;
            break;
        default:
            throw tpat_exception("Unknown sim system type");
    }

    // Create a system to integrate; we don't include a Jacobian (NULL)
    gsl_odeiv2_system sys = {eomFcn, NULL, static_cast<size_t>(ic_dim), eomParams};
    
    // Define ODE objects, define them conditionaly based on varStepSize
    gsl_odeiv2_step *s;
    gsl_odeiv2_control *c;
    gsl_odeiv2_evolve *e;
    gsl_odeiv2_driver *d;

    if(varStepSize){
        printVerb(verbose, "  variable step size, using Runge-Kutta Cash-Karp 4-5 method\n");
        // Allocate space for the stepping object; use the rkck algorithm (doesn't require driver)
        s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, ic_dim);
        // Define a control that will keep the error in the state y within the specified tolerances
        c = gsl_odeiv2_control_y_new (absTol, relTol);
        // Allocate space for the integrated solution to evolve in
        e = gsl_odeiv2_evolve_alloc(ic_dim);
    }else{
        printVerb(verbose, "  fixed step size, using Adams-Bashforth, Adams-Moulton method\n");
        // Allocate space for a driver; the msadams algorithm requires access to the driver
        double signed_dt = revTime ? -1*dtGuess : dtGuess;
        d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msadams, signed_dt, absTol, relTol);
        // Allocate space for the stepping object
        s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_msadams, ic_dim);
        gsl_odeiv2_step_set_driver(s, d);
    }

    // Save the initial state, time, and STM
    saveIntegratedData(y, t[0]);

    // Update all event functions with IC
    printVerb(verbose, "  sim will use %d event functions:\n", ((int)events.size()));
    for(int ev = 0; ev < ((int)events.size()); ev++){
        printVerb(verbose, "  >>%s\n", events.at(ev).getTypeStr());
        events.at(ev).updateDist(y, t[0]);
    }

    steps++;    // We've put in the IC, next step will be a new state
    int status; // integrator status
    bool killSim = false;   // flag to stop simulation

    if(t_dim == 2){
        double t0 = t[0], tf = t[1];            // start and finish times for integration
        double dt = tf > t0 ? dtGuess : -dtGuess;   // step size (initial guess)
        int sgn = tf > t0 ? 1 : -1;

        while (sgn*t0 < sgn*tf && !killSim){
            // apply integrator for one step; new state and time are stored in y and t0
            if(varStepSize){
                status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
            }else{
                status = gsl_odeiv2_driver_apply(d, &t0, tf, y);
            }

            if(status != GSL_SUCCESS){
                printErr("tpat_simulation_engine::integrate: GSL ERR: %s\n", gsl_strerror(status));
                throw tpat_diverge("Integration did not succeed");
            }

            killSim = locateEvents(y, t0);

            // Put newly integrated state and time into state vector
            if(!killSim)
                saveIntegratedData(y, t0);
            
            steps++;
        }
    }else{
        // Integrate each segment between the input times
        for (int j = 0; j < t_dim - 1; j++){
            // define start and end times; t0 will be updated by integrator
            double t0 = t[j], tf = t[j+1];
            double dt = tf > t0 ? dtGuess : -dtGuess;
            int sgn = tf > t0 ? 1 : -1;

            while(sgn*t0 < sgn*tf && !killSim){

                if(varStepSize){
                    status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
                }else{
                    status = gsl_odeiv2_driver_apply(d, &t0, tf, y);
                }

                if(status != GSL_SUCCESS){
                    printErr("tpat_simulation_engine::integrate: GSL ERR: %s\n", gsl_strerror(status));
                    throw tpat_diverge("Integration did not succeed");
                }

                killSim = locateEvents(y, t0);
            }

            if(killSim)
                break;

            // Add the newly integrated state and current time fo the state vector
            saveIntegratedData(y, t0);
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
    printVerbColor(verbose, GREEN, "  **Integration complete**\n  Total: %d data points\n", steps);

    // Summarize event occurrences
    for(int ev = 0; ev < (int)(eventOccurs.size())-1; ev+=2){
        printVerb(verbose, "  Event %d (%s) occured at step %d\n", eventOccurs[ev+1],
            events[eventOccurs[ev+1]].getTypeStr(), eventOccurs[ev]);
    }
    traj->setLength();
}//================================END of cr3bp_integrate

/**
 *  @brief Take an integrated state and the time, and save those variables into the 
 *  appropriate vectors in the trajectory data object
 *
 *  @param y a pointer to the 42-element state array computed by the integrator
 *  @param t the time at which the integrated state occurs
 */
void tpat_simulation_engine::saveIntegratedData(double *y, double t){

    // Grab pointers to the trajectory object's vectors
    std::vector<double>* state = traj->getState();       // hold the entire integrated state
    std::vector<double>* times = traj->getTime();        // hold all times along trajectory
    std::vector<tpat_matrix>* allSTM = traj->getSTM();   // hold all STM along trajectory

    // Save the position and velocity states
    for(int i = 0; i < 6; i++){
        state->push_back(y[i]);
    }

    // printf("t=%5.2f :: %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n", t, y[0], y[1], y[2], y[3], y[4], y[5]);

    // Save time
    times->push_back(t);

    // Save STM
    double stmElm[36];
    std::copy(y+6, y+42, stmElm);
    allSTM->push_back(tpat_matrix(6,6,stmElm));

    // Compute acceleration
    double dsdt[6] = {0};
    switch(sysData.getType()){
        case tpat_sys_data::CR3BP_SYS:
        {
            // Use the simple EOMs to compute the velocity/acceleration. 
            cr3bp_simple_EOMs(0, y, dsdt, eomParams);

            // Cast trajectory to a cr3bp_traj and then store a value for Jacobi Constant
            tpat_cr3bp_traj *cr3bpTraj = static_cast<tpat_cr3bp_traj*>(traj);
            tpat_cr3bp_sys_data *cr3bpData = static_cast<tpat_cr3bp_sys_data *>(&sysData);
            
            cr3bpTraj->getJC()->push_back(cr3bp_getJacobi(y, cr3bpData->getMu()));
            break;
        }
        case tpat_sys_data::BCR4BPR_SYS:
        {
            // Use simple EOMs to compute acceleration, stored in dsdt
            bcr4bpr_simple_EOMs(t, y, dsdt, eomParams);

            // Cast the trajectory to a BCR4BPR guy and save the dqdT data
            tpat_bcr4bpr_traj *bcrTraj = static_cast<tpat_bcr4bpr_traj*>(traj);
            bcrTraj->get_dqdT()->push_back(y[42]);
            bcrTraj->get_dqdT()->push_back(y[43]);
            bcrTraj->get_dqdT()->push_back(y[44]);
            bcrTraj->get_dqdT()->push_back(y[45]);
            bcrTraj->get_dqdT()->push_back(y[46]);
            bcrTraj->get_dqdT()->push_back(y[47]);
            break;
        }
        default:
            throw tpat_exception("Unsupported system type");
    }

    // Save the accelerations
    state->push_back(dsdt[3]);
    state->push_back(dsdt[4]);
    state->push_back(dsdt[5]);
}//=========================================

/**
 *  @brief Set the pointer for EOM Parameters for each type of system
 */
void tpat_simulation_engine::setEOMParams(){
    switch(sysData.getType()){
        case tpat_sys_data::CR3BP_SYS:
        {
            eomParams = static_cast<tpat_cr3bp_sys_data *>(&sysData);
            break;
        }
        case tpat_sys_data::BCR4BPR_SYS:
        {
            eomParams = static_cast<tpat_bcr4bpr_sys_data *>(&sysData);
            break;
        }
        default:
            throw tpat_exception("Unsupported system type");
    }
}//==============================================

/**
 *  @brief Locate event occurences as exactly as possible and determine if the simulation 
 *  should end because of the event.
 *
 *  Look through the list of events and check each one to see if it has occured. This 
 *  initial check is not highly accurate, but if an event is determined to have occured,
 *  a Newton-Raphson process is begun to locate the exact state and time of the event.
 *
 *  The initial check compares the current integrated state (passed in as <tt>y</tt>) with
 *  the previous integrated state (stored in the event object). If the sign (+/-) of the 
 *  distance to the event function changes, the the trajectory has triggered the event. This
 *  check is performed in the event objects <tt>crossedEvent()</tt> function.
 *
 *  Once an event has triggered the initial check, we create a 2-node nodeset beginning
 *  with the previous state and integrating for a time interval that should end near the 
 *  exact event state/time. The nodeset is constrained such that the first node cannot change
 *  and the final node is constrained to enforce the event condition. A correction engine
 *  is created and applied to the nodest to determine the exact location of the event.
 *
 *  Once the event is located (exactly), the full state (42 or 48 elements) is retrieved from
 *  the correction engine and saved to the trajectory data vectors like a regularly integrated
 *  state would be.
 *
 *  @param y the most recent state on the integrated arc.
 *  @param t the time associated with y
 *  @return whether or not the simulation should end (an event triggers killSim)
 *
 *  TODO: Implement a way to save information about which event is fired at which state
 */
bool tpat_simulation_engine::locateEvents(double *y, double t){
    int numPts = traj->getTime()->size();

    // Look through all events
    for(int ev = 0; ev < ((int)events.size()); ev++){
        // Don't trigger if only two points have been integrated
        if(events.at(ev).crossedEvent(y, t) && numPts > 1){

            printVerb(verbose, "  Event %d detected; searching for exact crossing\n", ev);
            if(verbose){ events.at(ev).printStatus(); }

            // Create a nodeset from the previous state (stored in the event) and
            // integrating forwards for half the time between this state and the last one
            double t0 = traj->getTime(-2);          // Time from the state before last
            double ti = traj->getTime(-1);          // Time from the previous state
            double tof = t - t0 - 0.5*(t - ti);     // Approx. TOF 

            // Copy 6-element IC into vector - Use the state from two iterations ago to avoid
            // numerical problems when the previous state is REALLY close to the event
            std::vector<double> generalIC = traj->getState(-2);

            if(verbose){
                printColor(BLUE, "Num Pts = %d\n", numPts);
                printColor(BLUE, "t(now) = %f\nt(prev) = %f\nt(prev-1) = %f\n", t, 
                    traj->getTime(-1), traj->getTime(-2));
                printColor(BLUE, "tof = %f\n", tof);
                printColor(BLUE, "ic = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", generalIC[0],
                    generalIC[1], generalIC[2], generalIC[3], generalIC[4], generalIC[5]);
                // waitForUser();
            }   

            tpat_correction_engine corrector;
            corrector.setVarTime(true);
            corrector.setTol(absTol);
            corrector.setVerbose(verbose);
            corrector.setFindEvent(true);   // apply special settings to minimize computations

            switch(sysData.getType()){
                case tpat_sys_data::CR3BP_SYS:
                {
                    // Get the address of the IC
                    double *ic = &(generalIC[0]);

                    // Copy system data object
                    tpat_cr3bp_sys_data crSysData(*static_cast<tpat_cr3bp_sys_data *>(&sysData));

                    // Create a nodeset for this particular type of system
                    printVerb(verbose, "  Creating nodeset for event location\n");
                    tpat_cr3bp_nodeset eventNodeset(ic, crSysData, tof, 2, tpat_nodeset::TIME);

                    // Constraint to keep first node unchanged
                    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, ic, 6);

                    // Constraint to enforce event
                    tpat_constraint eventCon(events.at(ev).getConType(),
                        events.at(ev).getConNode(), events.at(ev).getConData());

                    eventNodeset.addConstraint(fixFirstCon);
                    eventNodeset.addConstraint(eventCon);

                    if(verbose){ eventNodeset.print(); }

                    printVerb(verbose, "  Applying corrections process to locate event\n");
                    try{
                        corrector.correct_cr3bp(&eventNodeset);
                    }catch(tpat_diverge &e){
                        printErr("Unable to locate event; corrector diverged\n");
                        return false;
                    }catch(tpat_linalg_err &e){
                        printErr("LinAlg Err while locating event; bug in corrector!\n");
                        return false;
                    }

                    // Because we set findEvent to true, this output nodeset should contain
                    // the full (42 or 48 element) final state
                    tpat_cr3bp_nodeset correctedNodes = corrector.getCR3BPOutput();
                    std::vector<double> *nodes = correctedNodes.getNodes();

                    // event time is the TOF of corrected path + time at the state we integrated from
                    double eventTime = correctedNodes.getTOF(0) + t0;

                    // Use the data stored in nodes and save the state and time of the event occurence
                    saveIntegratedData(&(nodes->at(6)), eventTime);
                    break;
                }
                case tpat_sys_data::BCR4BPR_SYS:
                {
                    // **** Make sure you fix the epoch of the first node as well as the states
                    generalIC.push_back(t0);
                    double *ic = &(generalIC[0]);

                    // Copy system data object
                    tpat_bcr4bpr_sys_data bcSysData(*static_cast<tpat_bcr4bpr_sys_data *>(&sysData));

                    // Create a nodeset for this particular type of system
                    printVerb(verbose, "  Creating nodeset for event location\n");
                    tpat_bcr4bpr_nodeset eventNodeset(ic, bcSysData, t0,
                        tof, 2, tpat_nodeset::TIME);

                    // Constraint to keep first node unchanged
                    tpat_constraint fixFirstCon(tpat_constraint::STATE, 0, ic, 7);

                    // Constraint to enforce event
                    tpat_constraint eventCon(events.at(ev).getConType(),
                        events.at(ev).getConNode(), events.at(ev).getConData());

                    eventNodeset.addConstraint(fixFirstCon);
                    eventNodeset.addConstraint(eventCon);

                    if(verbose){ eventNodeset.print(); }

                    printVerb(verbose, "  Applying corrections process to locate event\n");
                    try{
                        corrector.correct_bcr4bpr(&eventNodeset);
                    }catch(tpat_diverge &e){
                        printErr("Unable to locate event; corrector diverged\n");
                        return false;
                    }catch(tpat_linalg_err &e){
                        printErr("LinAlg Err while locating event; bug in corrector!\n");
                        return false;
                    }

                    // Because we set findEvent to true, this output nodeset should contain
                    // the full (42 or 48 element) final state
                    tpat_bcr4bpr_nodeset correctedNodes = corrector.getBCR4BPROutput();
                    std::vector<double> *nodes = correctedNodes.getNodes();

                    // event time is the TOF of corrected path + time at the state we integrated from
                    double eventTime = correctedNodes.getTOF(0) + t0;

                    // Use the data stored in nodes and save the state and time of the event occurence
                    saveIntegratedData(&(nodes->at(6)), eventTime);
                    break;
                }
                default:
                    throw tpat_exception("Unsupported system type");
            }
            
            // Remember that this event has occured; step # is one less than the current size
            // of the trajectory's time vector
            int timeSize = traj->getTime()->size();
            eventOccurs.push_back(timeSize - 1);
            eventOccurs.push_back(ev);

            if(events.at(ev).stopOnEvent()){
                printVerbColor(verbose, GREEN, "**Completed Event Location, ending integration**\n");
                return true;    // Tell the simulation to stop
            }else{
                printVerbColor(verbose, GREEN, "**Completed Event Location, continuing integration**\n");
                return false;
            }
        }// end of If(hasCrossed)
        // Save the distance and current state to the event
        events.at(ev).updateDist(y, t);
    }// end of loop

    return false;
}//========================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  @brief Clean out the trajectory storage variable so a new simulation can be run and store its data
 */
void tpat_simulation_engine::cleanEngine(){
    printVerb(verbose, "Cleaning the engine...\n");
    delete traj;   // de-allocate the memory
    traj = 0;       // set pointer to 0 (null pointer)
    eomParams = 0;
    eventOccurs.clear();

    isClean = true;
}//====================================================

/**
 *  @brief Reset all variables and options
 *
 *  Completely resets the simulation engine, reverting all variables (including ones the user
 *  has modified with set() functions) to their default values.
 */
void tpat_simulation_engine::reset(){
    if(!isClean)
        cleanEngine();

    events.clear();
    eventOccurs.clear();
    revTime = false;
    verbose = false;
    varStepSize = true;
    absTol = 1e-12;
    relTol = 1e-14;
    dtGuess = 1e-6;
    numSteps = 1000;
    madeCrashEvents = false;
}//==========================================

/**
 *  Clear all events from the simulation, including any created by default.
 */
void tpat_simulation_engine::clearEvents(){
    printVerb(verbose, "Clearing all events...\n");
    events.clear();
}//==========================================