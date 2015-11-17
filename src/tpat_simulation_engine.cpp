/**
 *  @file tpat_simulation_engine.cpp
 *  @brief Performs numerical integration on a set of initial conditions in 
 *  any dynamic model and system
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
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_traj_bcr4bpr.hpp"
#include "tpat_body_data.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_model.hpp"
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
tpat_simulation_engine::tpat_simulation_engine(tpat_sys_data *data){
    events.clear();
    eventOccurs.clear();
    sysData = data;
    createCrashEvents();
    printVerb(verbose, "Created Simulation Engine for %s system\n", data->getTypeStr().c_str());
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

    // Copy the trajectory object using the correct casting
    switch(sysData->getType()){
        case tpat_sys_data::CR3BP_SYS:
            sysData = static_cast<tpat_sys_data_cr3bp *>(s.sysData);
            traj = new tpat_traj_cr3bp(*static_cast<tpat_traj_cr3bp *>(s.traj));
            break;
        case tpat_sys_data::CR3BP_LTVP_SYS:
            sysData = static_cast<tpat_sys_data_cr3bp_ltvp *>(s.sysData);
            traj = new tpat_traj_cr3bp_ltvp(*static_cast<tpat_traj_cr3bp_ltvp *>(s.traj));
            break;
        case tpat_sys_data::BCR4BPR_SYS:
            sysData = static_cast<tpat_sys_data_bcr4bpr *>(s.sysData);
            traj = new tpat_traj_bcr4bpr(*static_cast<tpat_traj_bcr4bpr *>(s.traj));
            break;
        default:
            throw tpat_exception("tpat_simulation_engine::copyEngine: Cannot copy engine with unknown type");
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
        for(int p = 0; p < sysData->getNumPrimaries(); p++){
            // Put primary index # into an array, create event
            double Pix = (double)p;
            tpat_event crashEvt(sysData, tpat_event::CRASH, 0, true, &Pix);
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
tpat_traj tpat_simulation_engine::getTraj() const {
    // Make a copy and return that
    tpat_traj temp (*traj);
    return temp;
}//==============================================================

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
tpat_traj_cr3bp tpat_simulation_engine::getCR3BP_Traj() const{
    if(sysData->getType() == tpat_sys_data::CR3BP_SYS){
        /* Use a static cast to convert the general trajectory pointer
         * into a specific CR3BP trajectory pointer, then dereference
         * and return a COPY of the trajectory
         */
        tpat_traj_cr3bp temp( *(static_cast<tpat_traj_cr3bp *>(traj) ) );
        return temp;
    }else{
        throw tpat_exception("tpat_simulation_engine::getCR3BP_Traj: Wrong system type");
    }
}//==============================================

/**
 *  @brief Retrieve the CR3BP LTVP trajectory. 
 *
 *  To avoid static casts in driver programs,
 *  we create several different getTraj() type functions that will perform
 *  the static cast and return the specific type of trajectory object rather
 *  than a generic one.
 *
 *  @return a CR3BP LTVP Trajectory object
 */
tpat_traj_cr3bp_ltvp tpat_simulation_engine::getCR3BP_LTVP_Traj() const{
    if(sysData->getType() == tpat_sys_data::CR3BP_LTVP_SYS){
        /* Use a static cast to convert the general trajectory pointer
         * into a specific CR3BP trajectory pointer, then dereference
         * and return a COPY of the trajectory
         */
        tpat_traj_cr3bp_ltvp temp( *(static_cast<tpat_traj_cr3bp_ltvp *>(traj) ) );
        return temp;
    }else{
        throw tpat_exception("tpat_simulation_engine::getCR3BP_LTVP_Traj: Wrong system type");
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
tpat_traj_bcr4bpr tpat_simulation_engine::getBCR4BPR_Traj() const{
    if(sysData->getType() == tpat_sys_data::BCR4BPR_SYS){
        // Make a copy and return it
        tpat_traj_bcr4bpr temp( *( static_cast<tpat_traj_bcr4bpr *>(traj) ) );
        return temp;
    }
    else{
        throw tpat_exception("tpat_simulation_engine::getBCR4BPR_Traj: Wrong system type");
    }
}//=====================================

/**
 *  @brief Retrieve a vector of all events being watched for the current simulation
 *  @return a vector of events
 */
std::vector<tpat_event> tpat_simulation_engine::getEvents() const { return events; }

/**
 *  @brief Retrieve a vector of all events that occured during the 
 *  most recent simulation
 *  @return a vector of event records; the event indices correspond to
 *  the indices of the events stored in this simulation engine
 *  @see getEvents()
 */
std::vector<eventRecord> tpat_simulation_engine::getEventRecords() const { return eventOccurs; }

/**
 *  @brief Retrieve a list of all events that fired at the last step
 *  of the simulation, potentially ending the run.
 *  @return a vector of events
 */
std::vector<tpat_event> tpat_simulation_engine::getEndEvents() const{
    if(traj != NULL && traj != 0){
        std::vector<tpat_event> endEvents;
        for(size_t i = 0; i < eventOccurs.size(); i++){
            if(eventOccurs[i].stepIx == (traj->getLength() - 1)){
                endEvents.push_back(events[eventOccurs[i].eventIx]);
            }
        }
        return endEvents;
    }else{
        throw tpat_exception("tpat_simulation_engine::getEndEvents: No trajectory object exists. Please run a sim before querying end events");
    }
}//============================================

/**
 *  @brief Add an event for this integration
 *  @param type the event type
 *  @param dir -1 for negative direction, +1 for positive, 0 for both
 *  @param stop whether or not to stop integration when this event occurs
 */
void tpat_simulation_engine::addEvent(tpat_event::event_t type, int dir, bool stop){
    tpat_event temp(sysData, type, dir, stop);
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
void tpat_simulation_engine::setSysData(tpat_sys_data *d){ sysData = d; }

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
void tpat_simulation_engine::setAbsTol(double t){
    absTol = t;
    if(absTol > 1)
        printWarn("tpat_simulation_engine::setAbsTol: tolerance is greater than 1... just FYI\n");
}

/**
 *	@brief Specify the absolute integration tolerance, non-dimensional units
 *	The default value is 1e-14
 *	@param t the tolerance
 */
void tpat_simulation_engine::setRelTol(double t){
    relTol = t;
    if(relTol > 1)
        printWarn("tpat_simulation_engine::setAbsTol: tolerance is greater than 1... just FYI\n");
}

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
 *  Only the absolute value of the TOF is considered; to integrate backwards in
 *  time, use the setRevTime() function.
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
 *  Only the absolute value of the TOF is considered; to integrate backwards in
 *  time, use the setRevTime() function.
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
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *  Only the absolute value of the TOF is considered; to integrate backwards in
 *  time, use the setRevTime() function.
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
 *  Only the absolute value of the TOF is considered; to integrate backwards in
 *  time, use the setRevTime() function.
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

	switch(sysData->getType()){
		case tpat_sys_data::CR3BP_SYS:
		{
            // Initialize trajectory
            printVerb(verbose, "  initializing CR3BP trajectory\n");
            tpat_sys_data_cr3bp* data = static_cast<tpat_sys_data_cr3bp *>(sysData);
            eomParams = data;
            traj = new tpat_traj_cr3bp(data);
			break;
        }
        case tpat_sys_data::CR3BP_LTVP_SYS:
        {
            printVerb(verbose, "  initializing CR3BP LTVP trajectory\n");
            tpat_sys_data_cr3bp_ltvp* data = static_cast<tpat_sys_data_cr3bp_ltvp *>(sysData);
            eomParams = data;
            traj = new tpat_traj_cr3bp_ltvp(data);
            break;
        }
		case tpat_sys_data::BCR4BPR_SYS:
        {
            printVerb(verbose, "  initializing BCR4BPR trajectory\n");
            tpat_sys_data_bcr4bpr* data = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
            eomParams = data;
            traj = new tpat_traj_bcr4bpr(data);
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
 *  is saved to object-wide storage vectors via the <tt>sim_saveIntegratedData()</tt> function.
 *
 *  @param ic a 6-element initial state for the trajectory
 *  @param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  @param t_dim the dimension of t
 */
void tpat_simulation_engine::integrate(double ic[], double t[], int t_dim){
    // Save tolerance for trajectory
    traj->setTol(absTol > relTol ? absTol : relTol);
    tpat_model *model = sysData->getModel();

    // Get the dimension of the state vector for integration
    int core = model->getCoreStateSize();
    int ic_dim = core + (!simpleIntegration)*(model->getSTMStateSize() + model->getExtraStateSize());
    printVerb(verbose, "  IC has %d initial states\n", ic_dim);

    // Construct the full IC from the state ICs plus the STM ICs and any other ICs for more complex systems
    std::vector<double> fullIC(ic_dim, 0);
    std::copy(ic, ic + core, &(fullIC.front()));

    // ASSUMPTION: STM follows immediately after core states; any extras come after the STM
    if(!simpleIntegration){
        fullIC.at(core + 0) = 1;       // STM initial condition: 6x6 identity matrix
        fullIC.at(core + 7) = 1;
        fullIC.at(core + 14) = 1;
        fullIC.at(core + 21) = 1;
        fullIC.at(core + 28) = 1;
        fullIC.at(core + 35) = 1;      // Elements 42-48 for BCR4BPR ICs are all zero
    }

    double *y = &(fullIC.front());      // array of states that is passed to the integrator

    // Choose EOM function based on system type and simplicity
    printVerb(verbose, "  using %s integration\n", simpleIntegration ? "simple (no STM)" : "full (+ STM)");
    int (*eomFcn)(double, const double[], double[], void*) = 
        simpleIntegration ? model->getSimpleEOM_fcn() : model->getFullEOM_fcn();     // Pointer for the EOM function

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
    model->sim_saveIntegratedData(y, t[0], traj);

    // Update all event functions with IC
    printVerb(verbose, "  sim will use %d event functions:\n", ((int)events.size()));
    for(int ev = 0; ev < ((int)events.size()); ev++){
        printVerb(verbose, "  >>%s\n", events.at(ev).getTypeStr());
        events.at(ev).updateDist(y, t[0]);
    }

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
                printErr("tpat_simulation_engine::integrate: t = %.4e, GSL ERR: %s\n", t0, gsl_strerror(status));
                throw tpat_diverge("Integration did not succeed");
            }

            killSim = locateEvents(y, t0);

            if(killSim)
                break;

            // Put newly integrated state and time into state vector
            model->sim_saveIntegratedData(y, t0, traj);
        }
    }else{
        // Integrate each segment between the input times
        for (int j = 0; j < t_dim - 1; j++){
            // define start and end times; t0 will be updated by integrator
            double t0 = t[j], tf = t[j+1];
            double dt = tf > t0 ? dtGuess : -dtGuess;
            int sgn = tf > t0 ? 1 : -1;

            while(sgn*t0 < sgn*tf && !killSim){
                // printf("Integrating at t = %6.4f\n", t0);
                if(varStepSize){
                    status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t0, tf, &dt, y);
                }else{
                    status = gsl_odeiv2_driver_apply(d, &t0, tf, y);
                }

                if(status != GSL_SUCCESS){
                    printErr("tpat_simulation_engine::integrate: t = %.4e, GSL ERR: %s\n", t0, gsl_strerror(status));
                    throw tpat_diverge("Integration did not succeed");
                }

                killSim = locateEvents(y, t0);
            }

            if(killSim)
                break;

            // Add the newly integrated state and current time fo the state vector
            model->sim_saveIntegratedData(y, t0, traj);
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
    printVerbColor(verbose, GREEN, "  **Integration complete**\n  Total: %d data points\n", traj->getLength()-1);

    // Summarize event occurrences
    for(size_t i = 0; i < eventOccurs.size(); i++){
        printVerb(verbose, " Event %d (%s) occured at step %d\n", eventOccurs[i].eventIx,
            events[eventOccurs[i].eventIx].getTypeStr(), eventOccurs[i].stepIx);
    }
}//===============================================END of cr3bp_integrate

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
 */
bool tpat_simulation_engine::locateEvents(const double *y, double t){
    int numPts = traj->getLength();
    tpat_model *model = sysData->getModel();
    
    // Look through all events
    for(int ev = 0; ev < ((int)events.size()); ev++){
        // Don't trigger if only two points have been integrated
        if(events.at(ev).crossedEvent(y, t) && numPts > 1){

            printVerb(verbose, "  Event %d detected at step %d; searching for exact crossing\n", ev, numPts - 1);
            events.at(ev).incrementCount();  // Update the counter for the event

            if(verbose){ events.at(ev).printStatus(); }

            // Create a nodeset from the previous state (stored in the event) and
            // integrating forwards for half the time between this state and the last one
            double t0 = traj->getTime(-2);          // Time from the state before last
            double ti = traj->getTime(-1);          // Time from the previous state
            double tof = t - t0 - 0.5*(t - ti);     // Approx. TOF 

            // Copy IC into vector - Use the state from two iterations ago to avoid
            // numerical problems when the previous state is REALLY close to the event
            std::vector<double> generalIC = traj->getState(-2);

            if(verbose){
                printColor(BLUE, "Step index = %d\n", numPts-1);
                printColor(BLUE, "t(now) = %f\nt(prev) = %f\nt(prev-1) = %f\n", t, 
                    traj->getTime(-1), traj->getTime(-2));
                printColor(BLUE, "State(now) = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", y[0],
                    y[1], y[2], y[3], y[4], y[5]);
                printColor(BLUE, "tof = %f\n", tof);
                printColor(BLUE, "ic = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", generalIC[0],
                    generalIC[1], generalIC[2], generalIC[3], generalIC[4], generalIC[5]);
            }   

            // Use correction to locate the event very accurately
            if(model->sim_locateEvent(events.at(ev), traj, &(generalIC[0]), t0, tof, verbose)){
                // Remember that this event has occured; step # is one less than the current size
                // of the trajectory
                int timeSize = traj->getLength();
                eventRecord rec(ev, timeSize - 1);
                eventOccurs.push_back(rec);

                // Update event state
                std::vector<double> state = traj->getState(-1);
                double lastT = traj->getTime(-1);
                events.at(ev).updateDist(&(state[0]), lastT);
                
                if(events.at(ev).stopOnEvent() && events.at(ev).getTriggerCount() >= events.at(ev).getStopCount()){
                    printVerbColor(verbose, GREEN, "**Completed Event Location, ending integration**\n");
                    // No need to remember the most recent point; it will be discarded, leaving
                    // the point from mult. shooting as the last
                    return true;    // Tell the simulation to stop
                }else{
                    printVerbColor(verbose, GREEN, "**Completed Event Location, continuing integration**\n");
                    events.at(ev).updateDist(y, t); // Remember the most recent point
                    return false;
                }
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
    delete traj;    // de-allocate the memory
    traj = 0;       // set pointer to 0 (null pointer)
    eomParams = 0;  // set pointer to 0 (null pointer)

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
    madeCrashEvents = false;
}//==========================================