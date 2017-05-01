/**
 *  \file SimEngine.cpp
 *  \brief Performs numerical integration on a set of initial conditions in 
 *  any dynamic model and system
 *  
 *  \author Andrew Cox
 *  \version May 25, 2016
 *  \copyright GNU GPL v3.0
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

#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <ctime>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <vector>

#include "SimEngine.hpp"

#include "AsciiOutput.hpp"
#include "Node.hpp"
#include "SysData_bc4bp.hpp"
#include "Traj_bc4bp.hpp"
#include "BodyData.hpp"
#include "Calculations.hpp"
#include "Common.hpp"
#include "MultShootEngine.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "DynamicsModel.hpp"
#include "Utilities.hpp"

namespace boostInt = boost::numeric::odeint;

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  \brief Construct a simulation engine for a specific dynamical system
 */
SimEngine::SimEngine(){
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Created Simulation Engine\n");
}//===========================================

/**
 *  \brief Copy constructor
 *  \param s a simulation engine 
 */
SimEngine::SimEngine(const SimEngine& s){
    copyMe(s);
}//=====================================

/**
 *  \brief Default destructor
 */
SimEngine::~SimEngine(){
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Destroying simulation engine...\n");
}//===========================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *  \brief Assignment operator; make this engine equal another by copying its data
 *  \param s another simulation engine
 */
SimEngine& SimEngine::operator =(const SimEngine& s){
    copyMe(s);
    return *this;
}//=====================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *  \return whether or not the simulation will leverage simple integration, i.e.,
 *  integration without propagating the STM or other extra parameters - just the
 *  states.
 */
bool SimEngine::usesSimpleInt() const { return bSimpleIntegration; }

/**
 *	\return whether or not the simulation will run time in reverse
 */
bool SimEngine::usesRevTime() const {return bRevTime;}

/**
 *  \return whether or not the engine uses variable step size
 */
bool SimEngine::usesVarStepSize() const { return bVarStepSize; }

/**
 *	\return the absolute tolerance for the engine, non-dimensional units
 */
double SimEngine::getAbsTol() const {return absTol;}

/**
 *  \return ID of the control law implemented during this simulation
 */
unsigned int SimEngine::getCtrlLaw() const {return ctrlLawID; }

/**
 *	\return the relative tolerance for the engine, non-dimensional units
 */
double SimEngine::getRelTol() const {return relTol;}

/**
 *  \return the number of steps the integrator will be forced to take.
 *  The integrator may take intermediate steps between those enforced
 *  by the algorithm, but only <tt>numSteps</tt> data points will be output.
 */
int SimEngine::getNumSteps() const { return numSteps; }

/**
 *  \brief Retrieve a vector of all events being watched for the current simulation
 *  \return a vector of events
 */
std::vector<Event> SimEngine::getEvents() const { return events; }

/**
 *  \brief Add an event for this integration
 *  \param evt an event
 *  \return the index of the event within the events vector; returns -1 
 *  if the event is a duplicate and has not been added to the vector.
 */
int SimEngine::addEvent(Event evt){
    // Make sure this event hasn't been added before
    for(unsigned int e = 0; e < events.size(); e++){
        if(events[e] == evt){
            printErr("SimEngine::addEvent: Event has already been added to the engine; ignoring this new event\n");
            return -1;
        }
    }

    // New event? Great, save it to the list
    events.push_back(evt);
    return events.size()-1;
}//======================================

/**
 *  \brief Determine whether or not default events are created for each simulation
 *  \return whether or not default events are created for each simulation
 */
bool SimEngine::makesDefaultEvents() const { return bMakeDefaultEvents; }

/**
 *	\brief Specify whether or not the engine should run in reverse time
 *	\param b whether or not the engine should run in reverse time
 */
void SimEngine::setRevTime(bool b){ bRevTime = b; }

/**
 *  \brief Specify whether or not the engine should use variable step size.
 *  \param b whether or not the engine should use variable step size
 */
void SimEngine::setVarStepSize(bool b){ bVarStepSize = b; }

/**
 *	\brief Specify the absolute integration tolerance, non-dimensional units.
 *	The default value is 1e-12
 *	\param t the tolerance
 */
void SimEngine::setAbsTol(double t){
    absTol = t;
    if(absTol > 1)
        astrohelion::printWarn("SimEngine::setAbsTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *	\brief Specify the relative integration tolerance, non-dimensional units
 *	The default value is 1e-14
 *	\param t the tolerance
 */
void SimEngine::setRelTol(double t){
    relTol = t;
    if(relTol > 1)
        astrohelion::printWarn("SimEngine::setRelTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *  \brief Specify the ID of the control law to be implemented during the simulation
 *  \param id control law ID
 */
void SimEngine::setCtrlLaw(unsigned int id){ ctrlLawID = id; }

/**
 *  \brief Tell the simulation engine whether or not to make default events at the
 *  beginning of the simulation
 * 
 *  \param b whether or not to create default events
 */
void SimEngine::setMakeDefaultEvents(bool b){ bMakeDefaultEvents = b; }

/**
 *  \brief Set the maximum computation time limit (seconds).
 *  \details The numerical integration is stopped once the specified number
 *  of seconds have ellapsed since the beginning of the integration. 
 * 
 *  \param t maximum allowable seconds for the numerical integration.
 */
void SimEngine::setMaxCompTime(int t){ maxCompTime = t; }

/**
 *  \brief Specify the number of steps the integrator must take during the 
 *  the integration. 
 *
 *  Only these points will be output to the 
 *  trajectory object, although the GSL driver may take steps in between
 *  those specified to maintain numerical accuracy.
 */
void SimEngine::setNumSteps(int n){ numSteps = n; }

/**
 *  \brief Tell the engine whether or not to use simple integration (i.e., no STM propagation)
 *  \param b whether or not to use simple integration
 */
void SimEngine::setSimpleInt(bool b){ bSimpleIntegration = b; }

/**
 *  \brief Set the variable step integration step function
 *  \param integ The type of integrator to use for variable step propagations
 */
void SimEngine::setVarStepInteg(Integ_tp integ){ 
    switch(integ){
        case Integ_tp::RKCK:
        case Integ_tp::RK8PD:
        case Integ_tp::BOOST_RKCK:
        case Integ_tp::BOOST_RKF:
        case Integ_tp::BOOST_BS:
            varStep_integ = integ;
            break;
        default:
            printErr("SimEngine::SetVarStepInteg: Integrator type is not suited for variable step propagation; ignoring change\n");
    }
}//====================================================

/**
 *  \brief Set the fixed step integration step function
 *  \param integ The type of integrator to use for fixed step propagations
 */
void SimEngine::setFixStepInteg(Integ_tp integ){
    switch(integ){
        case Integ_tp::MSADAMS:
            fixStep_integ = integ;
            break;
        default:
            printErr("SimEngine::SetFixStepInteg: Integrator type is not suited for fixed step propagation; ignoring change\n");
    }
}//====================================================

//-----------------------------------------------------
// 		Simulation Functions
//-----------------------------------------------------

/**
 *	\brief Run a simulation given a set of initial conditions and run time. 
 *
 *  It is assumed that t0 = 0
 *	\param ic an array containting the non-dimensional initial state; number
 *    of elements must match the number of <code>coreStates</code> specified 
 *    in the system dynamics model
 *	\param tof the total integration time, or time-of-flight (non-dim units)
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  \param traj pointer to a trajectory object to store the output trajectory
 */
void SimEngine::runSim(const double *ic, double tof, Traj *traj){
	runSim(ic, 0, tof, traj);
}//=======================================================

/**
 *  \brief Run a simulation given a set of initial conditions and run time
 *
 *  It is assumed that t0 = 0
 *  \param ic a vector containing the IC (nondimensional states)
 *  \param tof the total integration time, or time-of-flight (non-dim units)
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  \param traj pointer to a trajectory object to store the output trajectory
 */
void SimEngine::runSim(std::vector<double> ic, double tof, Traj *traj){
    runSim(ic, 0, tof, traj);
}//=======================================================

/**
 *  \brief Run a simulation given a set of initial conditions and run time
 *
 *  It is assumed that t0 = 0
 *  \param ic a vector containing the IC (nondimensional states)
 *  \param t0 the epoch time associated with the IC (non-dim units)
 *  \param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  \param traj pointer to a trajectory object to store the output trajectory
 *  \throws Exception if <tt>ic</tt> has fewer than 6 elements
 *  \throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Traj object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(std::vector<double> ic, double t0, double tof, Traj *traj){
    if(ic.size() >= static_cast<size_t>(traj->getSysData()->getDynamicsModel()->getCoreStateSize())){
        std::vector<double> tempIC = ic;
        runSim(&(tempIC[0]), t0, tof, traj);
    }else{
        printErr("IC size = %zu\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by coreStates in the Dynamics Model");
    }
}//=======================================================

/**
 *  \brief Run a simulation given a set of initial conditions and run time
 *
 *  It is assumed that t0 = 0
 *  \param ic a vector containing the IC (nondimensional states)
 *  \param stm0 initial STM
 *  \param t0 the epoch time associated with the IC (non-dim units)
 *  \param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  \param traj pointer to a trajectory object to store the output trajectory
 *  \throws Exception if <tt>ic</tt> has fewer than 6 elements
 *  \throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Traj object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(std::vector<double> ic, MatrixXRd stm0, double t0, double tof, Traj *traj){
    if(ic.size() >= static_cast<size_t>(traj->getSysData()->getDynamicsModel()->getCoreStateSize())){
        std::vector<double> tempIC = ic;
        runSim(&(tempIC[0]), stm0, t0, tof, traj);
    }else{
        printErr("IC size = %zu\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by coreStates in the Dynamics Model");
    }
}//=======================================================

/**
 *	\brief Run a simulation in the specified system starting with a set of initial conditions,
 *  at a specified initial time, and integrating for a specified time-of-flight
 *	\param ic an array of non-dimensional initial states; number
 *    of elements must match the number of <code>coreStates</code> specified 
 *    in the system dynamics model
 *	\param t0 the time at the start of the integration, non-dimensional units
 *	\param tof time-of-flight, non-dimensional time units
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  \param traj pointer to a trajectory object to store the output trajectory
 */
void SimEngine::runSim(const double *ic, double t0, double tof, Traj *traj){
    std::vector<double> t_span;
    // Compute the final time based on whether or not we're using reverse time integration
    double tf = bRevTime ? t0 - std::abs(tof) : t0 + std::abs(tof);
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  time will span from %.3e to %.3e\n", t0, tf);
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  (Reverse Time is %s)\n", bRevTime ? "ON" : "OFF");

    if(bVarStepSize){
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

    int core_dim = traj->getSysData()->getDynamicsModel()->getCoreStateSize();
    MatrixXRd stm0 = MatrixXRd::Identity(core_dim, core_dim);
    runSim(ic, stm0, t_span, traj);
}//====================================================

void SimEngine::runSim(const double *ic, MatrixXRd stm0, double t0, double tof, Traj *traj){
    std::vector<double> t_span;
    // Compute the final time based on whether or not we're using reverse time integration
    double tf = bRevTime ? t0 - std::abs(tof) : t0 + std::abs(tof);
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  time will span from %.3e to %.3e\n", t0, tf);
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  (Reverse Time is %s)\n", bRevTime ? "ON" : "OFF");

    if(bVarStepSize){
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

    runSim(ic, stm0, t_span, traj);
}//====================================================

/**
 *  \brief Run a simulation in the specified system starting with a set of initial conditions,
 *  at a specified initial time, and integrating for a specified time-of-flight
 *  \param ic an array of non-dimensional initial states; number
 *    of elements must match the number of <code>coreStates</code> specified 
 *    in the system dynamics model
 *  \param stm0 initial STM
 *  \param t_span a vector of times to include in the solution.
 *  \param traj pointer to a trajectory object to store the output trajectory
 */
void SimEngine::runSim(const double *ic, MatrixXRd stm0, std::vector<double> t_span, Traj *traj){
    astrohelion::printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "Running simulation...\n");
    if(!bIsClean){
        cleanEngine();
    }

    if(bMakeDefaultEvents)
        createDefaultEvents(traj->getSysData());

    eomParams = new EOM_ParamStruct(traj->getSysData(), ctrlLawID);

    // Run the simulation
    bIsClean = false;   // Technically, nothing has changed yet, but this flag should be false even if any part of integrate throws an exception
    integrate(ic, stm0, &(t_span.front()), t_span.size(), traj);
}//====================================================

//-----------------------------------------------------
// 		Numerical Integration
//-----------------------------------------------------

/**
 *  \brief Integrate the state EOMs and STM EOMs with additional integration as required by 
 *  specific systems.
 *
 *  This function uses values stored in object-wide variables to determine the direction time flows,
 *  whether or not to use simple integration, and whether or not to use variable step sizes.
 *
 *  \param ic an array containing the initial state for the trajectory; number
 *    of elements must match the number of <code>coreStates</code> specified 
 *    in the system dynamics model
 *  \param stm0 initial STM
 *  \param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  \param t_dim the dimension of t
 *  \param traj pointer to a trajectory object to store the output trajectory
 *  
 *  \throws DivergeException if the integrator fails and cannot proceed
 */
void SimEngine::integrate(const double *ic, MatrixXRd stm0, const double *t, int t_dim, Traj *traj){
    // delete any auto-generated events from previous integrations
    for(unsigned int i = 0; i < events.size(); i++){
        if(static_cast<int>(events[i].getType()) <= 0){
            events.erase(events.begin() + i);
            i--;
        }
    }

    startTimestamp = time(nullptr);

    // Save tolerance for trajectory
    traj->setTol(absTol > relTol ? absTol : relTol);
    const DynamicsModel *model = traj->getSysData()->getDynamicsModel();

    // Initialize all events with the correct system data pointer
    for(unsigned int i = 0; i < events.size(); i++){
        events[i].initialize(traj->getSysData());
    }

    // Get the dimension of the state vector for integration
    unsigned int core_dim = model->getCoreStateSize();
    unsigned int ic_dim = core_dim + (!bSimpleIntegration)*(core_dim*core_dim + model->getExtraStateSize());
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  IC has %u initial states\n", ic_dim);

    // Construct the full IC from the state ICs plus the STM ICs and any other ICs for more complex systems
    std::vector<double> fullIC(ic_dim, 0);
    std::copy(ic, ic + core_dim, &(fullIC.front()));

    if(stm0.rows() != core_dim || stm0.cols() != core_dim){
        printErr("STM rows = %d, cols = %d, core_dim = %u\n", stm0.rows(), stm0.cols(), core_dim);
        throw Exception("SimEngine::integrate: Initial STM size does not match the core state size specified by the Dynamic Model");
    }

    // ASSUMPTION: STM follows immediately after core states; any extras come after the STM
    if(!bSimpleIntegration){
        // Add STM to initial conditions
        for(unsigned int r = 0; r < static_cast<unsigned int>(core_dim); r++){
            for(unsigned int c = 0; c < static_cast<unsigned int>(core_dim); c++){
                fullIC.at(core_dim + r*core_dim + c) = stm0(r,c);
            }
        }
    }

    double *y = &(fullIC.front());      // array of states that is passed to the integrator

    // Choose EOM function based on system type and simplicity
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  using %s integration\n", bSimpleIntegration ? "simple (no STM)" : "full (+ STM)");
    int (*eomFcn)(double, const double[], double[], void*) = 
        bSimpleIntegration ? model->getSimpleEOM_fcn() : model->getFullEOM_fcn();     // Pointer for the EOM function

    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, 
        "  using control law: %s\n", traj->getSysData()->getControlLaw()->lawIDToString(ctrlLawID).c_str());

    /*
     * BOOST INTEGRATOR ADDITION
     */
    // if(bVarStepSize && (
    //     varStep_integ == Integ_tp::BOOST_RKCK ||
    //     varStep_integ == Integ_tp::BOOST_RKF ||
    //     varStep_integ == Integ_tp::BOOST_BS) ){

    //     if(t_dim != 2)
    //         throw Exception("SimEngine::integrate: t_dim must be 2 for BOOST integration right now!");

    //     boost_eom_wrapper boostEOM(eomFcn, eomParams);
    //     boost_observer boostObserver(model, traj, eomParams);

    //     double t0 = t[0], tf = t[1];                // start and finish times for integration; t0 will be updated by integrator
    //     double dt = tf > t0 ? dtGuess : -dtGuess;   // step size (initial guess)

    //     // Save the initial state, time, and STM
    //     // model->sim_saveIntegratedData(y, t[0], traj, eomParams);

    //     switch(varStep_integ){
    //         case Integ_tp::BOOST_RKCK:
    //             boostInt::integrate_adaptive(
    //                 boostInt::make_controlled< boostInt::runge_kutta_cash_karp54< std::vector<double> > >(absTol, relTol),
    //                 boostEOM, fullIC, t0, tf, dt, boostObserver);
    //             break;
    //         case Integ_tp::BOOST_RKF:
    //             boostInt::integrate_adaptive(
    //                 boostInt::make_controlled< boostInt::runge_kutta_fehlberg78< std::vector<double> > >(absTol, relTol),
    //                 boostEOM, fullIC, t0, tf, dt, boostObserver);
    //             break;
    //         case Integ_tp::BOOST_BS:
    //         {
    //             boostInt::bulirsch_stoer<std::vector<double> > stepper(absTol, relTol);
    //             integrate_adaptive(stepper, boostEOM, fullIC, t0, tf, dt, boostObserver);
    //             break;
    //         }    
    //     }

    //     return; // Skip the rest of the function and exit
    // }// END OF BOOST INTEGRATORS

    /* Create a system to integrate; we don't include a Jacobin (NULL)
     *  The parameter set eomParams can be modified 
     *  between integration steps (i.e., change model parameters), but the ode functions must be reset
     *  via <code>gsl_odeiv2_driver_reset</code>, <code>gsl_odeiv2_evolve_reset</code>, or
     *  <code>gsl_odeiv2_step_reset</code> before continuing with an updated parameter set
     */
    gsl_odeiv2_system odeSys = {eomFcn, NULL, static_cast<size_t>(ic_dim), eomParams};
    
    // Define ODE objects, define them conditionaly based on bVarStepSize
    gsl_odeiv2_step *s = NULL;
    gsl_odeiv2_control *c = NULL;
    gsl_odeiv2_evolve *e = NULL;
    gsl_odeiv2_driver *d = NULL;

    if(bVarStepSize){
        // Allocate space for the stepping object
        switch(varStep_integ){
            case Integ_tp::RK8PD:
                s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, ic_dim);
                astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  variable step size, using Runge-Kutta Dormand-Prince 8-9 method\n");
                break;
            case Integ_tp::RKCK:
                s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, ic_dim);
                astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  variable step size, using Runge-Kutta Cash-Karp 4-5 method\n");
                break;
            default:
                free_odeiv2(s, c, e, d);
                throw Exception("SimEngine::integrate: Integrator type is not suited for variable step propagation; aborting integration process");
        }
        
        // Define a control that will keep the error in the state y within the specified tolerances
        c = gsl_odeiv2_control_y_new (absTol, relTol);
        // Allocate space for the integrated solution to evolve in
        e = gsl_odeiv2_evolve_alloc(ic_dim);
    }else{
        astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  fixed step size, using Adams-Bashforth, Adams-Moulton method\n");
        
        double signed_dt = bRevTime ? -1*dtGuess : dtGuess;

        switch(fixStep_integ){
            case Integ_tp::MSADAMS:
                // Allocate space for a driver; the msadams algorithm requires access to the driver
                d = gsl_odeiv2_driver_alloc_y_new(&odeSys, gsl_odeiv2_step_msadams, signed_dt, absTol, relTol);
                // Allocate space for the stepping object
                s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_msadams, ic_dim);
                gsl_odeiv2_step_set_driver(s, d);
                break;
            default:
                free_odeiv2(s, c, e, d);
                throw Exception("SimEngine::integrate: Integrator type is not suited for fixed step propagation; aborting integration process");
        }
    }

    // Update all event functions with IC
    astrohelion::printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  sim will use %d event functions:\n", static_cast<int>(events.size()));
    for(unsigned int ev = 0; ev < events.size(); ev++){
        astrohelion::printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  [%02u] - %s\n", ev, events.at(ev).getTypeStr());
        events.at(ev).updateDist(y, core_dim, t[0]);
    }

    // Save the initial state, time, and STM
    Node node0(y, core_dim, t[0]); // Construct a basic node
    int id0 = model->sim_addNode(node0, y, t[0], traj, eomParams, Event_tp::NONE); // Pass the node to the Dynamic Model so specific modifications may be made
    Segment seg;            // Construct a segment that has an origin at the initial node
    seg.setOrigin(id0);     // Assign origin ID
    seg.appendState(y, ic_dim);     // Save the initial state in the segment
    seg.appendTime(t[0]);           // Save the initial time in the segment
    model->sim_addSeg(seg, y, t[0], traj, eomParams);

    int status;             // integrator status
    int propStepCount = 0;  // Count of propagated steps (different than number of time steps)
    bool killSim = false;   // flag to stop simulation
    double t_int = t[0];    // current propagation time

    if(t_dim == 2){
        double tf = t[1];                // start and finish times for integration; t_int will be updated by integrator
        double dt = tf > t_int ? dtGuess : -dtGuess;   // step size (initial guess)
        int sgn = tf > t_int ? 1 : -1;

        while (sgn*t_int < sgn*tf && !killSim){
            // apply integrator for one step; new state and time are stored in y and t_int
            if(bVarStepSize){
                status = gsl_odeiv2_evolve_apply(e, c, s, &odeSys, &t_int, tf, &dt, y);
            }else{
                status = gsl_odeiv2_driver_apply(d, &t_int, tf, y);
            }

            if(status != GSL_SUCCESS){
                if(verbosity >= Verbosity_tp::SOME_MSG){
                    astrohelion::printErr("SimEngine::integrate: t = %.4e, GSL ERR: %s\n", t_int, gsl_strerror(status));
                }

                free_odeiv2(s, c, e, d);

                // Create event and record that will indicate propagation ended with an error
                Event event_err(Event_tp::SIM_ERR, 0, false);

                // Save the last successful step the integrator was able to take
                Node nodeF(y, core_dim, t_int);    // Construct a basic node
                int idf = model->sim_addNode(nodeF, y, t_int, traj, eomParams, event_err.getType()); // Pass the node to the DynamicModel for more specific actions
                // Get reference to most recent segment (faster than copying the object!)
                Segment &lastSeg = traj->getSegRefByIx(-1);
                lastSeg.setTerminus(idf);                       // Update the terminus to the final node
                lastSeg.appendState(y, ic_dim);                 // Save the final state and time to the segment
                lastSeg.appendTime(t_int);
                lastSeg.storeTOF();
                lastSeg.setSTM(y+core_dim, core_dim*core_dim);

                throw DivergeException("SimEngine::integrate: Integration did not succeed");
            }

            // Stop the simulation if a simulation-ending event occurs
            killSim = locateEvents(y, t_int, traj, propStepCount);

            // Stop the simulation if the maximum computation time has passed
            killSim = killSim || (maxCompTime > 0 && (time(nullptr) - startTimestamp) > maxCompTime);

            if(killSim)
                break;

            // Save propagation state to segment
            Segment &lastSeg = traj->getSegRefByIx(-1);     // Get another reference (not reassignable)
            lastSeg = traj->getSegRefByIx(-1);  // Update in case locateEvents() added a node and segment
            lastSeg.appendState(y, ic_dim);
            lastSeg.appendTime(t_int);
            propStepCount++;
        }
    }else{
        // Integrate each segment between the input times
        for (int j = 0; j < t_dim - 1; j++){
            // define start and end times; t_int will be updated by integrator
            double t_int = t[j], tf = t[j+1];
            double dt = tf > t_int ? dtGuess : -dtGuess;
            int sgn = tf > t_int ? 1 : -1;

            if(j > 0){
                // Create a node at every time in the t[] array
                Node nodeI(y, core_dim, t_int);
                int nodeID_i = model->sim_addNode(nodeI, y, t_int, traj, eomParams, Event_tp::SIM_TOF);
                Segment &lastSeg = traj->getSegRefByIx(-1);     // Get reference (not reassignable)
                lastSeg.setTerminus(nodeID_i);  // Update the previous segment to terminate at the new node
                lastSeg.storeTOF();
                lastSeg.setSTM(y+core_dim, core_dim*core_dim);

                // Create a new segment for the next time interval
                Segment newSeg;
                newSeg.setOrigin(nodeID_i);
                newSeg.appendState(y, ic_dim);
                newSeg.appendTime(t_int);
                model->sim_addSeg(newSeg, y, t_int, traj, eomParams);
                lastSeg = traj->getSegRefByIx(-1);  // Update reference to most recent segment
            }
            while(sgn*t_int < sgn*tf && !killSim){
                // printf("Integrating at t = %6.4f\n", t_int);
                if(bVarStepSize){
                    status = gsl_odeiv2_evolve_apply(e, c, s, &odeSys, &t_int, tf, &dt, y);
                }else{
                    status = gsl_odeiv2_driver_apply(d, &t_int, tf, y);
                }

                if(status != GSL_SUCCESS){
                    if(verbosity >= Verbosity_tp::SOME_MSG){
                        astrohelion::printErr("SimEngine::integrate: t = %.4e, GSL ERR: %s\n", t_int, gsl_strerror(status));
                    }
                    
                    free_odeiv2(s, c, e, d);

                    // Save the last successful step the integrator was able to take
                    Node nodeF(y, core_dim, t_int);    // Construct a basic node
                    // Add the node with info that there was a simulation error here
                    int idf = model->sim_addNode(nodeF, y, t_int, traj, eomParams, Event_tp::SIM_ERR);
                    Segment &lastSeg = traj->getSegRefByIx(-1);     // Get another reference (not reassignable)
                    lastSeg.setTerminus(idf);                       // Update the terminus to the final node
                    lastSeg.appendState(y, ic_dim);                 // Save the final state and time to the segment
                    lastSeg.appendTime(t_int);
                    lastSeg.storeTOF();
                    lastSeg.setSTM(y+core_dim, core_dim*core_dim);

                    throw DivergeException("SimEngine::integrate: Integration did not succeed");
                }

                killSim = locateEvents(y, t_int, traj, propStepCount);

                // Stop the simulation if the maximum computation time has passed
                killSim = killSim || (maxCompTime > 0 && (time(nullptr) - startTimestamp) > maxCompTime);

                propStepCount++;

                // Save the most recent time and state to the segment
                Segment &lastSeg = traj->getSegRefByIx(-1);     // Get another reference (not reassignable)
                lastSeg.appendState(y, ic_dim);
                lastSeg.appendTime(t_int);
            }

            if(killSim)
                break;
        }
    }

    // Clean Up - some have been allocated, some have not
    free_odeiv2(s, c, e, d);

    // Trigger events that may have ended the propagation (other than default or user-defined events)
    if(!killSim){
        // Ended normally at the desired TOF
        // Create a final node and update the final segment
        Node nodeF(y, core_dim, t_int);
        Segment &lastSeg = traj->getSegRefByIx(-1);
        int nodeID_f = model->sim_addNode(nodeF, y, t_int, traj, eomParams, Event_tp::SIM_TOF);
        lastSeg.setTerminus(nodeID_f);  // Just update the terminus: final state and time should have already been appended
        lastSeg.storeTOF();
        lastSeg.setSTM(y+core_dim, core_dim*core_dim);
    }else if(maxCompTime > 0 && (time(nullptr) - startTimestamp) > maxCompTime){
        // Ended at the time-out
        // Create a final node and update the final segment
        Node nodeF(y, core_dim, t_int);
        Segment &lastSeg = traj->getSegRefByIx(-1);
        int nodeID_f = model->sim_addNode(nodeF, y, t_int, traj, eomParams, Event_tp::SIM_COMPTIME);
        lastSeg.setTerminus(nodeID_f);  // Just update the terminus: final state and time should have already been appended
        
        // Propagation ended without saving the final state
        lastSeg.appendState(y, ic_dim);
        lastSeg.appendTime(t_int);

        lastSeg.storeTOF();
        lastSeg.setSTM(y+core_dim, core_dim*core_dim);
    }

    // Check lengths of vectors and set the numPoints value in traj
    astrohelion::printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Integration complete**\nTotal: %d Nodes\n", traj->getNumNodes());

    for(unsigned int n = 0; n < traj->getNumNodes(); n++){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, " %s Event occured at node %u\n",
            Event::getEventTpStr(traj->getNodeRefByIx(n).getTriggerEvent()), n);
    }
}//====================================================END of cr3bp_integrate

/**
 *  \brief Free all GSL ODE pointers that have been initialized/instantiated
 * 
 *  \param s stepper object pointer
 *  \param c control object pointer
 *  \param e evolver object pointer
 *  \param d driver object pointer
 */
void SimEngine::free_odeiv2(gsl_odeiv2_step *s, gsl_odeiv2_control *c, gsl_odeiv2_evolve *e, gsl_odeiv2_driver *d){
    if(bVarStepSize){
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c);
    }else{
        gsl_odeiv2_driver_free(d);
    }
    gsl_odeiv2_step_free(s);
}//====================================================

/**
 *  \brief Locate event occurences as exactly as possible and determine if the simulation 
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
 *  \param y the most recent state on the integrated arc.
 *  \param t the time associated with y
 *  \param traj pointer to a trajectory object to store the output trajectory
 *  \param propStepCount the number of steps the propagation has taken so far. This is different
 *  from the number of time steps as many propagation steps may be taken between specified time steps
 *  \return whether or not the simulation should end (an event triggers killSim)
 */
bool SimEngine::locateEvents(const double *y, double t, Traj *traj, int propStepCount){
    const DynamicsModel *model = traj->getSysData()->getDynamicsModel();
    unsigned int core_dim = model->getCoreStateSize();
    unsigned int full_dim = core_dim + (!bSimpleIntegration)*(core_dim*core_dim + model->getExtraStateSize());

    // Look through all events
    for(unsigned int ev = 0; ev < events.size(); ev++){
        // Don't trigger if only two points have been integrated
        // Also don't trigger if the type is 0 or negative: these are managed by the simulation engine
        if(propStepCount > 1 && to_underlying(events[ev].getType()) > 0 && events[ev].crossedEvent(y, core_dim, t)){

            astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG,
                "  Event %d detected on segment %d; searching for exact crossing\n", ev, traj->getNumSegs());
            events[ev].incrementCount();  // Update the counter for the event

            if(verbosity >= Verbosity_tp::ALL_MSG){ events[ev].printStatus(); }

            // Create a nodeset from the previous state (stored in the event) and
            // integrate forwards for half the time between this state and the last one
            Segment &lastSeg = traj->getSegRefByIx(-1);
            double t0 = lastSeg.getTimeByIx(-2);          // Time from the state before last
            double ti = lastSeg.getTimeByIx(-1);          // Time from the previous state
            double tof = t - t0 - 0.5*(t - ti);         // Approx. TOF 

            // Copy IC into vector - Use the state from two iterations ago to avoid
            // numerical problems when the previous state is REALLY close to the event
            std::vector<double> generalIC = lastSeg.getStateByRow(-2, full_dim);

            if(verbosity >= Verbosity_tp::ALL_MSG){
                astrohelion::printColor(BLUE, "Step index = %d\n", propStepCount-1);
                astrohelion::printColor(BLUE, "t(now) = %f\nt(prev) = %f\nt(prev-1) = %f\n", t, ti, t0);
                astrohelion::printColor(BLUE, "State(now) = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", y[0],
                    y[1], y[2], y[3], y[4], y[5]);
                astrohelion::printColor(BLUE, "tof = %f\n", tof);
                astrohelion::printColor(BLUE, "ic = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", generalIC[0],
                    generalIC[1], generalIC[2], generalIC[3], generalIC[4], generalIC[5]);
            }   

            // Use correction to locate the event very accurately
            if(model->sim_locateEvent(events[ev], traj, &(generalIC[0]), t0, tof, eomParams, verbosity)){
                // Update event state from the most recent node (a new node was created at the event occurence)
                std::vector<double> state = traj->getStateByIx(-1);
                double lastT = traj->getTimeByIx(-1);
                events[ev].updateDist(&(state[0]), core_dim, lastT);
                
                // This condition is also checked in model->sim_locateEvent(): that function adds a new
                // segment if the simulation is continuing
                if(events[ev].stopOnEvent() && events[ev].getTriggerCount() >= events[ev].getStopCount()){
                    astrohelion::printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Completed Event Location, ending integration**\n");
                    // No need to remember the most recent point; it will be discarded, leaving
                    // the point from mult. shooting as the last
                    return true;    // Tell the simulation to stop
                }else{
                    astrohelion::printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Completed Event Location, continuing integration**\n");
                    events[ev].updateDist(y, core_dim, t); // Remember the most recent point
                    return false;
                }
            }
        }// end of If(hasCrossed)

        if(to_underlying(events[ev].getType()) > 0){
            // Save the distance and current state to the event
            events[ev].updateDist(y, core_dim, t);
            printVerbColor(verbosity >= Verbosity_tp::DEBUG, BLUE, "Updating event %d state to [%.4f, %.4f, %.4f,...]\n",
                ev, y[0], y[1], y[2]);
        }
    }// end of loop

    return false;
}//========================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  \brief Clean out the trajectory storage variable so a new simulation can be run and store its data
 *  \details This function does not reset any parameters the user has set
 */
void SimEngine::cleanEngine(){
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Cleaning the engine...\n");
    delete eomParams;
    eomParams = nullptr;  // set pointer to 0 (null pointer)

    for(unsigned int e = 0; e < events.size(); e++){
        events[e].reset();
    }

    bIsClean = true;
}//====================================================

/**
 *  Create default events for the system
 *  \param sysData pointer to the system data object for the simulation
 */
void SimEngine::createDefaultEvents(const SysData *sysData){
    if(!bMadeDefaultEvents){
        
        std::vector<Event> defEvents = sysData->getDynamicsModel()->sim_makeDefaultEvents(sysData);
        for(unsigned int i = 0; i < defEvents.size(); i++)
            addEvent(defEvents[i]);

        bMadeDefaultEvents = true;
    }
}//====================================================

/**
 *  \brief Reset all variables and options
 *
 *  Completely resets the simulation engine, reverting all variables (including ones the user
 *  has modified with set() functions) to their default values.
 */
void SimEngine::reset(){
    if(!bIsClean)
        cleanEngine();

    events.clear();
    bRevTime = false;
    verbosity = Verbosity_tp::NO_MSG;
    bVarStepSize = true;
    absTol = 1e-12;
    relTol = 1e-14;
    dtGuess = 1e-6;
    numSteps = 1000;
    bMadeDefaultEvents = false;
    bMakeDefaultEvents = true;
    maxCompTime = -1;
    ctrlLawID = 0;
}//====================================================

/**
 *  Clear all events from the simulation, including any created by default.
 */
void SimEngine::clearEvents(){
    astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Clearing all events...\n");
    events.clear();
    bMadeDefaultEvents = false;
}//====================================================

/**
 *  Copy data from an input engine to this one
 *  \param s an input simulation engine
 *  @throw Exception if <tt>s</tt> has an unknown system data type
 */
void SimEngine::copyMe(const SimEngine &s){
    Engine::copyBaseEngine(s);
    delete eomParams;
    eomParams = nullptr;  // Will get set again by the runSim() method
    bRevTime = s.bRevTime;
    verbosity = s.verbosity;
    bVarStepSize = s.bVarStepSize;
    bSimpleIntegration = s.bSimpleIntegration;
    absTol = s.absTol;
    relTol = s.relTol;
    dtGuess = s.dtGuess;
    numSteps = s.numSteps;
    events = s.events;
    bMakeDefaultEvents = s.bMakeDefaultEvents;
    bMadeDefaultEvents = s.bMadeDefaultEvents;
    maxCompTime = s.maxCompTime;
    startTimestamp = s.startTimestamp;
    ctrlLawID = s.ctrlLawID;
}//====================================================

}// END of Astrohelion namespace