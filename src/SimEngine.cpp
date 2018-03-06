/**
 *  @file SimEngine.cpp
 *  @brief Performs numerical integration on a set of initial conditions in 
 *  any dynamic model and system
 *  
 *  @author Andrew Cox
 *  @version May 25, 2016
 *  @copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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
#include "Arcset.hpp"
#include "BodyData.hpp"
#include "Calculations.hpp"
#include "Constraint.hpp"
#include "Common.hpp"
#include "DynamicsModel.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "MultShootEngine.hpp"
#include "Node.hpp"
#include "Utilities.hpp"

namespace boostInt = boost::numeric::odeint;

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Construct a simulation engine for a specific dynamical system
 */
SimEngine::SimEngine(){
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Created Simulation Engine\n");
}//===========================================

/**
 *  @brief Copy constructor
 *  @param s a simulation engine 
 */
SimEngine::SimEngine(const SimEngine& s){
    copyMe(s);
}//=====================================

/**
 *  @brief Default destructor
 */
SimEngine::~SimEngine(){
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Destroying simulation engine...\n");
    if(eomParams)
        delete(eomParams);
}//===========================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *  @brief Assignment operator; make this engine equal another by copying its data
 *  @param s another simulation engine
 */
SimEngine& SimEngine::operator =(const SimEngine& s){
    copyMe(s);
    return *this;
}//=====================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *  @return whether or not the simulation will leverage simple integration, i.e.,
 *  integration without propagating the STM or other extra parameters - just the
 *  states.
 */
bool SimEngine::usesSimpleInt() const { return bSimpleIntegration; }

/**
 *	@return whether or not the simulation will run time in reverse
 */
bool SimEngine::usesRevTime() const {return bRevTime;}

/**
 *  @return whether or not the engine uses variable step size
 */
bool SimEngine::usesVarStepSize() const { return bVarStepSize; }

/**
 *	@return the absolute tolerance for the engine, non-dimensional units
 */
double SimEngine::getAbsTol() const {return absTol;}

/**
 *	@return the relative tolerance for the engine, non-dimensional units
 */
double SimEngine::getRelTol() const {return relTol;}

/**
 *  @return the number of steps the integrator will be forced to take.
 *  The integrator may take intermediate steps between those enforced
 *  by the algorithm, but only `numSteps` data points will be output.
 */
int SimEngine::getNumSteps() const { return numSteps; }

/**
 *  @brief Retrieve a vector of all events being watched for the current simulation
 *  @return a vector of events
 */
std::vector<Event> SimEngine::getEvents() const { return events; }

/**
 *  @brief Add an event for this integration
 *  @param evt an event
 *  @return the index of the event within the events vector; returns -1 
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
 *  @brief Determine whether or not default events are created for each simulation
 *  @return whether or not default events are created for each simulation
 */
bool SimEngine::makesDefaultEvents() const { return bMakeDefaultEvents; }

/**
 *	@brief Specify whether or not the engine should run in reverse time
 *	@param b whether or not the engine should run in reverse time
 */
void SimEngine::setRevTime(bool b){ bRevTime = b; }

/**
 *  @brief Specify whether or not the engine should use variable step size.
 *  @param b whether or not the engine should use variable step size
 */
void SimEngine::setVarStepSize(bool b){ bVarStepSize = b; }

/**
 *	@brief Specify the absolute integration tolerance, non-dimensional units.
 *	The default value is 1e-12
 *	@param t the tolerance
 */
void SimEngine::setAbsTol(double t){
    absTol = t;
    if(absTol > 1)
        printWarn("SimEngine::setAbsTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *	@brief Specify the relative integration tolerance, non-dimensional units
 *	The default value is 1e-14
 *	@param t the tolerance
 */
void SimEngine::setRelTol(double t){
    relTol = t;
    if(relTol > 1)
        printWarn("SimEngine::setRelTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *  @brief Tell the simulation engine whether or not to make default events at the
 *  beginning of the simulation
 * 
 *  @param b whether or not to create default events
 */
void SimEngine::setMakeDefaultEvents(bool b){ bMakeDefaultEvents = b; }

/**
 *  @brief Set the maximum computation time limit (seconds).
 *  @details The numerical integration is stopped once the specified number
 *  of seconds have ellapsed since the beginning of the integration. 
 * 
 *  @param t maximum allowable seconds for the numerical integration.
 */
void SimEngine::setMaxCompTime(double t){ maxCompTime = t; }

/**
 *  @brief Specify the number of steps the integrator must take during the 
 *  the integration. 
 *
 *  @details Only these points will be output to the 
 *  trajectory object, although the GSL driver may take steps in between
 *  those specified to maintain numerical accuracy.
 *  
 *  @param n the number of steps. Must be at least 2. If n < 2, the
 *  simulation will proceed but will use 2 steps rather than the invalid
 *  number specified.
 */
void SimEngine::setNumSteps(int n){
    if(n < 2){
        numSteps = 2;
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::setNumSteps: Must use at least 2 steps; input number is too few, using numSteps = 2 and continuing.\n");
    }
    else{
        numSteps = n;
    }
}//====================================================

/**
 *  @brief Tell the engine whether or not to use simple integration (i.e., no STM propagation)
 *  @param b whether or not to use simple integration
 */
void SimEngine::setSimpleInt(bool b){ bSimpleIntegration = b; }

/**
 *  @brief Set the variable step integration step function
 *  @param integ The type of integrator to use for variable step propagations
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
 *  @brief Set the fixed step integration step function
 *  @param integ The type of integrator to use for fixed step propagations
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

//-------------------------------------------------------------------------------------------------
// 		Simulation Functions - Default to 2 Nodes
//-------------------------------------------------------------------------------------------------

/**
 *	@brief Run a simulation - initial epoch, STM, and control variables are assumed
 *
 *  @details It is assumed that t0 = 0, the STM is identity, and any control 
 *  variables are zero.
 *  
 *	@param ic an array containting the non-dimensional initial state; number
 *    of elements must match the number of <code>coreDim</code> specified 
 *    in the system dynamics model
 *	@param tof the total integration time, or time-of-flight (non-dim units)
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law to apply on all simulated segments
 */
void SimEngine::runSim(const double *ic, double tof, Arcset *arcset, ControlLaw *pLaw){
	// Define dummy value: t0 = 0; call the next level of complication
    runSim(ic, 0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation - initial epoch, STM, and control variables are assumed
 *
 *  @details It is assumed that t0 = 0, the STM is identity, and any control 
 *  variables are zero.
 *  
 *  @param ic a vector containing the IC (nondimensional states)
 *  @param tof the total integration time, or time-of-flight (non-dim units)
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law to apply on all simulated segments
 */
void SimEngine::runSim(std::vector<double> ic, double tof, Arcset *arcset, ControlLaw *pLaw){
    // Define dummy value: t0 = 0; call the next level of complication
    runSim(ic, 0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation - STM and control variables are assumed
 *  @details It is assumed that the STM is identity, and any control 
 *  variables are zero.
 * 
 *  @param ic a vector containing the IC (nondimensional states)
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law (nullptr by default)
 *  
 *  @throws Exception if `ic` has fewer than 6 elements
 *  @throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Arcset object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(std::vector<double> ic, double t0, double tof, Arcset *arcset, ControlLaw *pLaw){
    // Checks
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    unsigned int core_dim = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    if(ic.size() < core_dim){
        printErr("IC size = %zu is too small\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by coreDim in the Dynamics Model");
    }

    // Call the version that takes a pointer as the first argument; all other arguments are the same
    std::vector<double> icCopy = ic;
    runSim(&(icCopy.front()), t0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation - STM and control values are assumed
 *  @details It is assumed that the STM is identity, and any control 
 *  variables are zero.
 * 
 *  @param ic an array containing the IC (nondimensional states)
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law (nullptr by default)
 *  
 *  @throws Exception if `ic` has fewer than 6 elements
 *  @throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Arcset object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(const double* ic, double t0, double tof, Arcset *arcset, ControlLaw *pLaw){
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    // Define dummy values for Ctrl0
    std::vector<double> ctrl0;
    unsigned int ctrl_dim = 0;
    if(pLaw){
        ctrl_dim = pLaw->getNumStates();
        ctrl0.assign(ctrl_dim, 0);
    }

    // Define dummy values for STM0
    unsigned int core_dim = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    std::vector<double> stm0;
    createIdentity(stm0, core_dim + ctrl_dim);

    // Call the next level of complication
    runSim(ic, &(ctrl0.front()), &(stm0.front()), t0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simultion - STM is assumed
 *  @details It is assumed that the initial STM is identity
 * 
 *  @param ic an array containing the IC (nondimensional states)
 *  @param ctrl0 initial control state values
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law
 *  
 *  @throws Exception if `ic` has fewer than 6 elements
 *  @throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Arcset object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(std::vector<double> ic, std::vector<double> ctrl0, double t0, double tof, Arcset *arcset, ControlLaw *pLaw){
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    unsigned int core_dim = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    unsigned int ctrl_dim = pLaw ? pLaw->getNumStates() : 0;

    if(ic.size() < core_dim){
        printErr("IC size = %zu is too small\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by the Dynamics Model");
    }

    if(ctrl0.size() < ctrl_dim){
        printErr("Ctrl state size = %zu is too small\n", ctrl0.size());
        throw Exception("SimEngine::runSim: Ctrl state must have at least the number of states specifed by the Control Law");
    }

    // Define dummy values for STM0
    std::vector<double> stm0;
    createIdentity(stm0, core_dim + ctrl_dim);

    std::vector<double> icCopy = ic;
    std::vector<double> ctrlCopy = ctrl0;

    // Call the next level of complication
    runSim(&(icCopy.front()), &(ctrlCopy.front()), &(stm0.front()), t0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation - No assumptions
 *  @details No assumptions are made about the integration; all parameters are passed in
 *  or specified by SimEngine member variables
 * 
 *  @param ic a vector containing the IC (nondimensional states)
 *  @param ctrl0 initial control states
 *  @param stm0 initial STM
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law
 *  
 *  @throws Exception if `ic` has fewer than 6 elements
 *  @throws DivergeException if the GSL integrators are make steps with acceptable error values.
 *    This usually occurs if a trajectory passes very near (or through) a primary. Note that all the
 *    data generated up to the integrator failure is saved in the Arcset object passed to
 *    the SimEngine regardless of the thrown exception(s).
 */
void SimEngine::runSim(std::vector<double> ic, std::vector<double> ctrl0, const MatrixXRd &stm0, double t0, double tof, Arcset *arcset, ControlLaw *pLaw){
    // Checks
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    unsigned int core_dim = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    unsigned int ctrl_dim = pLaw ? pLaw->getNumStates() : 0;

    if(ic.size() < core_dim){
        printErr("IC size = %zu is too small\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by the Dynamics Model");
    }

    if(ctrl0.size() < ctrl_dim){
        printErr("Ctrl state size = %zu is too small\n", ctrl0.size());
        throw Exception("SimEngine::runSim: Ctrl state must have at least the number of states specifed by the Control Law");
    }

    if(stm0.rows() != core_dim + ctrl_dim || stm0.cols() != core_dim + ctrl_dim){
        printErr("STM rows = %d, cols = %d, should = %u\n", stm0.rows(), stm0.cols(), core_dim + ctrl_dim);
        throw Exception("SimEngine::runSim: Initial STM size does not match the core state size + control state size");
    }
    

    // Copy arrays to pass by reference
    std::vector<double> icCopy = ic;
    std::vector<double> ctrlCopy = ctrl0;

    // Call the version that takes a pointers; all other arguments are the same
    runSim(&(icCopy.front()), &(ctrlCopy.front()), stm0.data(), t0, tof, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation - No assumptions
 *  @details No assumptions are made about the integration; all parameters are passed in
 *  or specified by SimEngine member variables
 * 
 *  @param ic a vector containing the IC (nondimensional states)
 *  @param ctrl0 initial control states
 *  @param stm0 initial STM
 *  @param t0 the epoch time associated with the IC (non-dim units)
 *  @param tof the total integration time, or time-of-flight (non-dim units).
 *    Only the absolute value of the TOF is considered; to integrate backwards in
 *    time, use the setRevTime() function.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to a control law
 */
void SimEngine::runSim(const double* ic, const double* ctrl0, const double* stm0, double t0, double tof, Arcset *arcset, ControlLaw *pLaw){
    std::vector<double> t_span;
    // Compute the final time based on whether or not we're using reverse time integration
    double tf = bRevTime ? t0 - std::abs(tof) : t0 + std::abs(tof);
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  time will span from %.3e to %.3e\n", t0, tf);
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  (Reverse Time is %s)\n", bRevTime ? "ON" : "OFF");

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

    // Construct the tspan vector; call the next level of complication
    runSim(ic, ctrl0, stm0, t_span, arcset, pLaw);
}//=======================================================

/**
 *  @brief Run a simulation in the specified system starting with a set of initial conditions,
 *  at a specified initial time, and integrating for a specified time-of-flight
 *  
 *  @details This version of runSim is the only function that directly calls the integrate() function;
 *  all other runSim() functions reroute to this
 *  
 *  @param ic an array of non-dimensional initial states; number
 *    of elements must match the number of <code>coreDim</code> specified 
 *    in the system dynamics model
 *  @param ctrl0 initial control state vector
 *  @param stm0 initial STM data in row-major order
 *  @param t_span a vector of times to include in the solution.
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  @param pLaw pointer to the control law to apply during propagation
 */
void SimEngine::runSim(const double *ic, const double *ctrl0, const double *stm0, std::vector<double> t_span, Arcset *arcset, ControlLaw *pLaw){
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    if(!arcset->getSysData()->getDynamicsModel()->supportsControl(pLaw)){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Control law is not supported. Exiting without simulation\n");
        return;
    }

    printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, GREEN, "Running simulation...\n");
    if(!bIsClean){
        cleanEngine();
    }

    if(bMakeDefaultEvents)
        createDefaultEvents(arcset->getSysData());

    eomParams = new EOM_ParamStruct(arcset->getSysData(), pLaw);

    // Run the simulation
    bIsClean = false;   // Technically, nothing has changed yet, but this flag should be false even if any part of integrate throws an exception
    integrate(ic, ctrl0, stm0, &(t_span.front()), t_span.size(), arcset);
}//====================================================

//-------------------------------------------------------------------------------------------------
//      Simulation Functions - User Specifies Number of Nodes
//-------------------------------------------------------------------------------------------------

/**
 *  @brief Run a simulation with more than two nodes - initial epoch, STM, and control variables are assumed
 *  @details It is assumed that the initial epoch is 0, the STM is identity,
 *  and the control variables are zero.
 * 
 *  @param ic initial state vector
 *  @param tof nondimensional time-of-flight along the trajectory. The sign of <code>tof</code>
 *  is ignored; for reverse time, set the reverse time flag via setRevTime()
 *  @param numNodes Number of nodes (including the initial and final nodes) to place
 *  on the trajectory
 *  @param arcset Data structure in which to store the propagated trajectory
 *  @param pLaw control law to apply while propagating
 */
void SimEngine::runSim_manyNodes(std::vector<double> ic, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw){
    // Set t0 = 0, use more specific function
    runSim_manyNodes(ic, 0, tof, numNodes, arcset, pLaw);
}//====================================================

/**
 *  @brief Run a simulation with more than two nodes - initial epoch, STM, and control variables are assumed
 *  @details It is assumed that the initial epoch is 0, the STM is identity,
 *  and the control variables are zero.
 * 
 *  @param ic initial state vector
 *  @param tof nondimensional time-of-flight along the trajectory. The sign of <code>tof</code>
 *  is ignored; for reverse time, set the reverse time flag via setRevTime()
 *  @param numNodes Number of nodes (including the initial and final nodes) to place
 *  on the trajectory
 *  @param arcset Data structure in which to store the propagated trajectory
 *  @param pLaw control law to apply while propagating
 */
void SimEngine::runSim_manyNodes(const double *ic, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw){
    // Set t0 = 0, use more specific function
    runSim_manyNodes(ic, 0, tof, numNodes, arcset, pLaw);
}//====================================================

/**
 *  @brief Run a simulation with more than two nodes - STM and control variables are assumed
 *  @details It is assumed that the STM is identity, and the control variables are zero.
 * 
 *  @param ic initial state vector
 *  @param t0 nondimensional epoch associated with the initial state
 *  @param tof nondimensional time-of-flight along the trajectory. The sign of <code>tof</code>
 *  is ignored; for reverse time, set the reverse time flag via setRevTime()
 *  @param numNodes Number of nodes (including the initial and final nodes) to place
 *  on the trajectory
 *  @param arcset Data structure in which to store the propagated trajectory
 *  @param pLaw control law to apply while propagating
 */
void SimEngine::runSim_manyNodes(std::vector<double> ic, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw){
    // Checks
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    unsigned int core_dim = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    if(ic.size() < core_dim){
        printErr("IC size = %zu is too small\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by the Dynamics Model");
    }

    std::vector<double> tempIC = ic;
    runSim_manyNodes(&(tempIC.front()), t0, tof, numNodes, arcset, pLaw);
}//====================================================

/**
 *  @brief Run a simulation with more than two nodes - STM and control variables are assumed
 *  @details It is assumed that the STM is identity, and the control variables are zero.
 * 
 *  @param ic initial state vector
 *  @param t0 nondimensional epoch associated with the initial state
 *  @param tof nondimensional time-of-flight along the trajectory. The sign of <code>tof</code>
 *  is ignored; for reverse time, set the reverse time flag via setRevTime()
 *  @param numNodes Number of nodes (including the initial and final nodes) to place
 *  on the trajectory
 *  @param arcset Data structure in which to store the propagated trajectory
 *  @param pLaw control law to apply while propagating
 */
void SimEngine::runSim_manyNodes(const double *ic, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw){
    if(numNodes < 2){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim_manyNodes: Cannot create an arcset with fewer than two nodes\n");
        return;
    }
    int t_sign = bRevTime ? -1 : 1;

    std::vector<double> t_span(numNodes, 0);
    double t_step = std::abs(tof)/static_cast<double>(numNodes - 1);

    for(int s = 0; s < numNodes-1; s++){
        t_span[s] = t0 + t_sign*s*t_step;
    }
    t_span[numNodes-1] = t0 + t_sign*std::abs(tof);  // Final node at the desired TOF

    unsigned int coreSize = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    
    std::vector<double> ctrl0;
    unsigned int ctrlSize = 0;
    if(pLaw){
        ctrlSize = pLaw->getNumStates();
        ctrl0.assign(ctrlSize, 0);
    }

    std::vector<double> stm0;
    createIdentity(stm0, coreSize + ctrlSize);

    runSim(ic, &(ctrl0.front()), &(stm0.front()), t_span, arcset, pLaw);
}//====================================================

/**
 *  @brief Run a simulation with more than two nodes - initial STM is assumed to be identity
 *  @details [long description]
 * 
 *  @param ic initial state vector
 *  @param ctrl0 initial control state vector
 *  @param t0 nondimensional epoch associated with the initial state
 *  @param tof nondimensional time-of-flight along the trajectory. The sign of <code>tof</code>
 *  is ignored; for reverse time, set the reverse time flag via setRevTime()
 *  @param numNodes Number of nodes (including the initial and final nodes) to place
 *  on the trajectory
 *  @param arcset Data structure in which to store the propagated trajectory
 *  @param pLaw control law to apply while propagating
 */
void SimEngine::runSim_manyNodes(std::vector<double> ic, std::vector<double> ctrl0, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw){
    // Checks
    if(numNodes < 2){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::runSim_manyNodes: Cannot create an arcset with fewer than two nodes\n");
        return;
    }

    unsigned int coreSize = arcset->getSysData()->getDynamicsModel()->getCoreStateSize();
    unsigned int ctrlSize = pLaw ? pLaw->getNumStates() : 0;

    
    if(ic.size() < coreSize){
        printErr("IC size = %zu is too small\n", ic.size());
        throw Exception("SimEngine::runSim: IC must have at least the number of states specified by the Dynamics Model");
    }

    if(ctrl0.size() < ctrlSize){
        printErr("Control state size = %zu is too small\n", ctrl0.size());
        throw Exception("SimEngine::runSim: Control state must have at least the number of states specified by the Control Law");
    }

    // Space out nodes in time
    int t_sign = bRevTime ? -1 : 1;

    std::vector<double> t_span(numNodes, 0);
    double t_step = std::abs(tof)/static_cast<double>(numNodes - 1);

    for(int s = 0; s < numNodes-1; s++){
        t_span[s] = t0 + t_sign*s*t_step;
    }
    t_span[numNodes-1] = t0 + t_sign*std::abs(tof);  // Final node at the desired TOF

    // Create some dummy variables
    std::vector<double> stm0;
    createIdentity(stm0, coreSize + ctrlSize);

    std::vector<double> icCopy = ic;
    std::vector<double> ctrlCopy = ctrl0;

    runSim(&(icCopy.front()), &(ctrlCopy.front()), &(stm0.front()), t_span, arcset, pLaw);
}//====================================================

//-----------------------------------------------------
// 		Numerical Integration
//-----------------------------------------------------

/**
 *  @brief Integrate the state EOMs and STM EOMs with additional integration as required by 
 *  specific systems.
 *
 *  @details This function uses values stored in member variables to determine the direction time flows,
 *  whether or not to use simple integration, and whether or not to use variable step sizes.
 *
 *  @param ic an array containing the initial state for the trajectory; number
 *    of elements must match the number of <code>coreDim</code> specified 
 *    in the system dynamics model
 *  @param ctrl0 initial control state; number of elements must match the number
 *    of states specified in the control law.
 *  @param stm0 initial STM in row-major order; number of elements is expected to
 *    be consistent with a square matrix with side length equal to the number of
 *    <code>coreDim</code>.
 *  @param t an array of times to integrate over; may contain 2 elements (t0, tf), or a range of times
 *  @param t_dim the dimension of t
 *  @param arcset pointer to a trajectory object to store the output trajectory
 *  
 *  @throws DivergeException if the integrator fails and cannot proceed
 */
void SimEngine::integrate(const double *ic, const double *ctrl0, const double *stm0, const double *t, unsigned int t_dim, Arcset *arcset){
    if(arcset == nullptr){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "SimEngine::integrate: Arcset pointer is nullptr; exiting to avoid memory leaks\n");
        return;
    }

    // delete any auto-generated events from previous integrations
    for(unsigned int i = 0; i < events.size(); i++){
        if(static_cast<int>(events[i].getType()) <= 0){
            events.erase(events.begin() + i);
            i--;
        }
    }

    time(&startTimestamp);  // Get the starting time

    // Save tolerance for trajectory
    arcset->setTol(absTol > relTol ? absTol : relTol);
    const DynamicsModel *model = arcset->getSysData()->getDynamicsModel();

    // Initialize all events with the correct system data pointer
    for(unsigned int i = 0; i < events.size(); i++){
        events[i].initialize(arcset->getSysData());
    }

    /**
     *  Data is stored in the state vector according to the following structure
     *  
     *  q = [core_state; control_state; stm_elements; extra_state]
     *  
     *      core_state      -   a (core_dim x 1) vector that contains the "core state," e.g.,
     *                          the position, velocity, mass of the spacecraft
     *      control_state   -   a (ctrl_dim x 1) vector that contains control state information
     *      stm_elements    -   represents a (core_dim + ctrl_dim x core_dim + ctrl_dim) state 
     *                          transition matrix in row-major order
     *      extra_state     -   a (extra_dim x 1) vector that contains "extra states," e.g.,
     *  
     *  If simpleIntegration is enabled, the STM and extra states are not included in the
     *  integration and their values in the Segment state array are filled by zeros 
     */
    const unsigned int core_dim = model->getCoreStateSize();            // The number of states in the most basic EOM propagation
    const unsigned int ctrl_dim = eomParams->pCtrlLaw ? eomParams->pCtrlLaw->getNumStates() : 0;  // The number of independent control variables
    const unsigned int extra_dim = model->getExtraStateSize();          // The number of extra states
    

    const unsigned int stm_dim = pow(core_dim + ctrl_dim, 2);           // Size of the state-transition matrix
    const unsigned int ic_dim = core_dim + ctrl_dim + (!bSimpleIntegration)*(stm_dim + extra_dim);     // Number of states used for this propagation
    const unsigned int full_dim = core_dim + ctrl_dim + stm_dim + extra_dim;       // Max number of states for the most complex EOM propagation for this model
    
    // dummy extra state variables to save if not all the states are propagated
    const std::vector<double> extraStates(full_dim - ic_dim, 0);

    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  IC has %u initial states\n", ic_dim);

    // Construct the full IC from the state ICs plus the STM ICs and any other ICs for more complex systems
    std::vector<double> fullIC(ic_dim, 0);
    std::copy(ic, ic + core_dim, &(fullIC.front()));

    // Add the control states to the vector
    std::copy(ctrl0, ctrl0+ctrl_dim, &(fullIC[core_dim]));

    // NOTE: STM follows immediately after core and control states; any extras come after the STM
    if(!bSimpleIntegration){
        std::copy(stm0, stm0 + stm_dim, &(fullIC[core_dim + ctrl_dim]));
    }        

    double *y = &(fullIC.front());      // array of states that is passed to the integrator

    // Choose EOM function based on system type and simplicity
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  using %s integration\n", bSimpleIntegration ? "simple (no STM)" : "full (+ STM)");
    int (*eomFcn)(double, const double[], double[], void*) = 
        bSimpleIntegration ? model->getSimpleEOM_fcn() : model->getFullEOM_fcn();     // Pointer for the EOM function

    printVerb(verbosity >= Verbosity_tp::ALL_MSG, 
        "  using control law: %s\n", eomParams->pCtrlLaw ? eomParams->pCtrlLaw->getTypeString().c_str() : "NONE");

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
    //     boost_observer boostObserver(model, arcset, eomParams);

    //     double t0 = t[0], tf = t[1];                // start and finish times for integration; t0 will be updated by integrator
    //     double dt = tf > t0 ? dtGuess : -dtGuess;   // step size (initial guess)

    //     // Save the initial state, time, and STM
    //     // model->sim_saveIntegratedData(y, t[0], arcset, eomParams);

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

    /* Create a system to integrate; we don't include a Jacobin (nullptr)
     *  The parameter set eomParams can be modified 
     *  between integration steps (i.e., change model parameters), but the ode functions must be reset
     *  via <code>gsl_odeiv2_driver_reset</code>, <code>gsl_odeiv2_evolve_reset</code>, or
     *  <code>gsl_odeiv2_step_reset</code> before continuing with an updated parameter set
     */
    gsl_odeiv2_system odeSys = {eomFcn, nullptr, static_cast<size_t>(ic_dim), eomParams};
    
    // Define ODE objects, define them conditionaly based on bVarStepSize
    gsl_odeiv2_step *s = nullptr;
    gsl_odeiv2_control *c = nullptr;
    gsl_odeiv2_evolve *e = nullptr;
    gsl_odeiv2_driver *d = nullptr;

    if(bVarStepSize){
        // Allocate space for the stepping object
        switch(varStep_integ){
            case Integ_tp::RK8PD:
                s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, ic_dim);
                printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  variable step size, using Runge-Kutta Dormand-Prince 8-9 method\n");
                break;
            case Integ_tp::RKCK:
                s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, ic_dim);
                printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  variable step size, using Runge-Kutta Cash-Karp 4-5 method\n");
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
        printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  fixed step size, using Adams-Bashforth, Adams-Moulton method\n");
        
        double signed_dt = bRevTime ? -1*dtGuess : dtGuess;

        switch(fixStep_integ){
            case Integ_tp::MSADAMS:
                // Allocate space for a driver; the msadams algorithm requires access to the driver
                d = gsl_odeiv2_driver_alloc_y_new(&odeSys, gsl_odeiv2_step_msadams, signed_dt, absTol, relTol);
                gsl_odeiv2_driver_set_nmax(d, maxDriverSteps);
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
    printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  sim will use %d event functions:\n", static_cast<int>(events.size()));
    for(unsigned int ev = 0; ev < events.size(); ev++){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  [%02u] - %s\n", ev, events.at(ev).getTypeStr());
        events.at(ev).updateDist(y, core_dim, t[0]);
    }

    // Save the initial state, time, and STM
    Node node0(y, core_dim, t[0]); // Construct a basic node
    int id0 = model->sim_addNode(node0, y, t[0], arcset, eomParams, Event_tp::NONE); // Pass the node to the Dynamic Model so specific modifications may be made
    
    // Construct a segment with origin at the initial node, undefined terminus, dummy TOF in the correct direction
    Segment seg(id0, Linkable::INVALID_ID, t[1] - t[0]);
    seg.appendState(y, ic_dim);     // Save the initial state in the segment
    seg.setStateWidth(ic_dim);
    seg.appendState(extraStates);   // Save dummy values for extra states that are not propagated
    seg.appendTime(t[0]);           // Save the initial time in the segment
    model->sim_addSeg(seg, y, t[0], arcset, eomParams);

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
                // printColor(GREEN, "Driver took %d steps for %.2f nd time\n", d->n, t_int);
            }

            if(status != GSL_SUCCESS){
                reportPropErrs(status, t_int);
                free_odeiv2(s, c, e, d);

                // Save the last successful step the integrator was able to take
                Node nodeF(y, core_dim, t_int);                 // Construct a basic node
                Segment &lastSeg = arcset->getSegRefByIx(-1);   // Get reference to most recent segment (faster than copying the object!)

                int idf = model->sim_addNode(nodeF, y, t_int, arcset, eomParams, Event_tp::SIM_ERR); // Pass the node to the DynamicModel for more specific actions
                arcset->getNodeRef(idf).addLink(lastSeg.getID());

                lastSeg.setTerminus(idf);                       // Update the terminus to the final node
                lastSeg.appendState(y, ic_dim);                 // Save the final state and time to the segment
                lastSeg.appendState(extraStates);               // Save dummy values for extra states that are not propagated
                lastSeg.appendTime(t_int);
                lastSeg.updateTOF();                             // Compute the actual TOF from the data and store in dedicated TOF variable
                lastSeg.setSTM(y+core_dim+ctrl_dim, stm_dim);   // Set the STM to be the most recent one

                throw DivergeException("SimEngine::integrate: Integration did not succeed");
            }

            // Stop the simulation if a simulation-ending event occurs
            killSim = locateEvents(y, t_int, arcset, propStepCount);

            // Stop the simulation if the maximum computation time has passed
            killSim = killSim || (maxCompTime > 0 && (time(nullptr) - startTimestamp) > maxCompTime);

            if(killSim)
                break;

            // Save propagation state to segment
            Segment &lastSeg = arcset->getSegRefByIx(-1);       // Get another reference (not reassignable)
            lastSeg.appendState(y, ic_dim);                     // Save newest state and time
            lastSeg.appendState(extraStates);                   // Save dummy values for extra states that are not propagated
            lastSeg.appendTime(t_int);
            propStepCount++;
        }
    }else{
        // Integrate each segment between the input times
        for (unsigned int j = 0; j < t_dim - 1; j++){
            // define start and end times; t_int will be updated by integrator
            t_int = t[j];
            double tf = t[j+1];
            double dt = tf > t_int ? dtGuess : -dtGuess;
            int sgn = tf > t_int ? 1 : -1;

            if(j > 0){
                // Create a node at every time in the t[] array
                Node nodeI(y, core_dim, t_int);
                Segment &lastSeg = arcset->getSegRefByIx(-1);     // Get reference (not reassignable)

                int nodeID_i = model->sim_addNode(nodeI, y, t_int, arcset, eomParams, Event_tp::SIM_TOF);
                arcset->getNodeRef(nodeID_i).addLink(lastSeg.getID());

                lastSeg.setTerminus(nodeID_i);  // Update the previous segment to terminate at the new node
                lastSeg.updateTOF();
                lastSeg.setSTM(y+core_dim+ctrl_dim, stm_dim);

                // Create a new segment for the next time interval - origin is correct, terminus is undetermined, tof is approximate
                Segment newSeg(nodeID_i, Linkable::INVALID_ID, tf - t_int);
                newSeg.appendState(y, ic_dim);
                newSeg.setStateWidth(ic_dim);
                lastSeg.appendState(extraStates);   // Save dummy values for extra states that are not propagated
                newSeg.appendTime(t_int);
                model->sim_addSeg(newSeg, y, t_int, arcset, eomParams);
            }
            while(sgn*t_int < sgn*tf && !killSim){
                // printf("Integrating at t = %6.4f\n", t_int);
                if(bVarStepSize){
                    status = gsl_odeiv2_evolve_apply(e, c, s, &odeSys, &t_int, tf, &dt, y);
                }else{
                    status = gsl_odeiv2_driver_apply(d, &t_int, tf, y);
                }

                if(status != GSL_SUCCESS){
                    reportPropErrs(status, t_int);
                    free_odeiv2(s, c, e, d);

                    // Save the last successful step the integrator was able to take
                    Node nodeF(y, core_dim, t_int);    // Construct a basic node
                    Segment &lastSeg = arcset->getSegRefByIx(-1);     // Get another reference (not reassignable)

                    // Add the node with info that there was a simulation error here
                    int idf = model->sim_addNode(nodeF, y, t_int, arcset, eomParams, Event_tp::SIM_ERR);
                    arcset->getNodeRef(idf).addLink(lastSeg.getID());

                    lastSeg.setTerminus(idf);                       // Update the terminus to the final node
                    lastSeg.appendState(y, ic_dim);                 // Save the final state and time to the segment
                    lastSeg.appendState(extraStates);               // Save dummy values for extra states that are not propagated
                    lastSeg.appendTime(t_int);
                    lastSeg.updateTOF();
                    lastSeg.setSTM(y+core_dim+ctrl_dim, stm_dim);

                    throw DivergeException("SimEngine::integrate: Integration did not succeed");
                }

                killSim = locateEvents(y, t_int, arcset, propStepCount);

                // Stop the simulation if the maximum computation time has passed
                killSim = killSim || (maxCompTime > 0 && difftime(time(nullptr), startTimestamp) > maxCompTime);

                propStepCount++;

                // Save the most recent time and state to the segment
                Segment &lastSeg = arcset->getSegRefByIx(-1);   // Get another reference (not reassignable)
                lastSeg.appendState(y, ic_dim);
                lastSeg.appendState(extraStates);               // Save dummy values for extra states that are not propagated
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
        Segment &lastSeg = arcset->getSegRefByIx(-1);

        // Add the node and add a link to the previous segment
        int nodeID_f = model->sim_addNode(nodeF, y, t_int, arcset, eomParams, Event_tp::SIM_TOF);
        arcset->getNodeRef(nodeID_f).addLink(lastSeg.getID());

        lastSeg.setTerminus(nodeID_f);  // Just update the terminus: final state and time should have already been appended
        lastSeg.updateTOF();
        lastSeg.setSTM(y+core_dim+ctrl_dim, stm_dim);
    }else if(maxCompTime > 0 && difftime(time(nullptr), startTimestamp) > maxCompTime){
        // Ended at the time-out
        // Create a final node and update the final segment
        Node nodeF(y, core_dim, t_int);
        Segment &lastSeg = arcset->getSegRefByIx(-1);

        // Add the node and add a link to the previous segment
        int nodeID_f = model->sim_addNode(nodeF, y, t_int, arcset, eomParams, Event_tp::SIM_COMPTIME);
        arcset->getNodeRef(nodeID_f).addLink(lastSeg.getID());
    
        // Propagation ended without saving the final state
        lastSeg.setTerminus(nodeID_f);
        lastSeg.appendState(y, ic_dim);
        lastSeg.appendState(extraStates);   // Save dummy values for extra states that are not propagated
        lastSeg.appendTime(t_int);
        lastSeg.updateTOF();
        lastSeg.setSTM(y+core_dim+ctrl_dim, stm_dim);
    }

    // Check lengths of vectors and set the numPoints value in arcset
    printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Integration complete**\n");
    if(verbosity >= Verbosity_tp::ALL_MSG){
        arcset->print();
        arcset->printInChrono();
    }

    for(unsigned int n = 0; n < arcset->getNumNodes(); n++){
        printVerb(verbosity >= Verbosity_tp::SOME_MSG, " %s Event occured at node %u\n",
            Event::getEventTpStr(arcset->getNodeRefByIx(n).getTriggerEvent()), n);
    }
}//====================================================END of cr3bp_integrate

/**
 *  @brief Free all GSL ODE pointers that have been initialized/instantiated
 * 
 *  @param s stepper object pointer
 *  @param c control object pointer
 *  @param e evolver object pointer
 *  @param d driver object pointer
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
 *  @brief Locate event occurences as exactly as possible and determine if the simulation 
 *  should end because of the event.
 *
 *  Look through the list of events and check each one to see if it has occured. This 
 *  initial check is not highly accurate, but if an event is determined to have occured,
 *  a Newton-Raphson process is begun to locate the exact state and time of the event.
 *
 *  The initial check compares the current integrated state (passed in as `y`) with
 *  the previous integrated state (stored in the event object). If the sign (+/-) of the 
 *  distance to the event function changes, the the trajectory has triggered the event. This
 *  check is performed in the event objects `crossedEvent()` function.
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
 *  @param pArcset pointer to a trajectory object to store the output trajectory
 *  @param propStepCount the number of steps the propagation has taken so far. This is different
 *  from the number of time steps as many propagation steps may be taken between specified time steps
 *  @return whether or not the simulation should end (an event triggers killSim)
 */
bool SimEngine::locateEvents(const double *y, double t, Arcset *pArcset, int propStepCount){
    const DynamicsModel *model = pArcset->getSysData()->getDynamicsModel();
    unsigned int core_dim = model->getCoreStateSize();

    // Look through all events
    bool continueSim = true;
    for(unsigned int ev = 0; ev < events.size(); ev++){
        // Don't trigger if only two points have been integrated
        // Also don't trigger if the type is 0 or negative: these are managed by the simulation engine
        if(propStepCount > 1 && to_underlying(events[ev].getType()) > 0 && 
            events[ev].crossedEvent(y, core_dim, t, bRevTime ? -1 : 1)){

            printVerb(verbosity >= Verbosity_tp::ALL_MSG,
                "  Event %d detected on segment %d; searching for exact crossing\n", ev, pArcset->getNumSegs());
            events[ev].incrementCount();  // Update the counter for the event

            if(verbosity >= Verbosity_tp::ALL_MSG){ events[ev].printStatus(); }

            // Use correction to locate the event very accurately
            if(locateEvent_multShoot(y, t, ev, pArcset)){
                // This condition is also checked in model->sim_locateEvent(): that function adds a new
                // segment if the simulation is continuing
                if(events[ev].stopOnEvent() && events[ev].getTriggerCount() >= events[ev].getStopCount()){
                    printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Completed Event Location, ending integration**\n");
                    // No need to remember the most recent point; it will be discarded, leaving
                    // the point from mult. shooting as the last
                    continueSim = continueSim && false;    // Tell the simulation to stop
                }else{
                    printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, GREEN, "**Completed Event Location, continuing integration**\n");
                    // Save the most recent point (the one after the event that triggered the direction change in events[ev])
                    events[ev].updateDist(y, core_dim, t);
                    continueSim = continueSim && true;  // Simulation will continue unless another event has occurred
                }
            }
        }else{
            if(to_underlying(events[ev].getType()) > 0){
                // Save the distance and current state to the event
                events[ev].updateDist(y, core_dim, t);
            }
        }
    }// end of loop

    return !continueSim;
}//========================================

/**
 *  @brief Use a multiple shooting correction algorithm to accurately locate an event crossing
 *
 *  The simulation engine calls this function if and when it determines that an event 
 *  has been crossed. To accurately locate the event, we employ differential corrections
 *  and find the exact event occurence in space and time.
 *
 *  @param y full state array at the current integration step
 *  @param t time at the current integration step
 *  @param evtIx index of the event being located
 *  @param pArcset a pointer to the arcset the event occurred on
 *
 *  @return wether or not the event has been located. If it has, a new Node
 *  has been appended to the arcset. If propagation will continue, a new
 *  segment has also been created
 */
bool SimEngine::locateEvent_multShoot(const double *y, double t, int evtIx, Arcset* pArcset){

    const DynamicsModel *model = pArcset->getSysData()->getDynamicsModel();
    const unsigned int core_dim = model->getCoreStateSize();
    const unsigned int ctrl_dim = eomParams->pCtrlLaw ? eomParams->pCtrlLaw->getNumStates() : 0;

    // Create a nodeset from the previous state (stored in the event) and
    // integrate forwards for half the time between this state and the last one
    Segment &lastSeg = pArcset->getSegRefByIx(-1);

    double t0 = 0, tof = 0, ti = NAN;
    std::vector<double> arcIC(lastSeg.getStateWidth()), arcFC(lastSeg.getStateWidth());
    if(lastSeg.getNumTimes() == 1){
        t0 = lastSeg.getTimeByIx(0);
        tof = 0.5*(t - t0);
        arcIC = lastSeg.getStateByRow(0);
        arcFC = std::vector<double>(y, y + arcIC.size());
    }else{
        t0 = lastSeg.getTimeByIx(-2);           // Time from the state before last
        double ti = lastSeg.getTimeByIx(-1);    // Time from the previous state
        tof = t - t0 - 0.5*(t - ti);            // Approx. TOF 
        
        // Copy IC into vector - Use the state from two iterations ago to avoid
        // numerical problems when the previous state is REALLY close to the event
        arcIC = lastSeg.getStateByRow(-2);
        arcFC = lastSeg.getStateByRow(-1);
    }

    if(verbosity >= Verbosity_tp::ALL_MSG){
        // printColor(BLUE, "Step index = %d\n", propStepCount-1);
        printColor(BLUE, "t(now) = %f\nt(prev) = %f\nt(prev-1) = %f\n", t, ti, t0);
        printColor(BLUE, "State(now) = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", y[0],
            y[1], y[2], y[3], y[4], y[5]);
        printColor(BLUE, "tof = %f\n", tof);
        printColor(BLUE, "ic = [%9.4e %9.4e %9.4e %9.4e %9.4e %9.4e]\n", arcIC[0],
            arcIC[1], arcIC[2], arcIC[3], arcIC[4], arcIC[5]);
    }   

    // **************************************************
    // Use correction to locate the event very accurately
    // **************************************************
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  Creating arcset for event location\n");
    
    Arcset eventArcset(pArcset->getSysData()), correctedSet(pArcset->getSysData());

    // Construct a simple arcset from previous states we know work well
    Node n0(&(arcIC.front()), core_dim, t0);
    model->sim_addNode(n0, &(arcIC.front()), t0, &eventArcset, eomParams, Event_tp::NONE);

    // Determine whether or not we're using an "ENDSEG" type constraint
    bool useEndSeg;
    switch(events[evtIx].getConType()){
        case Constraint_tp::ENDSEG_STATE:
        case Constraint_tp::ENDSEG_JC:
        case Constraint_tp::ENDSEG_APSE:
        case Constraint_tp::ENDSEG_DIST:
        case Constraint_tp::ENDSEG_MIN_DIST:
        case Constraint_tp::ENDSEG_MAX_DIST:
            useEndSeg = true;
            break;
        default:
            useEndSeg = false;
            break;
    }

    Segment tempSeg(0, Linkable::INVALID_ID, tof);
    tempSeg.appendState(arcIC);
    tempSeg.appendState(arcFC);
    tempSeg.setStateWidth(arcIC.size());

    if(useEndSeg){
        // No need to add a terminal node, just add the segment
        model->sim_addSeg(tempSeg, &(arcFC[0]), t0, &eventArcset, eomParams);

        // Constrain the end of the segment
        Constraint eventCon(events[evtIx].getConType(), 0, events[evtIx].getConData());
        eventArcset.addConstraint(eventCon);
    }else{
        // Add a final node
        Node nf(&(arcFC.front()), core_dim, t0+tof);
        model->sim_addNode(nf, &(arcFC.front()), t0+tof, &eventArcset, eomParams, Event_tp::SIM_TOF);

        // Add a segment to link the two
        tempSeg.setTerminus(1);
        model->sim_addSeg(tempSeg, &(arcFC[0]), t0, &eventArcset, eomParams);

        // Constrain the final node
        Constraint eventCon(events[evtIx].getConType(), 1, events[evtIx].getConData());
        eventArcset.addConstraint(eventCon);
    }

    // Remove the initial state from the design vector; it is fixed.
    // This reduces the size of the Jacobian matrix and expedites the shooting process
    Constraint rmInitState(Constraint_tp::RM_STATE, 0, nullptr, 0);
    eventArcset.addConstraint(rmInitState);

    // Remove the initial epoch too, if applicable
    if(model->supportsCon(Constraint_tp::EPOCH)){
        Constraint rmInitEpoch(Constraint_tp::RM_EPOCH, 0, nullptr, 0);
        eventArcset.addConstraint(rmInitEpoch);
    }

    if(verbosity >= Verbosity_tp::ALL_MSG){ eventArcset.print(); }

    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  Applying corrections process to locate event\n");
    MultShootEngine corrector;
    corrector.setTOFType(MSTOF_tp::VAR_FREE);
    corrector.setTol(pArcset->getTol());
    corrector.setVerbosity(verbosity);
    corrector.setFindEvent(true);   // apply special settings to minimize computations
    corrector.setFullFinalProp(false);
    
    // Because we set findEvent to true, this output nodeset should contain
    // the full (42 or 48 element) final state
    try{
        corrector.multShoot(&eventArcset, &correctedSet);
    }catch(DivergeException &e){
        if(verbosity >= Verbosity_tp::SOME_MSG)
            printErr("Unable to locate %s event; corrector diverged\n", events[evtIx].getTypeStr());
        return false;
    }catch(LinAlgException &e){
        if(verbosity >= Verbosity_tp::SOME_MSG)
            printErr("LinAlg Err while locating %s event; bug in corrector!\n", events[evtIx].getTypeStr());
        return false;
    }

    double eventTime = t0 + correctedSet.getTotalTOF();
    // Get the full (including STM, etc) state from the end of the propagated segment
    std::vector<double> state = correctedSet.getSegRefByIx(-1).getStateByRow(-1);
    // Update event state from the most recent node (a new node was created at the event occurence)
    events[evtIx].updateDist(&(state.front()), core_dim, eventTime);

    // TODO - eventually all events should use this ENDSEG type behavior
    Node evtNode;
    int id;
    if(useEndSeg){
        // The final node was deleted; create one to represent the event
        evtNode = Node(&(state.front()), core_dim, eventTime);
        id = model->sim_addNode(evtNode, &(state.front()), eventTime, pArcset, eomParams, events[evtIx].getType());
    }else{
        // The final node corresponds with the event
        evtNode = correctedSet.getNodeByIx(-1);
        evtNode.clearConstraints();     // Get rid of any cosntraints used to isolate the event
        // Add the node; pass nullptr for the state vector so that none of the state data in evtNode is overwritten (it should be complete)
        id = model->sim_addNode(evtNode, nullptr, eventTime, pArcset, eomParams, events[evtIx].getType());
    }

    pArcset->getNodeRef(id).addLink(lastSeg.getID());    // Link the new node to the previous segment

    lastSeg.setTerminus(id);
    lastSeg.appendState(state);
    lastSeg.appendTime(eventTime);
    lastSeg.updateTOF();
    lastSeg.setSTM(&(state[core_dim+ctrl_dim]), (core_dim+ctrl_dim)*(core_dim+ctrl_dim));

    // Create a new segment if the propagation is going to continue
    if(!(events[evtIx].stopOnEvent() && events[evtIx].getTriggerCount() >= events[evtIx].getStopCount())){
        // Initialize with origin ID, undetermined terminus, and a dummy TOF with the same sign as the previous segment
        Segment newSeg(id, Linkable::INVALID_ID, lastSeg.getTOF());
        newSeg.appendState(state);
        newSeg.setStateWidth(state.size());
        newSeg.appendTime(eventTime);
        model->sim_addSeg(newSeg, &(state.front()), eventTime, pArcset, eomParams);
    }

    return true;    // If the code reached this point, the event was definitely located
}//=======================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  @brief Clean out the trajectory storage variable so a new simulation can be run and store its data
 *  @details This function does not reset any parameters the user has set
 */
void SimEngine::cleanEngine(){
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Cleaning the engine...\n");
    
    if(eomParams)
        delete eomParams;
    eomParams = nullptr;  // set pointer to 0 (nullptr pointer)

    for(unsigned int e = 0; e < events.size(); e++){
        events[e].reset();
    }

    bIsClean = true;
}//====================================================

/**
 *  Create default events for the system
 *  @param sysData pointer to the system data object for the simulation
 */
void SimEngine::createDefaultEvents(const SysData *sysData){
    if(!bMadeDefaultEvents){
        
        std::vector<Event> defEvents = sysData->getDynamicsModel()->sim_makeDefaultEvents(sysData);
        for(unsigned int i = 0; i < defEvents.size(); i++)
            addEvent(defEvents[i]);

        bMadeDefaultEvents = true;
    }
}//====================================================

void SimEngine::reportPropErrs(int gsl_status, double t_int){
    printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED,
        "SimEngine::integrate: t = %.4e, GSL ERR: %s\n", t_int, gsl_strerror(gsl_status));

    switch(gsl_status){
        case GSL_EMAXITER:
            printVerbColor(verbosity >= Verbosity_tp::NO_MSG, RED,
                "  Details: SimEngine::integrate: maximum number of driver steps have occurred\n");
            break;
        case GSL_ENOPROG:
            printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED,
                "  Details: SimEngine::integrate: driver step size has dropped below the minimum value; no progress to be made\n");
            break;
        case GSL_EBADFUNC:
            printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED,
                "  Details: SimEngine::integrate: 'bad function' error thrown\n");
            break;  // Driver must be reset before using again, but current behavior ends the integration and frees driver, so not necessary
        default: 
            break; // Nothing to do!
    }
}//====================================================
/**
 *  @brief Reset all variables and options
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
}//====================================================

/**
 *  Clear all events from the simulation, including any created by default.
 */
void SimEngine::clearEvents(){
    printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Clearing all events...\n");
    events.clear();
    bMadeDefaultEvents = false;
}//====================================================

/**
 *  Copy data from an input engine to this one
 *  @param s an input simulation engine
 *  @throw Exception if `s` has an unknown system data type
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
}//====================================================

}// END of Astrohelion namespace