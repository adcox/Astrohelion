/**
 * 	Simulation Engine
 *
 *	This object handles all simulation tasks.
 *
 */

#include "adtk_calculations.hpp"
#include "adtk_matrix.hpp"
#include "adtk_simulation_engine.hpp"
#include "adtk_cr3bp_traj.hpp"

// #include <fstream>
// #include <cstdio>

#include <cstdlib>
#include <iostream>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_simulation_engine::adtk_simulation_engine(){
	revTime = false;
	verbose = false;
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
        cout << "This trajectory was not propagated in CR3BP";
        throw;
    }
}//==============================================

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
 *	@param tf the total integration time
 */
void adtk_simulation_engine::runSim(double ic[], double tf){
	this->runSim(ic, 0, tf);
}//=======================================================

/**
 *	
 *	@param ic
 *	@param t0
 *	@param tf
 */
void adtk_simulation_engine::runSim(double ic[], double t0, double tf){
	// Instantiate empty vector for state data
	vector<double> simData;

	// Create time span using revTime to adjust limits
	double t_span[2] = {t0, 0};
	t_span[1] = revTime ? t0 - tf : t0 + tf;

	switch(sysData->getType()){
		case adtk_sys_data::CR3BP_SYS:
			//do cr3bp simulation
			simData = cr3bp_integrate(ic, t_span, sysData->getMu(), 2);
			break;
		case adtk_sys_data::BCR4BPR_SYS:
			//do bcr4bp simulation
		default:
			cout << "Cannnot simulate with system type: " << sysData->getTypeStr() << endl;
			return;
	}

	// Compute number of integration steps, break data into different vectors for state and time
	int steps = (int)(simData.size()/42);
	if(verbose)
		cout << "Integration took " << steps << " steps"<<endl;

    switch(sysData->getType()){
        case adtk_sys_data::CR3BP_SYS:
            traj = new adtk_cr3bp_traj(steps);
            break;
        case adtk_sys_data::BCR4BPR_SYS:
            //do stuff for BCR4BP
            break;
        default:
            break;
    }

    // Get pointers to the state, time, and STM elements
    // state vector contains 9 elements in each row: pos, vel, accel
    vector<double>* state = traj->getState();
    vector<double>* t = traj->getTime();
    vector<adtk_matrix>* allSTM = traj->getSTM();

	// For now, save data to file
	// ofstream dataFile;
	// dataFile.open("data.csv", ios::out);

	for(int n = 0; n < steps; n++){

		state->at(n*9) = simData[n*43];
		state->at(n*9 + 1) = simData[n*43 + 1];
		state->at(n*9 + 2) = simData[n*43 + 2];
		state->at(n*9 + 3) = simData[n*43 + 3];
		state->at(n*9 + 4) = simData[n*43 + 4];
		state->at(n*9 + 5) = simData[n*43 + 5];

        t->at(n) = simData[(n+1)*43 - 1];   // Time is appended to the end of the 42 states

		// Compute acceleration
		double s[6] = {0}, dsdt[6] = {0};
		copy(simData.begin()+n*43, simData.begin()+n*43+6, s);

		switch(sysData->getType()){
			case adtk_sys_data::CR3BP_SYS:
				// Use the simple EOMs to compute the velocity/acceleration. 
				// Note that the time t is not used in computations, and mu may need to be a pointer...
				cr3bp_simple_EOMs(s, dsdt, sysData->getMu());

                // TODO Compute Jacobi Constant - put that function in adtk_calculations
				break;
			case adtk_sys_data::BCR4BPR_SYS:
				// Compute acceleration for BCR4BP
				break;
			default:
				cout << "Unknown sim type " << sysData->getTypeStr() << endl;
		}

		// Copy acceleration into state vector
		state->at(n*9 + 6) = dsdt[3];
		state->at(n*9 + 7) = dsdt[4];
		state->at(n*9 + 8) = dsdt[5];

		// char dataStr[256];
		// sprintf(dataStr, "%.20f %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f",
		// 	t[n], state[n][0], state[n][1], state[n][2], state[n][3], state[n][4], state[n][5], state[n][6], 
		// 	state[n][7], state[n][8]);
		// dataFile << dataStr << endl;

		// Copy STM elements into a sub-array, turn it into a matrix, insert into vector of STMs
		double stmElm[36];
		copy(simData.begin()+n*43+7, simData.begin()+(n+1)*43-2, stmElm);
		adtk_matrix STM(6,6, stmElm);
        allSTM->at(n) = STM;

	}

	// dataFile.close();
}

//-----------------------------------------------------
// 		Numerical Integration
//-----------------------------------------------------

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
vector<double> adtk_simulation_engine::cr3bp_integrate(double ic[], double t[], double mu, int t_dim){

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

    int steps = 0;                          // count number of integration steps
    vector<double> state;                   // hold the entire integrated state
    double *y = new double[ic_dim];         // hold ONE integrated state; resets every step 

    // Define the step type (or, which integrator are we using?)
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_msadams;

	// Create a system to integrate; we don't include a Jacobian (NULL)
    gsl_odeiv2_system sys = {cr3bp_EOMs, NULL, ic_dim, &mu};
    
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

    // put initial state vector [ic, t(0)] into the state vector
    for(int i = 0; i < ic_dim+1; i++){
        state.push_back(i < ic_dim ? y[i] : t[0]);
    }
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
    gsl_odeiv2_driver_free(d);
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    delete(y);

    cout << "Size of state: " << state.size() << endl;
    state.shrink_to_fit();
    cout << "Size of state: " << state.size() << endl;
    
    return state;   // return the array of states
}//END of cr3bp_integrate