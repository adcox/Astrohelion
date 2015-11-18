/**
 * Run a series of benchmark tests for computational time
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_event.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <cstdlib>
#include <time.h>
#include <sys/time.h>

double getCPUTime(){ return (double)clock() / CLOCKS_PER_SEC; }

void runVanillaSims();
void runEventSims();
void runCorrector();

int main(){
	double t0 = 0, tf = 0;

	printf("Running Benchmark Tests...\n");

	printf("  runVanillaSims()...\n");
	t0 = getCPUTime();
	runVanillaSims();
	tf = getCPUTime();
	double dt_runVanillaSims = tf - t0;

	printf("  runEventSims()...\n");
	t0 = getCPUTime();
	runEventSims();
	tf = getCPUTime();
	double dt_runEventSims = tf - t0;

	printf("  runCorrector()...\n");
	t0 = getCPUTime();
	runCorrector();
	tf = getCPUTime();
	double dt_runCorrector = tf - t0;

	printf("\n**************************\n");
	printf(" Benchmark Results\n");
	printf("\n**************************\n");
	printf("runVanillaSims() -> %.12f seconds\n", dt_runVanillaSims);
	printf("runEventSims() -> %.12f seconds\n", dt_runEventSims);
	printf("runCorrector() -> %.12f seconds\n", dt_runCorrector);

	return EXIT_SUCCESS;
}//=============================================================


/**
 *  @brief Run a bunch of simulations with no extra stopping conditions
 */
void runVanillaSims(){

	tpat_sys_data_cr3bp sys("earth", "moon");

	// EM 3:5 Resonant Orbit
	double IC[] = {0.604740830000000, 0, 0, 0, 1.007264110400000, 0};
	double tof = 30.97096176;

	tpat_simulation_engine sim(&sys);
	
	for(int n = 0; n < 100; n++){
		sim.runSim(IC, tof);	
	}
}//=============================================================

/**
 *  @brief Run a bunch of simulations with events that locate
 *  XZ and YZ plane crossings; these events do not stop the integration
 *  so the results should be comparable to runVanillaSims()
 */
void runEventSims(){
	tpat_sys_data_cr3bp sys("earth", "moon");

	// EM 3:5 Resonant Orbit
	double IC[] = {0.604740830000000, 0, 0, 0, 1.007264110400000, 0};
	double tof = 30.97096176;

	tpat_simulation_engine sim(&sys);
	sim.addEvent(tpat_event::YZ_PLANE, 0, false);
	sim.addEvent(tpat_event::XZ_PLANE, 0, false);
	for(int n = 0; n < 100; n++){
		sim.runSim(IC, tof);	
	}
}//=============================================================

void runCorrector(){
	tpat_sys_data_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	// Create a nodeset
	tpat_nodeset_cr3bp nodeset(IC, &sys, tof, 15);
	tpat_traj_cr3bp traj = tpat_traj_cr3bp::fromNodeset(nodeset);

	// Constraint node 7 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	tpat_correction_engine corrector;
	
	for(int n = 0; n < 100; n++){
		tpat_nodeset_cr3bp temp = nodeset;
		corrector.correct(&nodeset);
	}
}//=============================================================