/**
 *	Compute the periapsis and stability indices for the specified family
 */

#include "tpat.hpp"

#include "tpat_body_data.hpp"
#include "tpat_calculations.hpp"
#include "tpat_event.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include "matio.h"
#include <vector>

/** Where is the family data file located? **/
const char *inputFam = "../share/families_pac_checked/EM_L1_Halo.mat";

/** Where to save the output data? **/
const char *outputData = "data/EML1_Halo_PeriStability.mat";

int main(){
	// Load the family
	tpat_family_cr3bp halos(inputFam);

	// Count members, initialize storage
	// int N = halos.getNumMembers();
	int N = 700;
	std::vector<double> periDist(N,0);
	std::vector<double> stabIx(N,0);

	// Get positions of all primaries in the system
	std::vector<double> primPos = halos.getSysData().getModel()->getPrimPos(0, halos.getSysDataPtr());
	tpat_body_data P2data(halos.getSysData().getPrimary(1));

	// Set up sim engine to propagate until a periapsis
	tpat_simulation_engine sim(halos.getSysDataPtr());
	sim.clearEvents();	// remove crash events
	double evtData[] {1};
	tpat_event periEvt(halos.getSysDataPtr(), tpat_event::APSE, 1, true, evtData);
	sim.addEvent(periEvt);

	// Propagate all family members until periapsis and record their distance
	for(int n = 0; n < N; n++){
		tpat_family_member_cr3bp m = halos.getMember(n);
		sim.runSim(m.getIC(), m.getTOF()*1.1);	// Integrate a little farther than one rev to check the IC
		std::vector<eventRecord> eventsOccured = sim.getEventRecords();

		if(eventsOccured.size() == 0){
			printWarn("Simulation did not encounter a periapsis... saving NAN instead\n");
			periDist[n] = NAN;
		}else{
			std::vector<tpat_event> allEvents = sim.getEvents();
			for(size_t e = 0; e < eventsOccured.size(); e++){
				tpat_event event = allEvents[eventsOccured[e].eventIx];
				if(event.getType() == tpat_event::APSE){
					tpat_traj_cr3bp traj = sim.getCR3BP_Traj();
					std::vector<double> state = traj.getState(eventsOccured[e].stepIx);
					double dx = state[0] - primPos[1*3 + 0];
					double dy = state[1] - primPos[1*3 + 1];
					double dz = state[2] - primPos[1*3 + 2];
					periDist[n] = std::sqrt(dx*dx + dy*dy + dz*dz) * halos.getSysData().getCharL() - P2data.getRadius();
				}
			}
		}

		std::vector<cdouble> eigVals = m.getEigVals();
		stabIx[n] = getStabilityIndex(eigVals);
		// stabIx[n] = std::real( 0.5*(eigVals[0] + eigVals[1]) );

		printf("Orbit %03d: Stability = %.4f, Peri Dist = %.4f km\n", n, stabIx[n], periDist[n]);
	}

	printf("\n**********************\nSaving data to MAT file\n**********************\n");
	mat_t *matfp = Mat_CreateVer(outputData, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		size_t dims[] {(size_t)N,1};
		matvar_t *periVar = Mat_VarCreate("PeriAlt", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(periDist[0]), MAT_F_DONT_COPY_DATA);
		saveVar(matfp, periVar, "PeriAlt", MAT_COMPRESSION_NONE);

		matvar_t *stabilityVar = Mat_VarCreate("StabilityIx", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(stabIx[0]), MAT_F_DONT_COPY_DATA);
		saveVar(matfp, stabilityVar, "StabilityIx", MAT_COMPRESSION_NONE);
	}

	Mat_Close(matfp);
}