/**
 *	Generate a grid of initial conditions inside a manifold tube and see where they go near the moon!
 */

#include "tpat_calculations.hpp"
#include "tpat_event.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <vector>

// *********************
// Function Declarations
// *********************
std::vector<double> generateManifoldTube(tpat_sys_data_cr3bp*, double, int);
std::vector<double> getRange(int, int, double, std::vector<double>);
std::vector<double> createICGrid(tpat_sys_data_cr3bp*, std::vector<double>, double, int, int, int, int);

// *********************
// Function Definitions
// *********************

int main(){
	tpat_sys_data_cr3bp emSys("earth", "moon");

	for(double C = 3.14; C > 3.01; C -= 0.01){
		printf("Computing manifolds...\n");
		std::vector<double> crossings = generateManifoldTube(&emSys, C, 400);
		char crossFile[64];
		sprintf(crossFile, "data/crossings_C%.3f.mat", C);
		saveMatrixToFile(crossFile, "crossings", crossings, crossings.size()/6, 6);

		printf("Creating Grid of ICs...\n");
		std::vector<double> allIC = createICGrid(&emSys, crossings, C, 100, 50, 50, 50);

		char gridFile[64];
		sprintf(gridFile, "data/gridIC_C%.3f.mat", C);
		saveMatrixToFile(gridFile, "gridICs", allIC, allIC.size()/6, 6);
	}
	return EXIT_SUCCESS;
}//=========================================================

/**
 *  @brief Compute a manifold tube from Earth-Moon L2 Southern Halo family
 * 
 *  @param emSys pointer to Earth-Moon system data object
 *  @param JC desired Jacobi constant for the halo/manifolds
 *  @param numManifolds number of manifolds
 *  @return a nx6 row-major-order vector of states on a hyperplane crossing
 */
std::vector<double> generateManifoldTube(tpat_sys_data_cr3bp *emSys, double JC, int numManifolds){
	
	tpat_family_cr3bp haloFam("../share/families_natParam_checked/EM_L2_NHalo.mat");

	std::vector<tpat_family_member_cr3bp> halos = haloFam.getMemberByJacobi(JC);

	std::vector<double> crossings;

	if(halos.size() < 1){
		printErr("Could not find any halos with a JC value of %.4f... Exiting\n", JC);
		return crossings;
	}

	std::vector<double> haloIC = halos[0].getIC();
	// Flip the Halo to be a southern guy
	haloIC[2] = -haloIC[2];	// z
	haloIC[5] = -haloIC[5];	// z-dot

	// Integrate ICs to get a full orbit
	tpat_simulation_engine sim(emSys);
	sim.runSim(haloIC, halos[0].getTOF());
	tpat_traj_cr3bp aHalo = sim.getCR3BP_Traj();

	// Use a very short integration time to get the ICs for a bunch of halo manifolds
	std::vector<tpat_traj_cr3bp> manifolds = getManifolds(MAN_S_P, &aHalo, numManifolds, 1e-5);
	double LPt[3];
	cr3bp_getEquilibPt(*emSys, 2, 1e-14, LPt);

	// Stop integration at an XZ plane just outside the L2 point
	double evtData[] = {LPt[1] + 0.3, NAN, NAN, NAN, NAN, NAN};
	tpat_event endEvt(emSys, tpat_event::XZ_PLANE, 0, true, evtData);	// Stop at y = just outside L2

	sim.setRevTime(true);	// running backwards in time for Stable manifolds
	sim.addEvent(endEvt);

	crossings.reserve(manifolds.size()*6);
	std::vector<double> ics;
	ics.reserve(manifolds.size()*7);
	std::vector<double> blank {NAN, NAN, NAN, NAN, NAN, NAN};

	for(size_t m = 0; m < manifolds.size(); m++){

		std::vector<double> IC = manifolds[m].getState(0);
		sim.runSim(IC, 4*PI);

		std::vector<tpat_event> endEvents = sim.getEndEvents();
		tpat_traj_cr3bp traj = sim.getCR3BP_Traj();

		ics.insert(ics.end(), IC.begin(), IC.end());
		ics.push_back(traj.getTime(-1));

		if(std::find(endEvents.begin(), endEvents.end(), endEvt) != endEvents.end()){
			// printf("Succesfully found plane crossing for manifold #%02zu\n", m);
			std::vector<double> lastState = traj.getState(-1);
			crossings.insert(crossings.end(), lastState.begin(), lastState.end());
		}else{
			crossings.insert(crossings.end(), blank.begin(), blank.end());
			printf("Did not find plane crossing for manifold #%02zu\n", m);
		}
	}

	saveMatrixToFile("data/ManCross_IC.mat", "ICs", ics, ics.size()/7, 7);
	return crossings;
}//=========================================================

/**
 *  @brief Determine the range of a dependent variable at a specific value of an independent variable
 *  on a manifold crossing contour
 *  @details [long description]
 * 
 *  @param indVarIx Index of the independent variable
 *  @param depVarIx Index of the dependent variable
 *  @param indVarVal Value of the dependent variable
 *  @param crossings an nx6 row-major-order listing of the manifold crossings at a map
 *  @return the minimum and maximum range of the dependent variable at the specified value
 *  of the independent variable
 */
std::vector<double> getRange(int indVarIx, int depVarIx, double indVarVal, std::vector<double> crossings){
	/* Store four points: {i, d} where i is the independent variable, d is the dependent variable
	 * "left" indicates that i < indVarVal
	 * "right" indicates that i > indVarVal
	 * "max" indicates that d is maximized
	 * "min" inidciates that d is minimized
	 */ 
	double leftMax[] = {-1e20, -1e20};
	double leftMin[] = {-1e20, 1e20};
	double rightMax[] = {1e20, -1e20};
	double rightMin[] = {1e20, 1e20};
	int maxTries = 25;
	// Compute the distance from each crossing to the desired independent variable value
	std::vector<double> horizDist;
	for(size_t i = 0; i < crossings.size()/6.0; i++){
		horizDist.push_back(std::abs(indVarVal - crossings[i*6 + indVarIx]));
	}

	// Compute an average for the dependent variable
	double sumDepVals = 0;
	std::vector<double> distCpy = horizDist;
	for(int i = 0; i < maxTries; i++){
		// Find the closest point
		std::vector<double>::iterator min = std::min_element(distCpy.begin(), distCpy.end());
		int minIx = min - distCpy.begin();	// get its index

		sumDepVals += crossings[minIx*6 + depVarIx];	// Add the value of the dependent variable to the sum

		distCpy[minIx] = 1e20;	// make sure this point isn't used as minimum again
	}
	double meanDepVal = sumDepVals/((double)maxTries);	// Compute the average value of the dependent variable here

	bool LT = false, LB = false, RT = false, RB = false;	// Track which points we've found
	int count = 0;
	while(count < maxTries && (!LT || !LB || !RT || !RB)){
		std::vector<double>::iterator min = std::min_element(horizDist.begin(), horizDist.end());
		int minIx = min - horizDist.begin();
		double signedDist = indVarVal - crossings[minIx*6 + indVarIx];

		if(signedDist > 0){
			// Left Side
			if(crossings[minIx*6 + depVarIx] > meanDepVal && !LT){
				leftMax[0] = crossings[minIx*6 + indVarIx];
				leftMax[1] = crossings[minIx*6 + depVarIx];
				LT = true;
			}else if(crossings[minIx*6 + depVarIx] < meanDepVal && !LB){
				leftMin[0] = crossings[minIx*6 + indVarIx];
				leftMin[1] = crossings[minIx*6 + depVarIx];
				LB = true;
			}
		}else{
			// Right Side
			if(crossings[minIx*6 + depVarIx] > meanDepVal && !RT){
				rightMax[0] = crossings[minIx*6 + indVarIx];
				rightMax[1] = crossings[minIx*6 + depVarIx];
				RT = true;
			}else if(crossings[minIx*6 + depVarIx] < meanDepVal && !RB){
				rightMin[0] = crossings[minIx*6 + indVarIx];
				rightMin[1] = crossings[minIx*6 + depVarIx];
				RB = true;
			}
		}

		horizDist[minIx] = 1e20;	// Make sure this point isn't picked as a minimum again
		count++;
	}

	std::vector<double> range = {NAN, NAN};

	if( std::abs(leftMax[0]) == 1e20 || std::abs(rightMax[0]) == 1e20 ||
		std::abs(leftMin[0]) == 1e20 || std::abs(rightMin[0]) == 1e20){
		// // Print outs for debugging
		// printErr("getRange: Unable to locate four points\n");
		// printf("  leftMax = [%.4f, %.4f]\n", leftMax[0], leftMax[1]);
		// printf("  rightMax = [%.4f, %.4f]\n", rightMax[0], rightMax[1]);
		// printf("  leftMin = [%.4f, %.4f]\n", leftMin[0], leftMin[1]);
		// printf("  rightMin = [%.4f, %.4f]\n", rightMin[0], rightMin[1]);
	}else{
		// Interpolate between points to get a guess for actual bounds
		double slopeMax = (rightMax[1] - leftMax[1])/(rightMax[0] - leftMax[0]);
		double slopeMin = (rightMin[1] - leftMin[1])/(rightMin[0] - leftMin[0]);
		range[0] = leftMin[1] + slopeMin*(indVarVal - leftMin[0]);	// minimum
		range[1] = leftMax[1] + slopeMax*(indVarVal - leftMax[0]);	// maximum
	}

	return range;
}

std::vector<double> createICGrid(tpat_sys_data_cr3bp *emSys, std::vector<double> crossings, double C,
	int numXSteps, int numZSteps, int numVxSteps, int numVySteps){
	
	std::vector<double> allIC;
	allIC.reserve(numXSteps*numZSteps*numVxSteps*numVySteps);
	
	// Determine range of y for this manifold crossing
	std::vector<double> allX;
	for(size_t i = 0; i < crossings.size()/6; i++){
		allX.push_back(crossings[i*6 + 0]);
	}

	std::vector<double>::iterator minX = std::min_element(allX.begin(), allX.end());
	std::vector<double>::iterator maxX = std::max_element(allX.begin(), allX.end());

	// Bound the search in y
	double search_xStep = (*maxX - *minX)/((double)(numXSteps+1));
	double x = 0;
	for(int ix = 0; ix < numXSteps; ix++){
		x = *minX + (ix+1)*search_xStep;
		printf("  X step %03d/%03d\n", ix, numXSteps-1);

		// Get ranges for z, vx, and vy
		std::vector<double> zBounds = getRange(0,2,x,crossings);
		
		if(zBounds[0] != zBounds[0]){
			printErr("Could not find bounds for z at x = %.4f\n", x);
			continue;
		}

		double search_zStep = (zBounds[1] - zBounds[0])/((double)(numZSteps + 1));

		double ic[] = {x, crossings[1], 0, 0, 0, 0};
		for(int iz = 0; iz < numZSteps; iz++){
			// Choose Z
			ic[2] = zBounds[0] + (iz+1)*search_zStep;

			// Compute bounds for Vx using the chosen x (1) and z (2) coordinates
			std::vector<double> vxBounds_x = getRange(0,3,ic[0],crossings);
			std::vector<double> vxBounds_z = getRange(2,3,ic[2],crossings);

			if(vxBounds_x[0] != vxBounds_x[0] || vxBounds_z[0] != vxBounds_z[0]){
				printErr("Could not find bounds for vx at (x,z) = (%.4f, %.4f)\n", ic[0], ic[2]);
				continue;
			}

			double search_vxStart = std::max(vxBounds_x[0], vxBounds_z[0]);
			double search_vxEnd = std::min(vxBounds_x[1], vxBounds_z[1]);
			double search_vxStep = (search_vxEnd - search_vxStart)/((double)(numVxSteps + 1));

			for(int ivx = 0; ivx < numVxSteps; ivx++){
				ic[3] = search_vxStart + (ivx+1)*search_vxStep;

				// Compute bounds for Vy using the chosen x (1) and z (2) coordinates
				std::vector<double> vyBounds_x = getRange(0,4,ic[0],crossings);
				std::vector<double> vyBounds_z = getRange(2,4,ic[2],crossings);
				std::vector<double> vyBounds_vx = getRange(3,4,ic[3],crossings);

				if(vyBounds_x[0] != vyBounds_x[0] || vyBounds_z[0] != vyBounds_z[0] ||
					vyBounds_vx[0] != vyBounds_vx[0]){
					
					// printErr("Could not find bounds for vx at (x,z,vx) = (%.4f, %.4f, %.4f)\n", ic[1], ic[2], ic[3]);
					continue;
				}

				double search_vyStart = std::max(std::max(vyBounds_x[0], vyBounds_z[0]), vyBounds_vx[0]);
				double search_vyEnd = std::min(std::min(vyBounds_x[1], vyBounds_z[1]), vyBounds_vx[1]);
				double search_vyStep = (search_vyEnd - search_vyStart)/((double)(numVySteps + 1));

				for(int ivy = 0; ivy < numVySteps; ivy++){
					ic[4] = search_vyStart + (ivy+1)*search_vyStep;

					// Given the choice of z, determine acceptable range of values for vz
					std::vector<double> vzBounds_x = getRange(0,5,ic[0],crossings);
					// std::vector<double> vzBounds_z = vzBounds_x;
					// std::vector<double> vzBounds_vx = vzBounds_x;
					// std::vector<double> vzBounds_vy = vzBounds_x;
					std::vector<double> vzBounds_z = getRange(2,5,ic[2],crossings);
					std::vector<double> vzBounds_vx = getRange(3,5,ic[3],crossings);
					std::vector<double> vzBounds_vy = getRange(4,5,ic[4],crossings);

					if(vzBounds_z[0] != vzBounds_z[0] || vzBounds_x[0] != vzBounds_x[0] || 
						vzBounds_vx[0] != vzBounds_vx[0] || vzBounds_vy[0] != vzBounds_vy[0]){

						// printErr("Could not find vz range for y = %.4f\n", y);
						continue;
					}

					double maxZ = std::min(std::min(std::min(vzBounds_z[1], vzBounds_x[1]), 
						vzBounds_vx[1]), vzBounds_vy[1]);
					double minZ = std::max(std::max(std::max(vzBounds_z[0], vzBounds_x[0]), 
						vzBounds_vx[0]), vzBounds_vy[0]);

					ic[5] = cr3bp_getVel_withC(ic, emSys->getMu(), C, 5);

					// Don't allow NAN
					if(ic[5] != ic[5])
						continue;

					// Check to make sure the computed vz falls within the bounds
					if(ic[5] < maxZ && ic[5] > minZ){
						allIC.insert(allIC.end(), ic, ic+6);
					}else if(-ic[5] < maxZ && -ic[5] > minZ){
						ic[5] *= -1;
						allIC.insert(allIC.end(), ic, ic+6);
					}
				}
			}
		}
	}

	return allIC;
}