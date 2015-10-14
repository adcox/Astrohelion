/*
 *	Compute the location of the saddle point using only Sun and Earth, then compute it with Sun,
 *	Earth, and one additional body. Compare the difference to quantify how much a specific body perturbs
 *	the position of the saddle point.
 *	
 *	Analyzes perturbations from: Moon, Mars, Jupiter, Saturn, Venus, Mercury, Uranus, Neptune
 */
 
#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_utilities.hpp"

#include "cspice/SpiceUsr.h"
#include "gsl/gsl_linalg.h"

#include <cstdlib>

tpat_matrix getSPLoc_Ephm(double, std::vector<int>);

int main(){
	char timeKernel[] = "/Users/andrew/Documents/Purdue/Astrodynamics_Research/Code/C++/libTPAT/share/data_SPICE/naif0010.tls.pc";
	char planetKernel[] = "/Users/andrew/Documents/Purdue/Astrodynamics_Research/Code/C++/libTPAT/share/data_SPICE/de430.bsp";
	char massKernel[] = "/Users/andrew/Documents/Purdue/Astrodynamics_Research/Code/C++/libTPAT/share/data_SPICE/gm_de431.tpc";
	
	// Load the kernels
	furnsh_c(timeKernel);
    checkAndReThrowSpiceErr("Could not load SPICE time Kernel");
	furnsh_c(planetKernel);
	checkAndReThrowSpiceErr("Could not load SPICE planet Kernel");
    furnsh_c(massKernel);
    checkAndReThrowSpiceErr("Could not load SPICE mass Kernel");

    std::vector<int> allPerturbing = {301, 4, 5, 6, 299, 199, 7, 8};
    
    /*
     *	Now, do stuff to compute dates and data...
     */
    std::vector<int> bodies = {10, 399, 301};
    try{
	    tpat_matrix sp_pos = getSPLoc_Ephm(0, bodies);
	    printf("SP Position =\n");
	    sp_pos.print();
	}catch(tpat_diverge &e){
		printErr("Could not converge on SP location...\n");
	}

	double et0 = dateToEpochTime("2015/01/01");
	double et = 0;
	std::vector<double> data;
	for(double d = 0; d < 1500; d += 0.5){
		et = et0 + d*24*3600;	// Approx. Jan 2015 + specified days
		data.push_back(et);

		std::vector<int> bodies;
		tpat_matrix SE_spPos = getSPLoc_Ephm(et, bodies);	// If this throws an error (i.e. diverge), don't catch it

		for(size_t p = 0; p < allPerturbing.size(); p++){
			bodies.clear();
			bodies.insert(bodies.end(), allPerturbing.at(p));

			tpat_matrix sp_pos(3,1);
			try{
				sp_pos = getSPLoc_Ephm(et, bodies);

				// Compute distance between this perturbed SP location and the SE one
				data.push_back(norm(sp_pos - SE_spPos));
			}catch(tpat_diverge &e){
				printErr("SP location diverged at d = %.1f, p = %d\n", et, p);
				data.push_back(NAN);
			}

		}
	}

	// Save the data
	size_t cols = 1+allPerturbing.size();
	saveMatrixToFile("data/SPEffects.mat", "data", data, data.size()/cols, cols);

    // Unload the kernels
    unload_c(timeKernel);
    checkAndReThrowSpiceErr("Could not unload SPICE time Kernel");
    unload_c(planetKernel);
    checkAndReThrowSpiceErr("Could not unload SPICE planet Kernel");
    unload_c(massKernel);
    checkAndReThrowSpiceErr("Could not unload SPICE mass Kernel");
}//=====================================

/**
 *  @brief Compute the location of the saddle point using ephemeris data to locate all 
 *  primary masses
 * 
 *  @param et Epoch Time (seconds since J2000)
 *  @param bodyIx a vector of SPICE ID codes for the bodies to be considered when computing
 *  the location of the saddle point
 * 
 *  @return a three-element matrix representing the position of the saddle point relative to 
 *  the Sun in a J2000 Ecliptic reference frame 
 */
tpat_matrix getSPLoc_Ephm(double et, std::vector<int> bodyIx){

	// printf("Computing SP Location\n");
	// Get positions of all bodies relative to SSB
	ConstSpiceChar *obs = "SUN";		// observer: Sun
	ConstSpiceChar *abcorr = "NONE";	// Do not apply aberration corrections
	ConstSpiceChar *ref = "ECLIPJ2000";	// reference frame: J2000, Ecliptic plane
	ConstSpiceChar *GM_item = "GM";		// Retrieve GM for each body

	// Make sure sun and earth are first two
	bodyIx.insert(bodyIx.begin(), 399);
	bodyIx.insert(bodyIx.begin(), 10);

	std::vector<tpat_matrix> positions;
	std::vector<double> GMs;
	for(size_t b = 0; b < bodyIx.size(); b++){
		// skip Sun and Earth if they aren't the first two
		if(b > 1 && (bodyIx[b] == 10 || bodyIx[b] == 399))
			continue;
		
		SpiceDouble pos[3];			// position of the body
		SpiceDouble lt = 0;			// light time for the body; not needed
		SpiceChar targ[128];
		sprintf(targ, "%d", bodyIx[b]);	// turn the ID into a string

		// Get the position of the body
		spkpos_c(targ, et, ref, abcorr, obs, pos, &lt);	
		checkAndReThrowSpiceErr("Could not compute body position");

		SpiceInt numVal = 0;
		SpiceDouble GM = 0;

		bodvrd_c(targ, GM_item, 1, &numVal, &GM);
		checkAndReThrowSpiceErr("Could not get body GM value");

		if(numVal == 0){
			throw tpat_exception("No SPICE error, but zero GM values were returned for the body");
		}

		// Create a 3-element column vector for the position vector, save to vector of positions
		positions.push_back(tpat_matrix(3, 1, pos));
		GMs.push_back(GM);
		// printf("  Obtained data for %s\n", getNameFromSpiceID(bodyIx[b]).c_str());
		// printf("    Pos = [%.4f, %.4f, %.4f] km\n", pos[0], pos[1], pos[2]);
		// printf("    GM = %.4e\n", GM);
	}

	// Approximate SP location using only the sun and earth
	double mu3B = GMs[1]/(GMs[0] + GMs[1]);
	double x = (2*mu3B - 2*mu3B*mu3B + sqrt(mu3B*(1-mu3B)) - 1)/(2*mu3B - 1);
	tpat_matrix sp_pos = positions[0] + x*(positions[1] - positions[0]);

	double err = 1e20, okErr = 1e-12;
	int count = 0, maxIts = 20;

	while(err > okErr && count < maxIts){
		tpat_matrix accel(3,1);
		double dfxdx = 0, dfydy = 0, dfzdz = 0, dfxdy = 0, dfxdz = 0, dfydz = 0;

		for(size_t b = 0; b < positions.size(); b++){
			tpat_matrix p = positions[b];
			tpat_matrix r = sp_pos - p;	// Vector from SP to bodies
			double rMag = norm(r);

			// Compute acceleration due to this body
			accel += GMs[b] * r/pow(rMag,3);

			// Compute this body's contributions to the partials w.r.t. position states
			dfxdx += GMs[b] * (pow(rMag, -3) - 3*r.at(0)*r.at(0)/pow(rMag, 5));
			dfydy += GMs[b] * (pow(rMag, -3) - 3*r.at(1)*r.at(1)/pow(rMag, 5));
			dfzdz += GMs[b] * (pow(rMag, -3) - 3*r.at(2)*r.at(2)/pow(rMag, 5));
			dfxdy += GMs[b] * -3*r.at(0)*r.at(1)/pow(rMag,5);
			dfxdz += GMs[b] * -3*r.at(0)*r.at(2)/pow(rMag,5);
			dfydz += GMs[b] * -3*r.at(1)*r.at(2)/pow(rMag,5);
		}

		double J_data[] = {dfxdx, dfxdy, dfxdz, dfxdy, dfydy, dfydz, dfxdz, dfydz, dfzdz};
		tpat_matrix J(3,3,J_data);
		J *= -1;	// multiply by -1 for solution process

		// Solve Equation Jx = b, where b is acceleration and x is difference between
		// next and current SP position vector
		gsl_vector_view b = gsl_vector_view_array(accel.getDataPtr(), accel.getRows());

		// Allocate memory for intermediate vector w
		gsl_vector *w = gsl_vector_alloc(3);
		gsl_permutation *perm = gsl_permutation_alloc(J.getRows());
		int permSign = 0;
		int status = gsl_linalg_LU_decomp(J.getGSLMat(), perm, &permSign);
		if(status){
			printErr("GSL ERR: %s\n", gsl_strerror(status));
			gsl_permutation_free(perm);
			gsl_vector_free(w);
			throw tpat_linalg_err("Unable to decompose J into L and U");
		}
		status = gsl_linalg_LU_solve(J.getGSLMat(), perm, &(b.vector), w);
		if(status){
			printErr("GSL ERR: %s\n", gsl_strerror(status));
			gsl_permutation_free(perm);
			gsl_vector_free(w);
			throw tpat_linalg_err("Unable to invert J, likely singular");
		}

		// w, in this case, is x2 - x_1
		tpat_matrix diff(w, false);

		sp_pos += diff;		// Update SP position
		err = norm(accel);	// Update error
		printColor(YELLOW, "Iterations %02d: ||A|| = %.8e\n", count, err);

		count++;
	}

	if(err < okErr)
		return sp_pos;
	else
		throw tpat_diverge("Could not converge on SP position");
}