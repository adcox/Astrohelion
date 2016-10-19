#define BOOST_TEST_MODULE BC4BP_Def_Test

#include <boost/test/unit_test.hpp>
#include <cspice/SpiceUsr.h>

#include "Calculations.hpp"
#include "DynamicsModel_bc4bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

/**
 *  @brief [brief description]
 *  @details [long description]
 * 
 *  @param pBCSys SEM BC4BP System Data
 *  @param pSESys SE CR3BP System Data
 *  @param t time in nondimensional SEM time units
 *  @return inertial, dimensional position of Moon in ecliptic J2000 frame relative to Sun, km
 */
std::vector<double> getInertMoonPos(SysData_bc4bp *pBCSys, SysData_cr3bp *pSESys, double t){
	double bcPrimPos[9];
	DynamicsModel_bc4bp::getPrimaryPos(t, pBCSys, bcPrimPos);

	double sePrimPos[9];
	std::copy(bcPrimPos, bcPrimPos + 9, sePrimPos);

	// Shift BC4BP coordinates from EM barycenter to SEM barycenter
	for(unsigned int i = 0; i < 3; i++){
		sePrimPos[3*i] += (1.0/pBCSys->getK() - pBCSys->getMu());
	}

	// Copy the ones we want
	std::vector<double> seMoonPos = {sePrimPos[6], sePrimPos[7], sePrimPos[8], 0, 0, 0};

	// Dimensionalize and re-nondimensionalize for the SE system
	for(unsigned int i = 0; i < 3; i++){
		seMoonPos[i] *= pBCSys->getCharL()/pSESys->getCharL();
	}

	return cr3bp_rot2inert_state(seMoonPos, pSESys, (t*pBCSys->getCharT())/pSESys->getCharT(), pBCSys->getEpoch0(), 1);

	// Conversion from EM yields the EXACT same answer (to about machine precision)
	// SysData_cr3bp emSys("earth", "moon");
	// std::vector<double> moonPos_rot {1-emSys.getMu(), 0, 0, 0, 0, 0};
	// return cr3bp_rot2inert_state(moonPos_rot, &emSys, t*pBCSys->getCharT()/emSys.getCharT(), pBCSys->getEpoch0(), 0);
}//====================================================

std::vector<double> getSpiceMoonPos(double et){
	// Now, compute Moon location from ephemeris data
    ConstSpiceChar *abcorr = "none";
    ConstSpiceChar *ref = "ECLIPJ2000";
    ConstSpiceChar *obs = "Sun";	// Maybe need to try SSB?
    ConstSpiceChar *targ = "MOON";

    SpiceDouble lt = 0;
    SpiceDouble ephmState[6];

    spkezr_c(targ, et, ref, abcorr, obs, ephmState, &lt);
    std::vector<double> moonPos(ephmState, ephmState+6);

    // unload_c(timeKernel);
    // unload_c(deKernel);
    return moonPos;
}//====================================================

BOOST_AUTO_TEST_SUITE(Linkable)

BOOST_AUTO_TEST_CASE(PrimaryPositions){
	SysData_bc4bp bcSys("sun", "earth", "moon");
	SysData_cr3bp seSys("sun", "earth");

	// double epoch = dateToEpochTime("2016/01/01");
	// double epoch = SysData_bc4bp::REF_EPOCH;
	// double epoch = 0;
	// DynamicsModel_bc4bp::orientAtEpoch(epoch, &bcSys);
	double epoch = bcSys.getEpoch0();
	printf("Epoch0 = %.4f sec ET\n", epoch);

	// double t1 = 24*3600, t2 = 7*24*3600, t3 = 12*7*24*3600;		// Times to compute primary positions at, seconds past reference epoch
	double t1 = 0;
	double t2 = 1*24*3600;
	double t3 = 2*24*3600;
	double epoch1 = epoch + t1, epoch2 = epoch + t2, epoch3 = epoch + t3;	// Epochs to compute primary positions at, epoch time, seconds

	std::vector<double> moonPos_myConversion;
	std::vector<double> moonPos_spice;

	for(double t = 0; t < 100*24*3600; t += 4*3600){
		std::vector<double> me = getInertMoonPos(&bcSys, &seSys, t/bcSys.getCharT());
		std::vector<double> spice = getSpiceMoonPos(epoch + t);
		moonPos_myConversion.insert(moonPos_myConversion.end(), me.begin(), me.begin()+3);
		moonPos_spice.insert(moonPos_spice.end(), spice.begin(), spice.begin()+3);
	}

	// Convert Positions to Sun-centered inertial coordinates
	// std::vector<double> sciMoonPos1 = getInertMoonPos(&bcSys, &seSys, t1/bcSys.getCharT());
	// std::vector<double> sciMoonPos2 = getInertMoonPos(&bcSys, &seSys, t2/bcSys.getCharT());
	// std::vector<double> sciMoonPos3 = getInertMoonPos(&bcSys, &seSys, t3/bcSys.getCharT());

	// std::vector<double> ephmState1 = getSpiceMoonPos(epoch1);
	// std::vector<double> ephmState2 = getSpiceMoonPos(epoch2);
	// std::vector<double> ephmState3 = getSpiceMoonPos(epoch3);

	// printf("Epoch 1: %.4f sec\n", epoch1);
 //    printf("Moon Position 1 from BC4BP: [%.4f, %.4f, %.4f] km\n", sciMoonPos1[0], sciMoonPos1[1], sciMoonPos1[2]);
 //    printf("Moon Position 1 from SPICE: [%.4f, %.4f, %.4f] km\n", ephmState1[0], ephmState1[1], ephmState1[2]);

 //    printf("Epoch 2: %.4f sec\n", epoch2);
 //    printf("Moon Position 2 from BC4BP: [%.4f, %.4f, %.4f] km\n", sciMoonPos2[0], sciMoonPos2[1], sciMoonPos2[2]);
 //    printf("Moon Position 2 from SPICE: [%.4f, %.4f, %.4f] km\n", ephmState2[0], ephmState2[1], ephmState2[2]);

 //    printf("Epoch 3: %.4f sec\n", epoch3);
 //    printf("Moon Position 3 from BC4BP: [%.4f, %.4f, %.4f] km\n", sciMoonPos3[0], sciMoonPos3[1], sciMoonPos3[2]);
 //    printf("Moon Position 3 from SPICE: [%.4f, %.4f, %.4f] km\n", ephmState3[0], ephmState3[1], ephmState3[2]);

	saveMatrixToFile("moonPos_me.mat", "moonPos_custom", moonPos_myConversion, moonPos_myConversion.size()/3, 3);
	saveMatrixToFile("moonPos_spice.mat", "moonPos_spice", moonPos_spice, moonPos_spice.size()/3, 3);
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()