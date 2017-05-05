#define BOOST_TEST_MODULE BC4BP_Def_Test

#include <boost/test/unit_test.hpp>
#include <cspice/SpiceUsr.h>

#include "Calculations.hpp"
#include "DynamicsModel_bc4bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Node.hpp"
#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

/**
 *  \brief Get the inertial moon position 
 *  \details The Moon's location is obtained by computing the moon position in the Sun-Earth-Moon BC4BP,
 *  then shifting and scaling to the Sun-Earth system, then converting to Earth-Moon system, and finally,
 *  to the Earth-centered J2000 frame.
 * 
 *  \param pBCSys SEM BC4BP System Data
 *  \param pSESys SE CR3BP System Data
 *  \param t time in nondimensional SEM time units
 *  \return inertial, dimensional position of Moon in ecliptic J2000 frame relative to Sun, km
 */
std::vector<double> getInertMoonPos(SysData_bc4bp *pBCSys, SysData_cr3bp *pSESys, double t){
	double bcPrimPos[9] = {0};
	pBCSys->getDynamicsModel()->getPrimPos(t, pBCSys, -1, bcPrimPos);

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

	// return cr3bp_rot2inert_state(seMoonPos, pSESys, (t*pBCSys->getCharT())/pSESys->getCharT(), pBCSys->getEpoch0(), 2);

	// Convert to EM coordinates
	SysData_cr3bp emSys("earth", "moon");
	std::vector<double> emMoonPos = cr3bp_SE2EM_state(seMoonPos, (t*pBCSys->getCharT())/pSESys->getCharT(),
		pBCSys->getTheta0(), pBCSys->getPhi0(), pBCSys->getGamma(), emSys.getCharL(), emSys.getCharT(),
		pSESys->getCharL(), pSESys->getCharT(), pSESys->getMu());

	return cr3bp_rot2inert_state(emMoonPos, &emSys, (t*pBCSys->getCharT())/emSys.getCharT(), pBCSys->getEpoch0(), 1);

	// Conversion from EM yields the EXACT same answer (to about machine precision)
	// SysData_cr3bp emSys("earth", "moon");
	// std::vector<double> moonPos_rot {1-emSys.getMu(), 0, 0, 0, 0, 0};
	// return cr3bp_rot2inert_state(moonPos_rot, &emSys, t*pBCSys->getCharT()/emSys.getCharT(), pBCSys->getEpoch0(), 1);
}//====================================================

std::vector<double> getInertMoonPos_2(SysData_bc4bp *pBCSys, SysData_cr3bp *pSESys, double t){
	SysData_cr3bp emSys("earth", "moon");

	double bcPrimPos[9] = {0};
	pBCSys->getDynamicsModel()->getPrimPos(t, pBCSys, -1, bcPrimPos);

	std::vector<double> bcMoonPos = {bcPrimPos[6], bcPrimPos[7], bcPrimPos[8], 0, 0, 0};

	// Create a node
	Node semNode(bcMoonPos, t);
	Node semNode2(bcMoonPos, t + 3600/pBCSys->getCharT());

	Arcset_bc4bp semNodes(pBCSys);
	int id1 = semNodes.addNode(semNode);
	int id2 = semNodes.addNode(semNode2);
	semNodes.addSeg(Segment(id1, id2, semNode2.getEpoch() - semNode.getEpoch()));

	Arcset_cr3bp seNodes = bcr4bpr_SEM2SE(semNodes, pSESys);
	
	Arcset_cr3bp eci_fromSE = cr3bp_rot2inert(seNodes, pBCSys->getEpoch0(), 2);
	Arcset_cr3bp emNodes = cr3bp_SE2EM(seNodes, &emSys, pBCSys->getTheta0(), pBCSys->getPhi0(), pBCSys->getGamma());
	Arcset_cr3bp eci_fromEM = cr3bp_rot2inert(emNodes, pBCSys->getEpoch0(), 1);

	// return eci_fromSE.getStateByIx(0);
	return eci_fromEM.getStateByIx(0);
}//====================================================

std::vector<double> getSpiceMoonPos(double et){
	// Now, compute Moon location from ephemeris data
    ConstSpiceChar *abcorr = "none";
    ConstSpiceChar *ref = "ECLIPJ2000";
    // ConstSpiceChar *obs = "Sun";	// Maybe need to try SSB?
    ConstSpiceChar *obs = "EARTH";	// Compare to my code transformation to Earth-centered inertial
    ConstSpiceChar *targ = "MOON";

    SpiceDouble lt = 0;
    SpiceDouble ephmState[6];

    spkezr_c(targ, et, ref, abcorr, obs, ephmState, &lt);
    checkAndReThrowSpiceErr("getSpiceMoonPos spkezr_c error");
    std::vector<double> moonPos(ephmState, ephmState+6);

    return moonPos;
}//====================================================

BOOST_AUTO_TEST_SUITE(BC4BP_Primary_Positions)

BOOST_AUTO_TEST_CASE(PrimaryPositions){
	SysData_bc4bp bcSys("sun", "earth", "moon");
	SysData_cr3bp seSys("sun", "earth");

	double epoch = dateToEphemerisTime("2016/01/01");
	// double epoch = SysData_bc4bp::REF_EPOCH + 4*24*3600;
	DynamicsModel_bc4bp::orientAtEpoch(epoch, &bcSys);
	// double epoch = bcSys.getEpoch0();
	printf("Epoch0 = %.4f sec ET\n", epoch);
	printf("Theta0 = %.4f deg\n", bcSys.getTheta0()*180/PI);
	printf("Phi0 = %.4f deg\n", bcSys.getPhi0()*180/PI);
	printf("Gamma = %.4f deg\n", bcSys.getGamma()*180/PI);

	std::vector<double> moonPos_myConversion;
	std::vector<double> moonPos_spice;
	std::vector<double> time;

	for(double t = 0; t < 35*365*24*3600; t += 4*3600){
		std::vector<double> me = getInertMoonPos(&bcSys, &seSys, t/bcSys.getCharT());
		std::vector<double> spice = getSpiceMoonPos(epoch + t);

		BOOST_CHECK(std::abs(me[0] - spice[0]) < 1000);
		BOOST_CHECK(std::abs(me[1] - spice[1]) < 2000);
		BOOST_CHECK(std::abs(me[2] - spice[2]) < 300);
		
		// moonPos_myConversion.insert(moonPos_myConversion.end(), me.begin(), me.begin()+3);
		// moonPos_spice.insert(moonPos_spice.end(), spice.begin(), spice.begin()+3);
		// time.push_back(t);
	}

	// saveMatrixToFile("moonPos_me.mat", "moonPos_custom", moonPos_myConversion, moonPos_myConversion.size()/3, 3);
	// saveMatrixToFile("moonPos_spice.mat", "moonPos_spice", moonPos_spice, moonPos_spice.size()/3, 3);
	// saveMatrixToFile("moonPos_times.mat", "time", time, time.size(), 1);
}

BOOST_AUTO_TEST_SUITE_END()