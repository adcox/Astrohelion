#define BOOST_TEST_MODULE Trajectory

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj_bc4bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Traj_cr3bp_lt.hpp"

using namespace astrohelion;

BOOST_AUTO_TEST_CASE(CR3BP_Save_Load){
	SysData_cr3bp emData("earth", "moon");
	SimEngine sim;
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	Traj_cr3bp crTraj(&emData);
	sim.runSim(ic, T, &crTraj);

	// Query the acceleration so it is computed
	std::vector<double> a = crTraj.getAccelByIx(-1);
	crTraj.saveToMat("data/crTraj.mat");

	Traj_cr3bp crTemp(&emData);
	crTemp.readFromMat("data/crTraj.mat");

	// printf("Testing Save/Read functions on CR3BP Trajectory\n");
	BOOST_CHECK(crTraj.getStateByIx(-1) == crTemp.getStateByIx(-1));
	BOOST_CHECK(crTraj.getAccelByIx(-1) == crTemp.getAccelByIx(-1));
	BOOST_CHECK(crTraj.getTimeByIx(-1) == crTemp.getTimeByIx(-1));
	BOOST_CHECK(crTraj.getSTMByIx(-1) == crTemp.getSTMByIx(-1));
	BOOST_CHECK(crTraj.getJacobiByIx(-1) == crTemp.getJacobiByIx(-1));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_Save_Load){
	SysData_bc4bp semData("sun", "earth", "moon");
	SimEngine sim;
	double ic[] = {-0.745230328320519, 7.22625684942683e-04, 7.45549413286038e-05, -7.30710697247992e-06, -0.0148897145134465, -1.23266135281459e-06};
	double T = 313;	// SE L1 Period
	Traj_bc4bp bcTraj(&semData);
	sim.runSim(ic, T, &bcTraj);

	// Query the acceleration so it is computed
	std::vector<double> a = bcTraj.getAccelByIx(-1);
	bcTraj.saveToMat("data/bcTraj.mat");

	Traj_bc4bp bcTemp(&semData);
	bcTemp.readFromMat("data/bcTraj.mat");

	// printf("Testing Save/Read functions on BC4BP Trajectory\n");
	BOOST_CHECK(bcTraj.getStateByIx(-1) == bcTemp.getStateByIx(-1));
	BOOST_CHECK(bcTraj.getAccelByIx(-1) == bcTemp.getAccelByIx(-1));
	BOOST_CHECK(bcTraj.getTimeByIx(-1) == bcTemp.getTimeByIx(-1));
	BOOST_CHECK(bcTraj.getSTMByIx(-1) == bcTemp.getSTMByIx(-1));
	BOOST_CHECK(bcTraj.get_dqdTByIx(-1) == bcTemp.get_dqdTByIx(-1));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Save_Load){
	SysData_cr3bp_lt emData("earth", "moon", 12e-3, 1500, 14);
	SimEngine sim;
	sim.setCtrlLaw(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	Traj_cr3bp_lt ltTraj(&emData);
	sim.runSim(ic, T, &ltTraj);

	// Query the acceleration so it is computed
	std::vector<double> a = ltTraj.getAccelByIx(-1);
	ltTraj.saveToMat("data/lowthrustTraj.mat");

	Traj_cr3bp_lt ltTemp(&emData);
	ltTemp.readFromMat("data/lowthrustTraj.mat");

	// printf("Testing Save/Read functions on CR3BP Trajectory\n");
	BOOST_CHECK(ltTraj.getStateByIx(-1) == ltTemp.getStateByIx(-1));
	BOOST_CHECK(ltTraj.getAccelByIx(-1) == ltTemp.getAccelByIx(-1));
	BOOST_CHECK(ltTraj.getTimeByIx(-1) == ltTemp.getTimeByIx(-1));
	BOOST_CHECK(ltTraj.getSTMByIx(-1) == ltTemp.getSTMByIx(-1));
	BOOST_CHECK(ltTraj.getJacobiByIx(-1) == ltTemp.getJacobiByIx(-1));	
}//====================================================