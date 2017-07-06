#define BOOST_TEST_MODULE Event

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <ctime>

#include "Common.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_cr3bp_lt.hpp"

using namespace astrohelion;
using namespace boost::unit_test;

std::vector<int> PIx {0, 1};

BOOST_AUTO_TEST_SUITE(EventAccuracy)

BOOST_AUTO_TEST_CASE(CR3BP_YZ_PLANE){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_XY_PLANE){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_XZ_PLANE){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_CRASH){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_JC){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_DIST){
	//TODO - Write this test
	BOOST_CHECK(true);
}//====================================================

BOOST_DATA_TEST_CASE(CR3BP_APSE, data::make(PIx), p){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo
	double tof = 40*PI;

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.setMakeDefaultEvents(false);

	std::vector<double> apseData {static_cast<double>(p)};
	Event evt(Event_tp::APSE, 0, false, apseData);
	sim.addEvent(evt);

	Arcset_cr3bp traj(&sys);
	sim.runSim(ic, tof, &traj);

	double primPos[3];
	sys.getDynamicsModel()->getPrimPos(0, &sys, apseData[0], primPos);

	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(traj.getNodeRefByIx_const(n).getTriggerEvent() == evt.getType()){
			std::vector<double> q = traj.getStateByIx(n);
			double relPos[3] = {q[0] - primPos[0], q[1] - primPos[1], q[2] - primPos[2]};
			double relVel[3] = {q[3], q[4], q[5]};
			double rdot = relPos[0]*relVel[0] + relPos[1]*relVel[1] + relPos[2]*relVel[2];

			BOOST_CHECK_SMALL(rdot, sim.getAbsTol());
		}
	}
}//====================================================

BOOST_AUTO_TEST_SUITE_END()