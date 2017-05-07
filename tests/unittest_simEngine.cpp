#define BOOST_TEST_MODULE SimEngineOps

#include <boost/test/unit_test.hpp>

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

//************************************************************
//* Data Storage Tests
//************************************************************
BOOST_AUTO_TEST_SUITE(DataStorage)

BOOST_AUTO_TEST_CASE(CR3BP_Propagation){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo
	double tof = 4*PI;
	
	SimEngine engine;
	// engine.setVerbosity(Verbosity_tp::DEBUG);
	Arcset_cr3bp traj(&sys);
	engine.runSim(ic, tof, &traj);

	// Make sure the correct number of nodes and segments have been created
	BOOST_CHECK(traj.getNumNodes() == 2);
	BOOST_CHECK(traj.getNumSegs() == 1);

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}

	// Make sure the times match
	BOOST_CHECK(std::abs(traj.getTimeByIx(0) - traj.getSegByIx(0).getTimeByIx(0)) < engine.getAbsTol());
	BOOST_CHECK(std::abs(traj.getTimeByIx(-1) - traj.getSegByIx(-1).getTimeByIx(-1)) < engine.getAbsTol());
	BOOST_CHECK(std::abs(traj.getTotalTOF() - tof) < engine.getAbsTol());

	// Make sure the states match
	std::vector<double> initialNodeState = traj.getStateByIx(0);
	std::vector<double> finalNodeState = traj.getStateByIx(-1);
	std::vector<double> initialSegState = traj.getSegByIx(0).getStateByRow(0, 42);	// 42 (6 + 36) state variables per row
	std::vector<double> finalSegState = traj.getSegByIx(-1).getStateByRow(-1, 42);

	for(unsigned int i = 0; i < 6; i++){
		BOOST_CHECK(std::abs(initialNodeState[i] - initialSegState[i]) < engine.getAbsTol());
		BOOST_CHECK(std::abs(finalNodeState[i] - finalSegState[i]) < engine.getAbsTol());
	}
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* Event Function Tests
//************************************************************
BOOST_AUTO_TEST_SUITE(Events)
/**
 *  \brief Check that an event stops integration in a CR3BP simulation
 */
BOOST_AUTO_TEST_CASE(CR3BP_Event_Stop){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	SimEngine engine;
	// engine.setVerbosity(Verbosity_tp::DEBUG);
	Arcset_cr3bp traj(&sys);
	engine.addEvent(Event(Event_tp::XZ_PLANE, 0, true));
	// engine.setRevTime(true);
	engine.runSim(ic, 4*PI, &traj);

	std::vector<double> qf = traj.getStateByIx(-1);
	BOOST_CHECK(std::abs(qf[1]) < engine.getAbsTol());

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

/**
 *  \brief Check that an event does not stop integration
 *  in a CR3BP simulation and that the event is correctly
 *  located in the simulation
 */
BOOST_AUTO_TEST_CASE(CR3BP_Event_NoStop){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	SimEngine engine;
	// engine.setVerbosity(Verbosity_tp::DEBUG);
	Arcset_cr3bp traj(&sys);
	Event planeCross(Event_tp::XZ_PLANE, 0, false);
	engine.addEvent(planeCross);
	double tof = 4*PI;
	engine.runSim(ic, tof, &traj);

	BOOST_CHECK(std::abs(traj.getTimeByIx(-1) - tof) < engine.getAbsTol());
	BOOST_CHECK(std::abs(traj.getTotalTOF() - tof) < engine.getAbsTol());
	
	bool foundEvent = false;
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(traj.getNodeByIx(n).getTriggerEvent() == planeCross.getType()){
			foundEvent = true;
			std::vector<double> q = traj.getStateByIx(n);
			BOOST_CHECK(std::abs(q[1]) < engine.getAbsTol());
		}
	}
	BOOST_CHECK(foundEvent);

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Event_ManyRevs){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.912877132059092, 0, 0, 0, -0.367515774495231, 0};

	SimEngine engine;
	Arcset_cr3bp traj(&sys);
	Event planeCross(Event_tp::XZ_PLANE, 0, true);
	unsigned int stopCount = 4;
	planeCross.setStopCount(stopCount);
	engine.setMakeDefaultEvents(false);
	engine.addEvent(planeCross);
	engine.runSim(ic, 10*PI, &traj);

	// Make sure it ended on the event
	std::vector<double> qf = traj.getStateByIx(-1);
	BOOST_CHECK(std::abs(qf[1]) < engine.getAbsTol());

	// Make sure each event is valid
	bool foundEvent = false;
	unsigned int eventCount = 0;
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(traj.getNodeByIx(n).getTriggerEvent() == planeCross.getType()){
			foundEvent = true;
			eventCount++;
			std::vector<double> q = traj.getStateByIx(n);
			BOOST_CHECK(std::abs(q[1]) < engine.getAbsTol());
		}
	}
	BOOST_CHECK(foundEvent);
	BOOST_CHECK(eventCount == stopCount);

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Event_InALoop){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.912877132059092, 0, 0, 0, -0.367515774495231, 0};

	SimEngine engine;
	Event planeCross(Event_tp::XZ_PLANE, 0, true);
	unsigned int stopCount = 4;
	planeCross.setStopCount(stopCount);
	engine.setMakeDefaultEvents(false);
	engine.addEvent(planeCross);

	for(int it = 0; it < 5; it++){
		Arcset_cr3bp traj(&sys);
		engine.runSim(ic, 10*PI, &traj);

		// Make sure it ended on the event
		std::vector<double> qf = traj.getStateByIx(-1);
		BOOST_CHECK(std::abs(qf[1]) < engine.getAbsTol());

		// Make sure each event is valid
		bool foundEvent = false;
		unsigned int eventCount = 0;
		for(unsigned int n = 0; n < traj.getNumNodes(); n++){
			if(traj.getNodeByIx(n).getTriggerEvent() == planeCross.getType()){
				foundEvent = true;
				eventCount++;
				std::vector<double> q = traj.getStateByIx(n);
				BOOST_CHECK(std::abs(q[1]) < engine.getAbsTol());
			}
		}
		BOOST_CHECK(foundEvent);
		BOOST_CHECK(eventCount == stopCount);

		// Make sure linking has occurred correctly - Assumes all nodes
		// and segments are added in order
		for(unsigned int s = 0; s < traj.getNumSegs(); s++){
			BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
			BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
		}
		for(unsigned int n = 0; n < traj.getNumNodes(); n++){
			if(n == 0 || n == traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
			else{
				if(n < traj.getNumNodes()-1)
					BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

				if(n > 0)
					BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
			}
		}
	}
}//====================================================

/**
 *  \brief Check that an event stops integration in a BC4BP simulation
 */
BOOST_AUTO_TEST_CASE(BC4BP_Event_Stop){
	SysData_bc4bp sys("sun", "earth", "moon");
	double ic[] = {-0.753164160347879, -0.286207797697346, -0.094854106846973, -0.002178801479866, -0.010534646860849, 0.000667421933754};
	double t0 = 0;

	SimEngine engine;
	Arcset_bc4bp traj(&sys);
	engine.addEvent(Event(Event_tp::XY_PLANE, 0, true));
	engine.runSim(ic, t0, 120*PI, &traj);
	
	std::vector<double> qf = traj.getStateByIx(-1);
	BOOST_CHECK(std::abs(qf[2]) < engine.getAbsTol());

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Event_Stop){
	SysData_cr3bp_lt sys("earth", "moon", 12e-3, 1500, 14);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	SimEngine engine;
	Arcset_cr3bp_lt traj(&sys);
	engine.setCtrlLaw(ControlLaw_cr3bp_lt::CONST_C_2D_RIGHT);
	engine.addEvent(Event(Event_tp::XZ_PLANE, 0, true));
	engine.runSim(ic, 0, T, &traj);

	std::vector<double> qf = traj.getStateByIx(-1);
	BOOST_CHECK(std::abs(qf[1]) < engine.getAbsTol());

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

/**
 *  \brief Check that an event does not stop integration
 *  in a CR3BP simulation and that the event is correctly
 *  located in the simulation
 */
BOOST_AUTO_TEST_CASE(CR3BP_LT_Event_NoStop){
	SysData_cr3bp_lt sys("earth", "moon", 12e-3, 1500, 14);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1

	SimEngine engine;
	Arcset_cr3bp traj(&sys);
	Event planeCross(Event_tp::XZ_PLANE, 0, false);
	engine.addEvent(planeCross);
	engine.setCtrlLaw(ControlLaw_cr3bp_lt::CONST_C_2D_LEFT);

	double tof = 4*PI;
	engine.runSim(ic, tof, &traj);

	BOOST_CHECK(std::abs(traj.getTimeByIx(-1) - tof) < engine.getAbsTol());

	// Make sure each event is valid
	bool foundEvent = false;
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(traj.getNodeByIx(n).getTriggerEvent() == planeCross.getType()){
			foundEvent = true;
			std::vector<double> q = traj.getStateByIx(n);
			BOOST_CHECK(std::abs(q[1]) < engine.getAbsTol());
		}
	}
	BOOST_CHECK(foundEvent);

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Event_ManyRevs){
	SysData_cr3bp_lt sys("earth", "moon", 1e-3, 1500, 14);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1

	SimEngine engine;
	Arcset_cr3bp traj(&sys);
	Event planeCross(Event_tp::XZ_PLANE, 0, true);
	unsigned int stopCount = 4;
	planeCross.setStopCount(stopCount);
	engine.setMakeDefaultEvents(false);
	engine.addEvent(planeCross);
	engine.setCtrlLaw(ControlLaw_cr3bp_lt::CONST_C_2D_LEFT);
	engine.runSim(ic, 10*PI, &traj);

	// Make sure it ended on the event
	std::vector<double> qf = traj.getStateByIx(-1);
	BOOST_CHECK(std::abs(qf[1]) < engine.getAbsTol());

	// Make sure each event is valid
	bool foundEvent = false;
	unsigned int eventCount = 0;
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(traj.getNodeByIx(n).getTriggerEvent() == planeCross.getType()){
			foundEvent = true;
			eventCount++;
			std::vector<double> q = traj.getStateByIx(n);
			BOOST_CHECK(std::abs(q[1]) < engine.getAbsTol());
		}
	}
	BOOST_CHECK(foundEvent);
	BOOST_CHECK(eventCount == stopCount);

	// Make sure linking has occurred correctly - Assumes all nodes
	// and segments are added in order
	for(unsigned int s = 0; s < traj.getNumSegs(); s++){
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s));
		BOOST_CHECK(traj.getSeg(s).isLinkedTo(s+1));
	}
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		if(n == 0 || n == traj.getNumNodes()-1)
			BOOST_CHECK(traj.getNode(n).isLinkedTo(Linkable::INVALID_ID));
		else{
			if(n < traj.getNumNodes()-1)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n));

			if(n > 0)
				BOOST_CHECK(traj.getNode(n).isLinkedTo(n-1));
		}
	}
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* Creating many nodes
//************************************************************
BOOST_AUTO_TEST_SUITE(ManyNodes)

BOOST_AUTO_TEST_CASE(CR3BP_Forward){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	SimEngine engine;
	// engine.setVerbosity(Verbosity_tp::DEBUG);
	Arcset_cr3bp traj(&sys);
	engine.runSim_manyNodes(ic, 4.0, 5, &traj);

	// Make sure the correct number of nodes and segments are created
	BOOST_CHECK(traj.getNumNodes() == 5);
	BOOST_CHECK(traj.getNumSegs() == 4);

	// Make sure the nodes are equally spaced in time (should be at integer time values here)
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		BOOST_CHECK(std::abs(traj.getEpochByIx(n) - static_cast<double>(n)) < engine.getAbsTol());
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Reverse){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	SimEngine engine;
	// engine.setVerbosity(Verbosity_tp::DEBUG);
	Arcset_cr3bp traj(&sys);
	engine.setRevTime(true);
	engine.runSim_manyNodes(ic, 4.0, 5, &traj);

	// Make sure the correct number of nodes and segments are created
	BOOST_CHECK(traj.getNumNodes() == 5);
	BOOST_CHECK(traj.getNumSegs() == 4);

	// Make sure the nodes are equally spaced in time (should be at integer time values here)
	for(unsigned int n = 0; n < traj.getNumNodes(); n++){
		BOOST_CHECK(std::abs(traj.getEpochByIx(n) + static_cast<double>(n)) < engine.getAbsTol());
	}
}//==================================================== 

BOOST_AUTO_TEST_SUITE_END()


//************************************************************
//* Simulation Edge Cases
//************************************************************
BOOST_AUTO_TEST_SUITE(EdgeCases)

/**
 *  \brief Check to make sure that data is recorded to a trajectory
 *  object even if the GSL integrator crashes mid-run
 */
BOOST_AUTO_TEST_CASE(DataPreserved){
	SysData_cr3bp sys("sun", "earth");

	double ic[] = {1.0065, 0, 0, 3.19189119579733e-16, 0.0158375372644023, 0};

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	
	// GSL should crash
	Arcset_cr3bp traj(&sys);
	BOOST_CHECK_THROW(sim.runSim(ic, 6*PI*400, &traj), Exception);

	// But some data should still be written to traj
	BOOST_CHECK(traj.getNumNodes() > 1);
	BOOST_CHECK(traj.getNumSegs() > 0);
}//====================================================

BOOST_AUTO_TEST_CASE(DataPreserved_fixedStep){
	SysData_cr3bp sys("sun", "earth");

	double ic[] = {1.0065, 0, 0, 3.19189119579733e-16, 0.0158375372644023, 0};

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.setVarStepSize(false);
	sim.setNumSteps(2);

	// GSL should crash
	Arcset_cr3bp traj(&sys);
	BOOST_CHECK_THROW(sim.runSim(ic, 6*PI*400, &traj), Exception);

	// But some data should still be written to traj
	BOOST_CHECK(traj.getNumNodes() > 1);
	BOOST_CHECK(traj.getNumSegs() > 0);
}//====================================================

/**
 *  \brief Check to make sure the max computation time constraint
 *  functions properly
 */
BOOST_AUTO_TEST_CASE(Timeout){
	SysData_cr3bp sys("earth", "moon");
	double mu = sys.getMu();
	double ic[] = {1-mu,0,0,0,0,0};

	SimEngine sim;
	Arcset_cr3bp traj(&sys);
	sim.setMaxCompTime(3);

	int t0 = time(nullptr);
	sim.runSim(ic, 6*PI*400, &traj);
	int tf = time(nullptr);
	
	// Sim engine should exit after 3 seconds: our timer must be greater than 2 and less than 5
	// tf - t0 is usually 4 (it rounds to nearest second)
	BOOST_CHECK(tf - t0 > 2);
	BOOST_CHECK(tf - t0 < 5);
}//====================================================

BOOST_AUTO_TEST_SUITE_END()