#define BOOST_TEST_MODULE Arcset

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>

#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_cr3bp_lt.hpp"
#include "AsciiOutput.hpp"
#include "Constraint.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "MultShootEngine.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

bool dummy_predicate(Exception&);
bool dummy_predicate( Exception const& ex ) { return true; }

BOOST_AUTO_TEST_CASE(Concat_CR3BP){
	SysData_cr3bp sys("Saturn", "Titan");
	SysData_cr3bp emSys("earth", "moon");
	Arcset_cr3bp set1(&sys);
	Arcset_cr3bp set2(&sys);
	Arcset_cr3bp set3(&sys);
	Arcset_cr3bp set4(&emSys);

	double state1[] = {1,0,0,0,0,0};
	double state2[] = {2,1,0,0,0,0};
	double state3[] = {3,0,1,0,0,0};
	double state4[] = {4,0,0,1,0,0};

	set1.addNode(Node(state1, 6, 0));
	set1.addNode(Node(state2, 6, 1.1));
	set1.addSeg(Segment(0, 1, 1.1));

	set2.addNode(Node(state3, 6, 0));
	set2.addNode(Node(state4, 6, 2.2));
	set2.addSeg(Segment(0, 1, 2.2));

	set3.addNode(Node(state4, 6, 3.3));
	set3.addNode(Node(state3, 6, 4.4));
	set3.addNode(Node(state2, 6, 5.5));
	set3.addSeg(Segment(0, 1, 1.1));
	set3.addSeg(Segment(1, 2, 1.1));

	set4.addNode(Node(state4, 6, 0));

	Arcset_cr3bp sum1 = set1 + set2;
	BOOST_CHECK_EQUAL(sum1.getStateByIx(0)[0], 1);
	BOOST_CHECK_EQUAL(sum1.getStateByIx(1)[0], 3);
	BOOST_CHECK_EQUAL(sum1.getStateByIx(2)[0], 4);

	Arcset_cr3bp sum2 = set1;
	sum2 += set3;
	BOOST_CHECK_EQUAL(sum2.getStateByIx(0)[0], 1);
	BOOST_CHECK_EQUAL(sum2.getStateByIx(1)[0], 4);
	BOOST_CHECK_EQUAL(sum2.getStateByIx(2)[0], 3);
	BOOST_CHECK_EQUAL(sum2.getStateByIx(3)[0], 2);

	// Sum of different systems
	Arcset_cr3bp sum3 = set1;
	BOOST_CHECK_EXCEPTION(sum3 += set4, Exception, dummy_predicate);	
}//=======================================

BOOST_AUTO_TEST_CASE(CR3BP_NodesAtEvents){
	SysData_cr3bp sys("earth", "moon");
	double emDRO_ic[] = {0.66703088566639, 0, 0, 0, 0.763253816058075, 0};
	double emDRO_T = 5.18136624737627;

	// Second test case: Generate orbit, use createNodesAtEvent and check the functionality, TOF computation, etc.
	Arcset_cr3bp set2(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emDRO_ic, emDRO_T, 2, &set2);
	
	// set2.saveToMat("emDRO_2Nodes.mat");
	// cout << "CR3BP Arcset generated from ICs (saved to emDRO_2Nodes.mat):" << endl;
	BOOST_CHECK(set2.getNumNodes() == 2);
	BOOST_CHECK(set2.getTOFByIx(0) == emDRO_T);

	std::vector<double> xMoonData = {1 - sys.getMu()};
	Event xMoonEvt(Event_tp::YZ_PLANE, 0, true, xMoonData);
	Event xzPlaneEvt(Event_tp::XZ_PLANE, 0, true);
	std::vector<Event> events {xMoonEvt, xzPlaneEvt};

	set2.createNodesAtEvents(0, events);
	set2.putInChronoOrder();
	// set2.printInChrono();
	// set2.saveToMat("data/emDRO_newNodes.mat");
	// cout << "CR3BP createNodesAtEvents (saved to emDRO_newNodes.mat):" << endl;
	BOOST_CHECK(set2.getNumNodes() == 5);
	BOOST_CHECK(std::abs(set2.getStateByIx(1)[0] - xMoonData[0]) < 1e-10);
	BOOST_CHECK(std::abs(set2.getStateByIx(2)[1]) < 1e-10);
	BOOST_CHECK(std::abs(set2.getStateByIx(3)[0] - xMoonData[0]) < 1e-10);
	BOOST_CHECK(std::abs(set2.getTotalTOF() - emDRO_T) < 1e-10);
	// set2.print();
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_NodesAtEvent){
	SysData_bc4bp bcSys("sun", "earth", "moon");
	double qho_ic[] = {-0.86464955943628, -0.523239865136876, -0.0309591111054232, -0.00352683110021282, -0.00217207557203108, 0.00179392516522105};
	double qho_T0 = 100;
	double qho_Period = 360;

	Arcset_bc4bp set3(&bcSys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(qho_ic, qho_T0, qho_Period, 2, &set3);
	
	// cout << "BC4BP Arcset generated from ICs:" << endl;
	BOOST_CHECK(set3.getNumNodes() == 2);
	BOOST_CHECK(set3.getTOFByIx(0) == qho_Period);
	BOOST_CHECK(set3.getEpochByIx(0) == qho_T0);
	BOOST_CHECK(set3.getEpochByIx(1) == qho_T0 + qho_Period);
	// set3.print();

	Event sem_xzPlaneEvt(Event_tp::XZ_PLANE, 0, false);
	set3.createNodesAtEvent(0, sem_xzPlaneEvt);
	set3.putInChronoOrder();
	// set3.printInChrono();
	// set3.saveToMat("semQHO_newNodes.mat");
	// cout << "BC4BP createNodesAtEvent (saved to semQHO_newNodes.mat):" << endl;
	BOOST_CHECK(set3.getNumNodes() == 4);
	BOOST_CHECK(std::abs(set3.getStateByIx(1)[1]) < 1e-12);
	BOOST_CHECK(std::abs(set3.getEpochByIx(1) - qho_T0 - set3.getTOFByIx(0)) < 1e-12);
	BOOST_CHECK(std::abs(set3.getStateByIx(2)[1]) < 1e-12);
	BOOST_CHECK(std::abs(set3.getEpochByIx(2) - qho_T0 - set3.getTOFByIx(0)) < 1e-12 + set3.getTOFByIx(1));
	BOOST_CHECK(std::abs(set3.getTotalTOF() - qho_Period) < 1e-10);
	// set3.print();
}//====================================================

BOOST_AUTO_TEST_CASE(Arcset_MixedTime_Save_Load){
	SysData_cr3bp sys("earth", "moon");
	double state[42] = {0};
	for(unsigned int i = 0; i < 42; i++){
		if(i < 6)
			state[i] = i+1;	// state = [1, 2, 3, 4, 5, 6]
		else{
			if((i+1) % 7 == 0)
				state[i] = 1;	// Identity matrix for STM
		}
	}
	double errTol = 1e-14;

	const DynamicsModel *model = sys.getDynamicsModel();
	EOM_ParamStruct params(&sys, nullptr);
	Arcset_cr3bp set(&sys);

	Node node1(state, 6, 2.2), node2(state, 6, 0), node3(state, 6, 1.1), node4(state, 6, -2.2), node5(state, 6, -1.1);

	model->sim_addNode(node1, state, node1.getEpoch(), &set, &params, Event_tp::SIM_TOF);
	model->sim_addNode(node2, state, node2.getEpoch(), &set, &params, Event_tp::SIM_TOF);
	model->sim_addNode(node3, state, node3.getEpoch(), &set, &params, Event_tp::SIM_TOF);
	model->sim_addNode(node4, state, node4.getEpoch(), &set, &params, Event_tp::SIM_TOF);
	model->sim_addNode(node5, state, node5.getEpoch(), &set, &params, Event_tp::SIM_TOF);

	// set.addNode(Node(state, 6, 2.2));
	// set.addNode(Node(state, 6, 0));
	// set.addNode(Node(state, 6, 1.1));
	// set.addNode(Node(state, 6, -2.2));
	// set.addNode(Node(state, 6, -1.1));

	Segment seg1(4, 3, -1.1), seg2(1, 2, 1.1), seg3(1, 4, -1.1), seg4(2, 0, 1.1);
	seg1.setStateWidth(42);
	for(unsigned int i = 0; i < 2; i++){
		seg1.appendState(state, 42);
		seg2.appendState(state, 42);
		seg3.appendState(state, 42);
		seg4.appendState(state, 42);
	}
	seg2.setStateWidth(42);
	seg3.setStateWidth(42);
	seg4.setStateWidth(42);
	model->sim_addSeg(seg1, state, 0, &set, &params);
	model->sim_addSeg(seg2, state, 0, &set, &params);
	model->sim_addSeg(seg3, state, 0, &set, &params);
	model->sim_addSeg(seg4, state, 0, &set, &params);

	// set.addSeg(Segment(4, 3, -1.1));
	// set.addSeg(Segment(1, 2, 1.1));
	// set.addSeg(Segment(1, 4, -1.1));
	// set.addSeg(Segment(2, 0, 1.1));

	set.saveToMat("data/mixedTimeSet.mat");

	Arcset_cr3bp tempSet(&sys);
	std::vector<ControlLaw*> loadedLaws;
	tempSet.readFromMat("data/mixedTimeSet.mat", loadedLaws);

	BOOST_REQUIRE_EQUAL(set.getNumNodes(), tempSet.getNumNodes());
	for(unsigned int n = 0; n < set.getNumNodes(); n++){
		std::vector<double> state1 = set.getStateByIx(n);
		std::vector<double> state2 = tempSet.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(set.getEpochByIx(n) - tempSet.getEpochByIx(n), errTol);
	}

	BOOST_REQUIRE_EQUAL(set.getNumSegs(), tempSet.getNumSegs());
	for(unsigned int s = 0; s < set.getNumSegs(); s++){
		BOOST_CHECK_EQUAL(set.getSegRefByIx(s).getOrigin(), tempSet.getSegRefByIx(s).getOrigin());
		BOOST_CHECK_EQUAL(set.getSegRefByIx(s).getTerminus(), tempSet.getSegRefByIx(s).getTerminus());

		std::vector<double> states1 = set.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = tempSet.getSegRefByIx(s).getStateVector();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
		}

		std::vector<double> times1 = set.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = tempSet.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(set.getTOFByIx(s) - tempSet.getTOFByIx(s), errTol);

		MatrixXRd stm1 = set.getSTMByIx(s);
		MatrixXRd stm2 = tempSet.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK_EQUAL(set.getCtrlLawByIx(s), tempSet.getCtrlLawByIx(s));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Arcset_Save_Load){
	double errTol = 1e-14;
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};
	SysData_cr3bp emData("earth", "moon");
	Arcset_cr3bp crSet(&emData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, 2.77, 5, &crSet);

	// crSet.print();

	// Add a constraint
	std::vector<double> conData {1, 2, 3, 4, 5, 6};
	Constraint con(Constraint_tp::STATE, 0, conData);
	crSet.addConstraint(con);

	// Save the arcset
	crSet.saveToMat("data/crSet.mat");
	
	// Load the same arcset from the file
	Arcset_cr3bp crTemp(&emData);
	std::vector<ControlLaw*> loadedLaws;
	crTemp.readFromMat("data/crSet.mat", loadedLaws);

	BOOST_CHECK_EQUAL(loadedLaws.size(), 0);

	BOOST_CHECK_EQUAL(crSet.getNumNodes(), crTemp.getNumNodes());
	for(unsigned int n = 0; n < crSet.getNumNodes(); n++){
		std::vector<double> state1 = crSet.getStateByIx(n);
		std::vector<double> state2 = crTemp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(crSet.getEpochByIx(n) - crTemp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(crSet.getNumSegs(), crTemp.getNumSegs());
	for(unsigned int s = 0; s < crSet.getNumSegs(); s++){

		std::vector<double> states1 = crSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = crTemp.getSegRefByIx(s).getStateVector();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
		}

		std::vector<double> times1 = crSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = crTemp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(crSet.getTOFByIx(s) - crTemp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = crSet.getSTMByIx(s);
		MatrixXRd stm2 = crTemp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK_EQUAL(crSet.getCtrlLawByIx(s), crTemp.getCtrlLawByIx(s));
	}

	// Check to make sure that the constraint is saved/loaded properly
	BOOST_CHECK_EQUAL(crTemp.getNumCons(), 1);
	std::vector<Constraint> cons = crTemp.getNodeRefByIx(0).getConstraints();
	BOOST_CHECK_EQUAL(cons.size(), 1);
	BOOST_CHECK(cons[0] == con);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Arcset_SaveCurve_Load){
	double errTol = 1e-14;
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};
	SysData_cr3bp emData("earth", "moon");
	Arcset_cr3bp crSet(&emData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, 2.77, 5, &crSet);

	unsigned int coreDim = emData.getDynamicsModel()->getCoreStateSize();
	// crSet.print();

	crSet.saveToMat("data/crSet_curve.mat", Save_tp::SAVE_CURVE);
	Arcset_cr3bp crTemp(&emData);
	std::vector<ControlLaw*> loadedLaws;
	crTemp.readFromMat("data/crSet_curve.mat", loadedLaws);

	BOOST_CHECK_EQUAL(loadedLaws.size(), 0);

	BOOST_CHECK_EQUAL(crSet.getNumNodes(), crTemp.getNumNodes());
	for(unsigned int n = 0; n < crSet.getNumNodes(); n++){
		std::vector<double> state1 = crSet.getStateByIx(n);
		std::vector<double> state2 = crTemp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(crSet.getEpochByIx(n) - crTemp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(crSet.getNumSegs(), crTemp.getNumSegs());
	for(unsigned int s = 0; s < crSet.getNumSegs(); s++){

		std::vector<double> states1 = crSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = crTemp.getSegRefByIx(s).getStateVector();
		unsigned int stateDim = crTemp.getSegRefByIx(s).getStateWidth();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			if(i % stateDim < coreDim)
				BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
			else{
				// Only core states are saved/loaded; others should be zero
				BOOST_CHECK_SMALL(states2[i], errTol);
			}
		}

		std::vector<double> times1 = crSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = crTemp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(crSet.getTOFByIx(s) - crTemp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = crSet.getSTMByIx(s);
		MatrixXRd stm2 = crTemp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK_EQUAL(crSet.getCtrlLawByIx(s), crTemp.getCtrlLawByIx(s));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_Arcset_Save_Load){
	double errTol = 1e-14;
	SysData_bc4bp semData("sun", "earth", "moon");
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	Arcset_bc4bp bcSet(&semData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic2, 0, 40, 5, &bcSet);	

	bcSet.saveToMat("data/bcSet.mat");
	Arcset_bc4bp bcTemp(&semData);
	std::vector<ControlLaw*> loadedLaws;
	bcTemp.readFromMat("data/bcSet.mat", loadedLaws);
	
	BOOST_CHECK_EQUAL(loadedLaws.size(), 0);

	BOOST_CHECK_EQUAL(bcSet.getNumNodes(), bcTemp.getNumNodes());
	for(unsigned int n = 0; n < bcSet.getNumNodes(); n++){
		std::vector<double> state1 = bcSet.getStateByIx(n);
		std::vector<double> state2 = bcTemp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		std::vector<double> dqdT1 = bcSet.get_dqdTByIx(n);
		std::vector<double> dqdT2 = bcTemp.get_dqdTByIx(n);

		BOOST_CHECK_EQUAL(dqdT1.size(), dqdT2.size());
		for(unsigned int c = 0; c < dqdT1.size(); c++){
			BOOST_CHECK_SMALL(dqdT1[c] - dqdT2[c], errTol);
		}

		BOOST_CHECK_SMALL(bcSet.getEpochByIx(n) - bcTemp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(bcSet.getNumSegs(), bcTemp.getNumSegs());
	for(unsigned int s = 0; s < bcSet.getNumSegs(); s++){

		std::vector<double> states1 = bcSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = bcTemp.getSegRefByIx(s).getStateVector();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
		}

		std::vector<double> times1 = bcSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = bcTemp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(bcSet.getTOFByIx(s) - bcTemp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = bcSet.getSTMByIx(s);
		MatrixXRd stm2 = bcTemp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK_EQUAL(bcSet.getCtrlLawByIx(s), bcTemp.getCtrlLawByIx(s));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Arcset_Save_Load){
	double errTol = 1e-14;
	std::vector<double> ic {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};
	std::vector<double> ctrl0 {1.25, 0.1};

	SysData_cr3bp_lt ltData("earth", "moon", 14);
	std::vector<double> ltParams {0.3, 1500};
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, ltParams);

	Arcset_cr3bp_lt ltSet(&ltData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, ctrl0, 0, 2.8, 5, &ltSet, &control);
	
	// ltSet.print();
	ltSet.saveToMat("data/ltSet.mat");
	Arcset_cr3bp_lt temp(&ltData);
	std::vector<ControlLaw*> loadedLaws {};
	temp.readFromMat("data/ltSet.mat", loadedLaws);

	BOOST_CHECK_EQUAL(loadedLaws.size(), 1);

	BOOST_CHECK_EQUAL(ltSet.getNumNodes(), temp.getNumNodes());
	for(unsigned int n = 0; n < ltSet.getNumNodes(); n++){
		std::vector<double> state1 = ltSet.getStateByIx(n);
		std::vector<double> state2 = temp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(ltSet.getEpochByIx(n) - temp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(ltSet.getNumSegs(), temp.getNumSegs());
	for(unsigned int s = 0; s < ltSet.getNumSegs(); s++){

		std::vector<double> states1 = ltSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = temp.getSegRefByIx(s).getStateVector();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
		}

		std::vector<double> times1 = ltSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = temp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(ltSet.getTOFByIx(s) - temp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = ltSet.getSTMByIx(s);
		MatrixXRd stm2 = temp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK(*(ltSet.getCtrlLawByIx(s)) == *(temp.getCtrlLawByIx(s)));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}

}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Arcset_SaveCurve_Load){
	double errTol = 1e-14;
	std::vector<double> ic {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};
	std::vector<double> ctrl0 {1.25, 0.1};

	SysData_cr3bp_lt ltData("earth", "moon", 14);
	std::vector<double> ltParams {0.3, 1500};
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, ltParams);

	unsigned int coreDim = ltData.getDynamicsModel()->getCoreStateSize();

	Arcset_cr3bp_lt ltSet(&ltData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, ctrl0, 0, 2.8, 5, &ltSet, &control);
	
	// ltSet.print();
	ltSet.saveToMat("data/ltSet_curve.mat", Save_tp::SAVE_CURVE);
	Arcset_cr3bp_lt temp(&ltData);
	std::vector<ControlLaw*> loadedLaws;
	temp.readFromMat("data/ltSet_curve.mat", loadedLaws);

	BOOST_CHECK_EQUAL(loadedLaws.size(), 1);

	BOOST_CHECK_EQUAL(ltSet.getNumNodes(), temp.getNumNodes());
	for(unsigned int n = 0; n < ltSet.getNumNodes(); n++){
		std::vector<double> state1 = ltSet.getStateByIx(n);
		std::vector<double> state2 = temp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(ltSet.getEpochByIx(n) - temp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(ltSet.getNumSegs(), temp.getNumSegs());
	for(unsigned int s = 0; s < ltSet.getNumSegs(); s++){
		unsigned int ctrlDim = temp.getCtrlLawByIx(s)->getNumStates();

		std::vector<double> states1 = ltSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = temp.getSegRefByIx(s).getStateVector();
		unsigned int stateDim = temp.getSegRefByIx(s).getStateWidth();

		BOOST_CHECK_EQUAL(states1.size(), states2.size());
		for(unsigned int i = 0; i < states1.size(); i++){
			if(i % stateDim < coreDim + ctrlDim)
				BOOST_CHECK_SMALL(states1[i] - states2[i], errTol);
			else{
				// Only core states are saved/loaded; others should be zero
				BOOST_CHECK_SMALL(states2[i], errTol);
			}
		}

		std::vector<double> times1 = ltSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> times2 = temp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(times1.size(), times2.size());
		for(unsigned int i = 0; i < times1.size(); i++){
			BOOST_CHECK_SMALL(times1[i] - times2[i], errTol);
		}

		BOOST_CHECK_SMALL(ltSet.getTOFByIx(s) - temp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = ltSet.getSTMByIx(s);
		MatrixXRd stm2 = temp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK(*(ltSet.getCtrlLawByIx(s)) == *(temp.getCtrlLawByIx(s)));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}

}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Arcset_SaveFrame_Load){
	double errTol = 1e-14;
	std::vector<double> ic {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};
	std::vector<double> ctrl0 {1.25, 0.1};

	SysData_cr3bp_lt ltData("earth", "moon", 14);
	std::vector<double> ltParams {0.3, 1500};
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, ltParams);

	unsigned int coreDim = ltData.getDynamicsModel()->getCoreStateSize();

	Arcset_cr3bp_lt ltSet(&ltData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, ctrl0, 0, 2.8, 5, &ltSet, &control);
	
	// ltSet.print();
	ltSet.saveToMat("data/ltSet_frame.mat", Save_tp::SAVE_FRAME);
	Arcset_cr3bp_lt temp(&ltData);
	std::vector<ControlLaw*> loadedLaws;
	temp.readFromMat("data/ltSet_frame.mat", loadedLaws);

	BOOST_CHECK_EQUAL(loadedLaws.size(), 1);

	BOOST_CHECK_EQUAL(ltSet.getNumNodes(), temp.getNumNodes());
	for(unsigned int n = 0; n < ltSet.getNumNodes(); n++){
		std::vector<double> state1 = ltSet.getStateByIx(n);
		std::vector<double> state2 = temp.getStateByIx(n);

		BOOST_CHECK_EQUAL(state1.size(), state2.size());
		for(unsigned int c = 0; c < state1.size(); c++){
			BOOST_CHECK_SMALL(state1[c] - state2[c], errTol);
		}

		BOOST_CHECK_SMALL(ltSet.getEpochByIx(n) - temp.getEpochByIx(n), errTol);
	}

	BOOST_CHECK_EQUAL(ltSet.getNumSegs(), temp.getNumSegs());
	for(unsigned int s = 0; s < ltSet.getNumSegs(); s++){
		unsigned int ctrlDim = temp.getCtrlLawByIx(s)->getNumStates();

		std::vector<double> states1 = ltSet.getSegRefByIx(s).getStateVector();
		std::vector<double> states2 = temp.getSegRefByIx(s).getStateVector();

		BOOST_CHECK_EQUAL(temp.getSegRefByIx(s).getStateWidth(), ltSet.getSegRefByIx(s).getStateWidth());
		unsigned int stateDim = temp.getSegRefByIx(s).getStateWidth();

		BOOST_CHECK_EQUAL(states2.size(), 2*stateDim);

		std::vector<double> ref_q0 = ltSet.getSegRefByIx(s).getStateByRow(0);
		std::vector<double> ref_qf = ltSet.getSegRefByIx(s).getStateByRow(-1);
		std::vector<double> load_q0 = temp.getSegRefByIx(s).getStateByRow(0);
		std::vector<double> load_qf = temp.getSegRefByIx(s).getStateByRow(-1);

		BOOST_CHECK_EQUAL(ref_q0.size(), stateDim);
		BOOST_CHECK_EQUAL(ref_qf.size(), stateDim);
		BOOST_CHECK_EQUAL(load_q0.size(), stateDim);
		BOOST_CHECK_EQUAL(load_qf.size(), stateDim);

		for(unsigned int c = 0; c < stateDim; c++){
			if(c < coreDim + ctrlDim){
				BOOST_CHECK_SMALL(load_q0[c] - load_q0[c], errTol);
				BOOST_CHECK_SMALL(load_qf[c] - load_qf[c], errTol);
			}else{
				// Only core states are saved/loaded; others should be zero
				BOOST_CHECK_SMALL(load_q0[c], errTol);
				BOOST_CHECK_SMALL(load_qf[c], errTol);
			}
		}

		std::vector<double> ref_times = ltSet.getSegRefByIx(s).getTimeVector();
		std::vector<double> load_times = temp.getSegRefByIx(s).getTimeVector();
		BOOST_CHECK_EQUAL(load_times.size(), 2);
		BOOST_CHECK_SMALL(ref_times.front() - load_times[0], errTol);
		BOOST_CHECK_SMALL(ref_times.back() - load_times[1], errTol);

		BOOST_CHECK_SMALL(ltSet.getTOFByIx(s) - temp.getTOFByIx(s), errTol);

		MatrixXRd stm1 = ltSet.getSTMByIx(s);
		MatrixXRd stm2 = temp.getSTMByIx(s);
		BOOST_CHECK_EQUAL(stm1.rows(), stm2.rows());
		BOOST_CHECK_EQUAL(stm1.cols(), stm2.cols());
		for(unsigned int r = 0; r < stm1.rows(); r++){
			for(unsigned int c = 0; c < stm1.cols(); c++){
				BOOST_CHECK_SMALL(stm1(r,c) - stm2(r,c), errTol);
			}
		}

		BOOST_CHECK(*(ltSet.getCtrlLawByIx(s)) == *(temp.getCtrlLawByIx(s)));
	}

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}

}//====================================================




