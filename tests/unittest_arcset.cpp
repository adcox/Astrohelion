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

	set3.addNode(Node(state2, 6, 3.3));
	set3.addNode(Node(state3, 6, 4.4));
	set3.addNode(Node(state4, 6, 5.5));
	set3.addSeg(Segment(0, 1, 1.1));
	set3.addSeg(Segment(1, 2, 1.1));

	set4.addNode(Node(state4, 6, 0));

	Arcset_cr3bp sum1 = set1 + set2;
	BOOST_CHECK(sum1.getStateByIx(0)[0] == 1);
	BOOST_CHECK(sum1.getStateByIx(1)[0] == 2);
	BOOST_CHECK(sum1.getStateByIx(2)[0] == 4);

	Arcset_cr3bp sum2 = set1;
	sum2 += set3;
	BOOST_CHECK(sum2.getStateByIx(0)[0] == 1);
	BOOST_CHECK(sum2.getStateByIx(1)[0] == 2);
	BOOST_CHECK(sum2.getStateByIx(2)[0] == 3);
	BOOST_CHECK(sum2.getStateByIx(3)[0] == 4);

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

BOOST_AUTO_TEST_CASE(CR3BP_Nodeset_Save_Load){
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};
	SysData_cr3bp emData("earth", "moon");
	Arcset_cr3bp crSet(&emData);
	SimEngine sim;
	sim.runSim_manyNodes(ic, 2.77, 5, &crSet);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(&crSet, NULL);

	crSet.saveToMat("data/crSet.mat");
	Arcset_cr3bp crTemp(&emData);
	std::vector<ControlLaw*> loadedLaws;
	crTemp.readFromMat("data/crSet.mat", loadedLaws);

	BOOST_CHECK(crSet.getStateByIx(-1) == crTemp.getStateByIx(-1));
	BOOST_CHECK(crSet.getTOFByIx(-1) == crTemp.getTOFByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 0);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_Nodeset_Save_Load){
	SysData_bc4bp semData("sun", "earth", "moon");
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	Arcset_bc4bp bcSet(&semData);
	SimEngine sim;
	sim.runSim_manyNodes(ic2, 0, 40, 5, &bcSet);

	// Add a constraint
	// double data[] = {82.576, 0, 8.001, NAN, NAN, NAN, NAN};
	// double data[] = {5,5,5,NAN,NAN,NAN,NAN};
	double data[] = {1.5};
	Constraint bcCon1(Constraint_tp::MAX_DELTA_V, 0, data, 1);
	bcSet.addConstraint(bcCon1);

	int nodes[] = {2,3};
	std::vector<int> velCon(nodes, nodes+2);
	bcSet.allowDV_at(velCon);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(&bcSet, nullptr);

	bcSet.saveToMat("data/bcSet.mat");
	Arcset_bc4bp bcTemp(&semData);
	std::vector<ControlLaw*> loadedLaws;
	bcTemp.readFromMat("data/bcSet.mat", loadedLaws);

	BOOST_CHECK(bcSet.getStateByIx(-1) == bcTemp.getStateByIx(-1));
	BOOST_CHECK(bcSet.getTOFByIx(-1) == bcTemp.getTOFByIx(-1));
	BOOST_CHECK(bcSet.getEpochByIx(-1) == bcTemp.getEpochByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 0);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Nodeset_Save_Load){
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};

	MultShootEngine corrector;

	SysData_cr3bp_lt ltData("earth", "moon", 14);
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);

	Arcset_cr3bp_lt ltSet(&ltData);
	SimEngine sim;
	sim.runSim_manyNodes(ic, 2.77, 5, &ltSet, &control);

	corrector.multShoot(&ltSet, nullptr);

	ltSet.saveToMat("data/ltSet.mat");
	Arcset_cr3bp_lt temp(&ltData);
	std::vector<ControlLaw*> loadedLaws;
	temp.readFromMat("data/ltSet.mat", loadedLaws);

	BOOST_CHECK(ltSet.getStateByIx(-1) == temp.getStateByIx(-1));
	BOOST_CHECK(temp.getCtrlLawByIx(0) == ltSet.getCtrlLawByIx(0));
	BOOST_CHECK(*(ltSet.getCtrlLawByIx(0)) == control);
	BOOST_CHECK(ltSet.getTOFByIx(-1) == temp.getTOFByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 1);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}

}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Save_Load){
	SysData_cr3bp emData("earth", "moon");
	SimEngine sim;
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	Arcset_cr3bp crTraj(&emData);
	sim.runSim(ic, T, &crTraj);

	// Query the acceleration so it is computed
	std::vector<double> a = crTraj.getStateDerivByIx(-1);
	crTraj.getJacobiByIx(-1);
	crTraj.saveToMat("data/crTraj.mat");

	Arcset_cr3bp crTemp(&emData);
	std::vector<ControlLaw*> loadedLaws;
	crTemp.readFromMat("data/crTraj.mat", loadedLaws);

	// printf("Testing Save/Read functions on CR3BP Trajectory\n");
	BOOST_CHECK(crTraj.getStateByIx(-1) == crTemp.getStateByIx(-1));
	BOOST_CHECK(crTraj.getStateDerivByIx(-1) == crTemp.getStateDerivByIx(-1));
	BOOST_CHECK(crTraj.getTimeByIx(-1) == crTemp.getTimeByIx(-1));
	BOOST_CHECK(crTraj.getSTMByIx(-1) == crTemp.getSTMByIx(-1));
	BOOST_CHECK(crTraj.getJacobiByIx(-1) == crTemp.getJacobiByIx(-1));
	BOOST_CHECK(crTraj.getTOFByIx(-1) == crTemp.getTOFByIx(-1));
	BOOST_CHECK(crTraj.getCtrlLawByIx(-1) == crTemp.getCtrlLawByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 0);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_Save_Load){
	SysData_bc4bp semData("sun", "earth", "moon");
	SimEngine sim;
	double ic[] = {-0.745230328320519, 7.22625684942683e-04, 7.45549413286038e-05, -7.30710697247992e-06, -0.0148897145134465, -1.23266135281459e-06};
	double T = 313;	// SE L1 Period
	Arcset_bc4bp bcTraj(&semData);
	sim.runSim(ic, T, &bcTraj);

	// Query the acceleration so it is computed
	std::vector<double> a = bcTraj.getStateDerivByIx(-1);
	bcTraj.saveToMat("data/bcTraj.mat");

	Arcset_bc4bp bcTemp(&semData);
	std::vector<ControlLaw*> loadedLaws;
	bcTemp.readFromMat("data/bcTraj.mat", loadedLaws);

	// printf("Testing Save/Read functions on BC4BP Trajectory\n");
	BOOST_CHECK(bcTraj.getStateByIx(-1) == bcTemp.getStateByIx(-1));
	BOOST_CHECK(bcTraj.getStateDerivByIx(-1) == bcTemp.getStateDerivByIx(-1));
	BOOST_CHECK(bcTraj.getTimeByIx(-1) == bcTemp.getTimeByIx(-1));
	BOOST_CHECK(bcTraj.getSTMByIx(-1) == bcTemp.getSTMByIx(-1));
	BOOST_CHECK(bcTraj.get_dqdTByIx(-1) == bcTemp.get_dqdTByIx(-1));
	BOOST_CHECK(bcTraj.getTOFByIx(-1) == bcTemp.getTOFByIx(-1));
	BOOST_CHECK(bcTraj.getCtrlLawByIx(-1) == bcTemp.getCtrlLawByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 0);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Save_Load){
	SysData_cr3bp_lt emData("earth", "moon", 14);
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);

	SimEngine sim;
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	Arcset_cr3bp_lt ltTraj(&emData);
	sim.runSim(ic, T, &ltTraj, &control);

	// Query the acceleration so it is computed
	std::vector<double> a = ltTraj.getStateDerivByIx(-1);
	ltTraj.getJacobiByIx(-1);
	ltTraj.saveToMat("data/lowthrustTraj.mat");

	Arcset_cr3bp_lt ltTemp(&emData);
	std::vector<ControlLaw*> loadedLaws;
	ltTemp.readFromMat("data/lowthrustTraj.mat", loadedLaws);

	// printf("Testing Save/Read functions on CR3BP Trajectory\n");
	BOOST_CHECK(ltTraj.getStateByIx(-1) == ltTemp.getStateByIx(-1));
	BOOST_CHECK(ltTraj.getStateDerivByIx(-1) == ltTemp.getStateDerivByIx(-1));
	BOOST_CHECK(ltTraj.getTimeByIx(-1) == ltTemp.getTimeByIx(-1));
	BOOST_CHECK(ltTraj.getSTMByIx(-1) == ltTemp.getSTMByIx(-1));
	BOOST_CHECK(ltTraj.getJacobiByIx(-1) == ltTemp.getJacobiByIx(-1));
	BOOST_CHECK(ltTraj.getTOFByIx(-1) == ltTemp.getTOFByIx(-1));
	BOOST_CHECK(ltTraj.getCtrlLawByIx(-1) == ltTemp.getCtrlLawByIx(-1));
	BOOST_CHECK(loadedLaws.size() == 1);

	if(loadedLaws.size() > 0){
		for(auto law : loadedLaws)
			delete(law);
	}
}//====================================================



