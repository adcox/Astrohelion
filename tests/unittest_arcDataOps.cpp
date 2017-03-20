#define BOOST_TEST_MODULE ArcDataOps

#include <boost/test/unit_test.hpp>

#include "Exceptions.hpp"
#include "Linkable.hpp"
#include "Node.hpp"
#include "Nodeset.hpp"
#include "Segment.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

#include <iostream>

using namespace astrohelion;

SysData_cr3bp sys("earth", "moon");
double state1[] = {1, 0, 0, 0, 0, 0};
double state2[] = {2, 2, 0, 0, 0, 0};
double state3[] = {3, 0, 3, 0, 0, 0};
double state4[] = {4, 0, 0, 4, 0, 0};
double state5[] = {5, 0, 0, 0, 5, 0};

int ivID = Linkable::INVALID_ID;

//************************************************************
//* Linkable Tests
//************************************************************

BOOST_AUTO_TEST_SUITE(Linkable)

BOOST_AUTO_TEST_CASE(Creation){
	Node n(state1, 6, 10);
	BOOST_CHECK(n.getLink(0) == ivID);
	BOOST_CHECK(n.getLink(1) == ivID);

	// Add one link; second slot should still have INVALID_ID
	n.addLink(3);
	BOOST_CHECK(n.getLink(0) == 3);
	BOOST_CHECK(n.getLink(1) == ivID);

	// Add a second link
	n.addLink(7);
	BOOST_CHECK(n.getLink(0) == 3);
	BOOST_CHECK(n.getLink(1) == 7);	
}//====================================================

BOOST_AUTO_TEST_CASE(func_isLinkedTo){
	Node n(state1, 6, 10);
	n.addLink(3);
	n.addLink(7);

	BOOST_CHECK(n.isLinkedTo(3));
	BOOST_CHECK(n.isLinkedTo(7));
	BOOST_CHECK(!n.isLinkedTo(4));
}//====================================================

BOOST_AUTO_TEST_CASE(func_clearLinks){
	Node n(state1, 6, 10);
	n.addLink(3);
	n.addLink(7);
	n.clearLinks();

	BOOST_CHECK(n.getLink(0) == ivID);
	BOOST_CHECK(n.getLink(1) == ivID);
}//====================================================

BOOST_AUTO_TEST_CASE(duplicateLinks){
	Node n(state1, 6, 10);
	n.addLink(3);
	BOOST_CHECK_THROW(n.addLink(3), Exception);
}//====================================================

BOOST_AUTO_TEST_CASE(func_removeLink){
	// Test to make sure links are removed correctly
	Node n(state1, 6, 10);
	n.addLink(3);
	n.addLink(7);
	n.removeLink(3);

	BOOST_CHECK(n.getLink(0) == ivID);
	BOOST_CHECK(n.getLink(1) == 7);

	n.removeLink(3);
	BOOST_CHECK(n.getLink(0) == ivID);
	BOOST_CHECK(n.getLink(1) == 7);
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* Arcset Creation
//************************************************************
BOOST_AUTO_TEST_SUITE(ArcsetTests)

BOOST_AUTO_TEST_CASE(createArcset){
	Nodeset set(&sys);

	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	BOOST_REQUIRE(n1ID != ivID);
	BOOST_REQUIRE(n2ID != ivID);

	Segment s(n1ID, n2ID, 15);
	Segment s_bad1(n1ID, 999, 12);
	Segment s_bad2(n2ID, 987, 12);

	int sID = set.addSeg(s);

	BOOST_REQUIRE(sID != ivID);

	Node node1 = set.getNode(n1ID);
	Node node2 = set.getNode(n2ID);

	// Basic node initialization
	BOOST_CHECK(n1.getEpoch() == 10);
	BOOST_CHECK(n1.getLink(0) == ivID);
	BOOST_CHECK(n1.getLink(1) == ivID);

	// Basic segment initialization
	BOOST_CHECK(s.getTOF() == 15);
	BOOST_CHECK(s.getLink(0) == n1ID);
	BOOST_CHECK(s.getLink(1) == n2ID);

	// Nodes linked correctly
	BOOST_CHECK(node1.getLink(0) == sID);
	BOOST_CHECK(node1.getLink(1) == ivID);
	BOOST_CHECK(node2.getLink(0) == sID);
	BOOST_CHECK(node2.getLink(1) == ivID);

	// Test bad segments
	BOOST_CHECK_THROW(set.addSeg(s_bad1), Exception);
	BOOST_CHECK_THROW(set.addSeg(s_bad2), Exception);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteSeg){
	Nodeset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);

	int sID = set.addSeg(s);
	set.deleteSeg(sID);

	BOOST_CHECK(set.getNumSegs() == 0);
	BOOST_CHECK(set.getNode(n1ID).getLink(0) == ivID);
	BOOST_CHECK(set.getNode(n2ID).getLink(0) == ivID);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteFirstNode){
	Nodeset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	set.addSeg(s);

	set.deleteNode(0);	// Should work
}//====================================================

BOOST_AUTO_TEST_CASE(deleteLastNode){
	Nodeset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	set.addSeg(s);

	set.deleteNode(set.getNumNodes()-1); 	// Should work
}//====================================================

BOOST_AUTO_TEST_CASE(deleteMiddleNode_LinearForwardTime){
	Nodeset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, 1.1));
	set.addNode(Node(state3, 6, 2.2));
	set.addNode(Node(state4, 6, 3.3));
	set.addNode(Node(state5, 6, 4.4));
	set.addSeg(Segment(0, 1, 1.1));
	set.addSeg(Segment(1, 2, 1.1));
	set.addSeg(Segment(2, 3, 1.1));
	int sID_last = set.addSeg(Segment(3, 4, 1.1));

	// Delete middle node (LINEAR, FORWARD TIME)
	set.deleteNode(1);
	Segment seg = set.getSeg(sID_last+1);	// should retrieve the new segment
	BOOST_CHECK(seg.getOrigin() == 0);
	BOOST_CHECK(seg.getTerminus() == 2);
	BOOST_CHECK(std::abs(seg.getTOF() - 2.2) < 1e-4);

	// Checking updated nodeIDMap
	BOOST_CHECK(set.getNodeIx(0) == 0);
	BOOST_CHECK(set.getNodeIx(2) == 1);
	BOOST_CHECK(set.getNodeIx(3) == 2);
	BOOST_CHECK(set.getNodeIx(4) == 3);

	// Checking deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);

	// Checking updated segIDMap
	BOOST_CHECK(set.getSegIx(2) == 0);
	BOOST_CHECK(set.getSegIx(3) == 1);

	// Checking deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(1), Exception);

	// Delete another middle node
	set.deleteNode(2);
	Segment seg2 = set.getSeg(sID_last+2);
	BOOST_CHECK(seg2.getOrigin() == 0);
	BOOST_CHECK(seg2.getTerminus() == 3);
	BOOST_CHECK(std::abs(seg2.getTOF() - 3.3) < 1e-4);

	// Checking updated nodeIDMap
	BOOST_CHECK(set.getNodeIx(0) == 0);
	BOOST_CHECK(set.getNodeIx(3) == 1);
	BOOST_CHECK(set.getNodeIx(4) == 2);

	// Checking deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);
	BOOST_CHECK_THROW(set.getNodeIx(2), Exception);

	// Checking updated segIDMap
	BOOST_CHECK(set.getSegIx(3) == 0);
	BOOST_CHECK(set.getSegIx(5) == 1);

	// Checking deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(1), Exception);
	BOOST_CHECK_THROW(set.getSegIx(2), Exception);
	BOOST_CHECK_THROW(set.getSegIx(4), Exception);
}//====================================================

// Next: tryDeleteMiddleNode_revTime()

BOOST_AUTO_TEST_CASE(deleteMiddleNode_revTime){
	Nodeset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, -1.1));
	set.addNode(Node(state3, 6, -2.2));
	set.addNode(Node(state4, 6, -3.3));
	set.addNode(Node(state5, 6, -4.4));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	set.addSeg(Segment(2, 3, -1.1));
	int sID_last = set.addSeg(Segment(3, 4, -1.1));

	// Delete the middle node
	set.deleteNode(1);
	Segment seg = set.getSeg(sID_last+1);	// should retrieve the new segment
	BOOST_CHECK(seg.getOrigin() == 0);
	BOOST_CHECK(seg.getTerminus() == 2);
	BOOST_CHECK(std::abs(seg.getTOF() + 2.2) < 1e-4);

	// Check updated nodeIDMap
	BOOST_CHECK(set.getNodeIx(0) == 0);
	BOOST_CHECK(set.getNodeIx(2) == 1);
	BOOST_CHECK(set.getNodeIx(3) == 2);
	BOOST_CHECK(set.getNodeIx(4) == 3);

	// Check deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);
}//====================================================

BOOST_AUTO_TEST_SUITE_END()