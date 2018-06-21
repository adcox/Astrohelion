#define BOOST_TEST_MODULE ArcDataOps

#include <boost/test/unit_test.hpp>

#include "Exceptions.hpp"
#include "Linkable.hpp"
#include "Node.hpp"
#include "Arcset.hpp"
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

Arcset forwardSet(&sys);
Arcset revSet(&sys);

bool pieceVecsAreEqual(std::vector<ArcPiece> v1, std::vector<ArcPiece> v2){
	if(v1.size() != v2.size())
		return false;

	for(size_t i = 0; i < v1.size(); i++){
		if(v1[i] != v2[i])
			return false;
	}

	return true;
}//====================================================

void initForwardSet(){
	forwardSet = Arcset(&sys);
	forwardSet.addNode(Node(state1, 6, 0));
	forwardSet.addNode(Node(state2, 6, 1.1));
	forwardSet.addNode(Node(state3, 6, 3.3));

	Segment seg1(0, 1, 1.1);
	std::vector<double> t1 {0, 1.1};
	std::vector<double> q1(state1+0, state1+6);
	q1.insert(q1.end(), state2+0, state2+6);
	seg1.setTimeVector(t1);
	seg1.setStateVector(q1);

	Segment seg2(1, 2, 2.2);
	std::vector<double> t2 {1.1, 3.3};
	std::vector<double> q2(state2+0, state2+6);
	q2.insert(q2.end(), state3+0, state3+6);
	seg2.setTimeVector(t2);
	seg2.setStateVector(q2);

	forwardSet.addSeg(seg1);
	forwardSet.addSeg(seg2);
}//====================================================

void initRevSet(){
	revSet = Arcset(&sys);
	revSet.addNode(Node(state1, 6, 0));
	revSet.addNode(Node(state2, 6, -1.1));
	revSet.addNode(Node(state3, 6, -3.3));

	Segment seg1(0, 1, -1.1);
	std::vector<double> t1 {0, -1.1};
	std::vector<double> q1(state1+0, state1+6);
	q1.insert(q1.end(), state2+0, state2+6);
	seg1.setTimeVector(t1);
	seg1.setStateVector(q1);

	Segment seg2(1, 2, -2.2);
	std::vector<double> t2 {-1.1, -3.3};
	std::vector<double> q2(state2+0, state2+6);
	q2.insert(q2.end(), state3+0, state3+6);
	seg2.setTimeVector(t2);
	seg2.setStateVector(q2);

	revSet.addSeg(seg1);
	revSet.addSeg(seg2);
}//====================================================

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
BOOST_AUTO_TEST_SUITE(Arcset_Basics)

BOOST_AUTO_TEST_CASE(createArcset){
	Arcset set(&sys);

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

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* Arcset - Deleting Nodes and Segments
//************************************************************
BOOST_AUTO_TEST_SUITE(Arcset_DeleteNodesAndSegs)

BOOST_AUTO_TEST_CASE(deleteSeg){
	Arcset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);

	int sID = set.addSeg(s);
	set.deleteSeg(sID);

	BOOST_CHECK_EQUAL(set.getNumSegs(), 0);
	BOOST_CHECK_EQUAL(set.getNode(n1ID).getLink(0), ivID);
	BOOST_CHECK_EQUAL(set.getNode(n2ID).getLink(0), ivID);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteFirstNode){
	Arcset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	int s1ID = set.addSeg(s);

	BOOST_CHECK_NO_THROW(set.deleteNode(n1ID));

	// Check to make sure segments linking to this node were updated
	BOOST_CHECK_EQUAL(set.getSeg(s1ID).getOrigin(), ivID);
	BOOST_CHECK_EQUAL(set.getSeg(s1ID).getTerminus(), n2ID);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteLastNode){
	Arcset set(&sys);
	Node n1(state1, 6, 10);
	Node n2(state2, 6, 25);
	Node n3(state3, 6, 40);
	Node n4(state4, 6, 55);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);
	int n3ID = set.addNode(n3);
	int n4ID = set.addNode(n4);

	set.addSeg(Segment(n1ID, n2ID, 15));
	set.addSeg(Segment(n2ID, n3ID, 15));
	int s3ID = set.addSeg(Segment(n3ID, n4ID, 15));

	BOOST_CHECK_NO_THROW(set.deleteNode(n4ID));

	// Check to make sure segments linking to this node were updated
	BOOST_CHECK_EQUAL(set.getSeg(s3ID).getOrigin(), n3ID);
	BOOST_CHECK_EQUAL(set.getSeg(s3ID).getTerminus(), ivID);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteMiddleNode_LinearForwardTime){
	Arcset set(&sys);
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
	BOOST_CHECK_EQUAL(seg.getOrigin(), 0);
	BOOST_CHECK_EQUAL(seg.getTerminus(), 2);
	BOOST_CHECK_CLOSE(seg.getTOF(), 2.2, 1e-4);

	// Checking updated nodeIDMap
	BOOST_CHECK_EQUAL(set.getNodeIx(0), 0);
	BOOST_CHECK_EQUAL(set.getNodeIx(2), 1);
	BOOST_CHECK_EQUAL(set.getNodeIx(3), 2);
	BOOST_CHECK_EQUAL(set.getNodeIx(4), 3);

	// Checking deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);

	// Checking updated segIDMap
	BOOST_CHECK_EQUAL(set.getSegIx(2), 0);
	BOOST_CHECK_EQUAL(set.getSegIx(3), 1);

	// Checking deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(1), Exception);

	// Delete another middle node
	set.deleteNode(2);
	Segment seg2 = set.getSeg(sID_last+2);
	BOOST_CHECK_EQUAL(seg2.getOrigin(), 0);
	BOOST_CHECK_EQUAL(seg2.getTerminus(), 3);
	BOOST_CHECK_CLOSE(seg2.getTOF(), 3.3, 1e-4);

	// Checking updated nodeIDMap
	BOOST_CHECK_EQUAL(set.getNodeIx(0), 0);
	BOOST_CHECK_EQUAL(set.getNodeIx(3), 1);
	BOOST_CHECK_EQUAL(set.getNodeIx(4), 2);

	// Checking deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);
	BOOST_CHECK_THROW(set.getNodeIx(2), Exception);

	// Checking updated segIDMap
	BOOST_CHECK_EQUAL(set.getSegIx(3), 0);
	BOOST_CHECK_EQUAL(set.getSegIx(5), 1);

	// Checking deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(1), Exception);
	BOOST_CHECK_THROW(set.getSegIx(2), Exception);
	BOOST_CHECK_THROW(set.getSegIx(4), Exception);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteMiddleNode_revTime){
	Arcset set(&sys);
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
	BOOST_CHECK_EQUAL(seg.getOrigin(), 0);
	BOOST_CHECK_EQUAL(seg.getTerminus(), 2);
	BOOST_CHECK_CLOSE(seg.getTOF(), -2.2, 1e-4);

	// Check updated nodeIDMap
	BOOST_CHECK_EQUAL(set.getNodeIx(0), 0);
	BOOST_CHECK_EQUAL(set.getNodeIx(2), 1);
	BOOST_CHECK_EQUAL(set.getNodeIx(3), 2);
	BOOST_CHECK_EQUAL(set.getNodeIx(4), 3);

	// Check deleted node
	BOOST_CHECK_THROW(set.getNodeIx(1), Exception);
}//====================================================

BOOST_AUTO_TEST_CASE(deleteMiddleNode_doubleSource1){
	Arcset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, -1.1));
	set.addNode(Node(state3, 6, -2.2));
	set.addNode(Node(state4, 6, 1.1));
	set.addNode(Node(state5, 6, 2.2));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	set.addSeg(Segment(0, 3, 1.1));
	int sID_last = set.addSeg(Segment(3, 4, 1.1));

	// Deleting the node has the correct behavior
	set.deleteNode(0);
	Segment seg = set.getSeg(sID_last+1);
	BOOST_CHECK_EQUAL(seg.getOrigin(), 1);
	BOOST_CHECK_EQUAL(seg.getTerminus(), 3);
	BOOST_CHECK_CLOSE(seg.getTOF(), 2.2, 1e-4);

	// Correct, updated nodeIDMap
	BOOST_CHECK_EQUAL(set.getNodeIx(1), 0);
	BOOST_CHECK_EQUAL(set.getNodeIx(2), 1);
	BOOST_CHECK_EQUAL(set.getNodeIx(3), 2);
	BOOST_CHECK_EQUAL(set.getNodeIx(4), 3);

	// Make sure the deleted node is, in fact, deleted
	BOOST_CHECK_THROW(set.getNodeIx(0), Exception);

	// Correct, updated segIDMap
	BOOST_CHECK_EQUAL(set.getSegIx(1), 0);
	BOOST_CHECK_EQUAL(set.getSegIx(3), 1);
	BOOST_CHECK_EQUAL(set.getSegIx(4), 2);

	// Check deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(2), Exception);
}//==================================================


/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the forward time segment has no terminal point
 */
BOOST_AUTO_TEST_CASE(deleteMiddleNode_doubleSource2){
	Arcset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, -1.1));
	set.addNode(Node(state3, 6, -2.2));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	int sID_last = set.addSeg(Segment(0, ivID, 1.1));

	set.deleteNode(0);
	Segment seg = set.getSeg(sID_last+1);

	BOOST_CHECK_EQUAL(seg.getOrigin(), 1);
	BOOST_CHECK_EQUAL(seg.getTerminus(), ivID);
	BOOST_CHECK_CLOSE(seg.getTOF(), 2.2, 1e-4);
}//==================================================

/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the reverse time segment has no terminal point
 */
BOOST_AUTO_TEST_CASE(deleteMiddleNode_doubleSource3){
	Arcset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, 1.1));
	set.addNode(Node(state3, 6, 2.2));
	set.addSeg(Segment(0, ivID, -1.1));
	set.addSeg(Segment(0, 1, 1.1));
	int sID_last = set.addSeg(Segment(1, 2, 1.1));

	set.deleteNode(0);
	Segment seg = set.getSeg(sID_last+1);

	BOOST_CHECK_EQUAL(seg.getOrigin(), 1);
	BOOST_CHECK_EQUAL(seg.getTerminus(), ivID);
	BOOST_CHECK_CLOSE(seg.getTOF(), -2.2, 1e-4);
}//==================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* getChronoOrder() Function
//************************************************************
BOOST_AUTO_TEST_SUITE(Arcset_getChronoOrder)

BOOST_AUTO_TEST_CASE(ForwardTime){
	Arcset set1(&sys);
	set1.addNode(Node(state1, 6, 0));
	set1.addNode(Node(state2, 6, 1.1));
	set1.addNode(Node(state3, 6, 2.2));
	set1.addSeg(Segment(0, 1, 1.1));
	set1.addSeg(Segment(1, 2, 1.1));

	std::vector<ArcPiece> set1_ans;
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	std::vector<ArcPiece> set1_pieces = set1.getChronoOrder();
	// set1.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set1_pieces, set1_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(ReverseTime){
	Arcset set2(&sys);
	set2.addNode(Node(state1, 6, 0));
	set2.addNode(Node(state2, 6, -1.1));
	set2.addNode(Node(state3, 6, -2.2));
	set2.addSeg(Segment(0, 1, -1.1));
	set2.addSeg(Segment(1, 2, -1.1));

	std::vector<ArcPiece> set2_ans;
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));

	std::vector<ArcPiece> set2_pieces = set2.getChronoOrder();
	// set2.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set2_pieces, set2_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledForwardTime){
	// Shuffled forward time
	Arcset set3(&sys);
	set3.addNode(Node(state2, 6, 1.1));
	set3.addNode(Node(state1, 6, 0));
	set3.addNode(Node(state3, 6, 2.2));
	set3.addSeg(Segment(1, 0, 1.1));
	set3.addSeg(Segment(0, 2, 1.1));
	
	std::vector<ArcPiece> set3_ans;
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	std::vector<ArcPiece> set3_pieces = set3.getChronoOrder();
	// set3.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set3_pieces, set3_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledReverseTime){
	Arcset set4(&sys);
	set4.addNode(Node(state2, 6, -1.1));
	set4.addNode(Node(state1, 6, 0));
	set4.addNode(Node(state3, 6, -2.2));
	set4.addSeg(Segment(1, 0, -1.1));
	set4.addSeg(Segment(0, 2, -1.1));
	
	std::vector<ArcPiece> set4_ans;
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));

	std::vector<ArcPiece> set4_pieces = set4.getChronoOrder();
	// set4.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set4_pieces, set4_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(MixedTime){
	Arcset set5(&sys);
	set5.addNode(Node(state1, 6, 2.2));
	set5.addNode(Node(state2, 6, 0));
	set5.addNode(Node(state3, 6, 1.1));
	set5.addNode(Node(state4, 6, -2.2));
	set5.addNode(Node(state5, 6, -1.1));
	set5.addSeg(Segment(4, 3, -1.1));
	set5.addSeg(Segment(1, 2, 1.1));
	set5.addSeg(Segment(1, 4, -1.1));
	set5.addSeg(Segment(2, 0, 1.1));
	
	std::vector<ArcPiece> set5_ans;
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));

	std::vector<ArcPiece> set5_pieces = set5.getChronoOrder();
	// set5.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set5_pieces, set5_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(MixedTime_LinkedSegs){
	// Mixed Time Set, Linked Segments
	Arcset set6(&sys);
	set6.addNode(Node(state1, 6, 0));
	set6.addNode(Node(state2, 6, 1.1));
	// set6.addNode(Node(state3, 2.2));
	// set6.addNode(Node(state4, -1.1));
	set6.addNode(Node(state5, 6, -2.2));
	set6.addSeg(Segment(0, 1, 1.1));
	set6.addSeg(Segment(1, -1, 1.1));
	set6.addSeg(Segment(2, -1, -1.1));

	double segLinkData[] = {2, 2, 2, NAN, NAN, NAN};
	Constraint segLinkCon(Constraint_tp::SEG_CONT_PV, 1, segLinkData, 6);
	set6.addConstraint(segLinkCon);
	// set6.print();

	std::vector<ArcPiece> set6_ans;
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	std::vector<ArcPiece> set6_pieces = set6.getChronoOrder();
	// set6.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set6_pieces, set6_ans));
}//====================================================

BOOST_AUTO_TEST_CASE(MixedTime_LinkedSegs2){
	// Mixed Time Set, Linked Segments, reversed constraint IDs
	Arcset set7(&sys);
	set7.addNode(Node(state1, 6, 0));
	set7.addNode(Node(state2, 6, 1.1));
	set7.addNode(Node(state5, 6, -2.2));
	set7.addSeg(Segment(0, 1, 1.1));
	set7.addSeg(Segment(1, -1, 1.1));
	set7.addSeg(Segment(2, -1, -1.1));

	double segLinkData2[] = {1, 1, 1, NAN, NAN, NAN};
	Constraint segLinkCon2(Constraint_tp::SEG_CONT_PV, 2, segLinkData2, 6);
	set7.addConstraint(segLinkCon2);

	std::vector<ArcPiece> set7_ans;
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	std::vector<ArcPiece> set7_pieces = set7.getChronoOrder();
	// set7.printInChrono();

	BOOST_CHECK(pieceVecsAreEqual(set7_pieces, set7_ans));
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Arcset_sortChrono)

BOOST_AUTO_TEST_CASE(ForwardTime){
	Arcset set1(&sys);
	set1.addNode(Node(state1, 6, 0));
	set1.addNode(Node(state2, 6, 1.1));
	set1.addNode(Node(state3, 6, 2.2));
	set1.addSeg(Segment(0, 1, 1.1));
	set1.addSeg(Segment(1, 2, 1.1));

	set1.sortChrono();
	
	// Node Order
	BOOST_CHECK_EQUAL(set1.getNodeByIx(0).getID(), 0); 
	BOOST_CHECK_EQUAL(set1.getNodeIx(0), 0);
	BOOST_CHECK_EQUAL(set1.getNodeByIx(1).getID(), 1); 
	BOOST_CHECK_EQUAL(set1.getNodeIx(1), 1);
	BOOST_CHECK_EQUAL(set1.getNodeByIx(2).getID(), 2);
	BOOST_CHECK_EQUAL(set1.getNodeIx(2), 2);

	// Segment Order
	BOOST_CHECK_EQUAL(set1.getSegByIx(0).getID(), 0);
	BOOST_CHECK_EQUAL(set1.getSegIx(0), 0);
	BOOST_CHECK_EQUAL(set1.getSegByIx(1).getID(), 1);
	BOOST_CHECK_EQUAL(set1.getSegIx(1), 1);
}//====================================================

BOOST_AUTO_TEST_CASE(ReverseTime){
	Arcset set2(&sys);
	set2.addNode(Node(state1, 6, 0));
	set2.addNode(Node(state2, 6, -1.1));
	set2.addNode(Node(state3, 6, -2.2));
	set2.addSeg(Segment(0, 1, -1.1));
	set2.addSeg(Segment(1, 2, -1.1));

	set2.sortChrono();
	BOOST_CHECK_EQUAL(set2.getNodeByIx(0).getID(), 2);
	BOOST_CHECK_EQUAL(set2.getNodeIx(2), 0);
	BOOST_CHECK_EQUAL(set2.getNodeByIx(1).getID(), 1);
	BOOST_CHECK_EQUAL(set2.getNodeIx(1), 1);
	BOOST_CHECK_EQUAL(set2.getNodeByIx(2).getID(), 0);
	BOOST_CHECK_EQUAL(set2.getNodeIx(0), 2);

	BOOST_CHECK_EQUAL(set2.getSegByIx(0).getID(), 1);
	BOOST_CHECK_EQUAL(set2.getSegIx(1), 0);
	BOOST_CHECK_EQUAL(set2.getSegByIx(1).getID(), 0);
	BOOST_CHECK_EQUAL(set2.getSegIx(0), 1);
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledForwardTime){
	Arcset set3(&sys);
	set3.addNode(Node(state2, 6, 1.1));
	set3.addNode(Node(state1, 6, 0));
	set3.addNode(Node(state3, 6, 2.2));
	set3.addSeg(Segment(1, 0, 1.1));
	set3.addSeg(Segment(0, 2, 1.1));
	
	set3.sortChrono();
	
	BOOST_CHECK_EQUAL(set3.getNodeByIx(0).getID(), 1);
	BOOST_CHECK_EQUAL(set3.getNodeIx(1), 0);
	BOOST_CHECK_EQUAL(set3.getNodeByIx(1).getID(), 0);
	BOOST_CHECK_EQUAL(set3.getNodeIx(0), 1);
	BOOST_CHECK_EQUAL(set3.getNodeByIx(2).getID(), 2);
	BOOST_CHECK_EQUAL(set3.getNodeIx(2), 2);

	BOOST_CHECK_EQUAL(set3.getSegByIx(0).getID(), 0);
	BOOST_CHECK_EQUAL(set3.getSegIx(0), 0);
	BOOST_CHECK_EQUAL(set3.getSegByIx(1).getID(), 1);
	BOOST_CHECK_EQUAL(set3.getSegIx(1), 1);
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledReverseTime){
	Arcset set4(&sys);
	set4.addNode(Node(state2, 6, -1.1));
	set4.addNode(Node(state1, 6, 0));
	set4.addNode(Node(state3, 6, -2.2));
	set4.addSeg(Segment(1, 0, -1.1));
	set4.addSeg(Segment(0, 2, -1.1));
	
	set4.sortChrono();
	
	BOOST_CHECK_EQUAL(set4.getNodeByIx(0).getID(), 2);
	BOOST_CHECK_EQUAL(set4.getNodeIx(2), 0);
	BOOST_CHECK_EQUAL(set4.getNodeByIx(1).getID(), 0);
	BOOST_CHECK_EQUAL(set4.getNodeIx(0), 1);
	BOOST_CHECK_EQUAL(set4.getNodeByIx(2).getID(), 1);
	BOOST_CHECK_EQUAL(set4.getNodeIx(1), 2);

	BOOST_CHECK_EQUAL(set4.getSegByIx(0).getID(), 1);
	BOOST_CHECK_EQUAL(set4.getSegIx(1), 0);
	BOOST_CHECK_EQUAL(set4.getSegByIx(1).getID(), 0);
	BOOST_CHECK_EQUAL(set4.getSegIx(0), 1);
}//====================================================

BOOST_AUTO_TEST_CASE(MixedTime){
	Arcset set5(&sys);
	set5.addNode(Node(state1, 6, 2.2));
	set5.addNode(Node(state2, 6, 0));
	set5.addNode(Node(state3, 6, 1.1));
	set5.addNode(Node(state4, 6, -2.2));
	set5.addNode(Node(state5, 6, -1.1));
	set5.addSeg(Segment(4, 3, -1.1));
	set5.addSeg(Segment(1, 2, 1.1));
	set5.addSeg(Segment(1, 4, -1.1));
	set5.addSeg(Segment(2, 0, 1.1));
	
	set5.sortChrono();
	BOOST_CHECK_EQUAL(set5.getNodeByIx(0).getID(), 3);
	BOOST_CHECK_EQUAL(set5.getNodeIx(3), 0);
	BOOST_CHECK_EQUAL(set5.getNodeByIx(1).getID(), 4);
	BOOST_CHECK_EQUAL(set5.getNodeIx(4), 1);
	BOOST_CHECK_EQUAL(set5.getNodeByIx(2).getID(), 1);
	BOOST_CHECK_EQUAL(set5.getNodeIx(1), 2);
	BOOST_CHECK_EQUAL(set5.getNodeByIx(3).getID(), 2);
	BOOST_CHECK_EQUAL(set5.getNodeIx(2), 3);
	BOOST_CHECK_EQUAL(set5.getNodeByIx(4).getID(), 0);
	BOOST_CHECK_EQUAL(set5.getNodeIx(0), 4);

	BOOST_CHECK_EQUAL(set5.getSegByIx(0).getID(), 0);
	BOOST_CHECK_EQUAL(set5.getSegIx(0), 0);
	BOOST_CHECK_EQUAL(set5.getSegByIx(1).getID(), 2);
	BOOST_CHECK_EQUAL(set5.getSegIx(2), 1);
	BOOST_CHECK_EQUAL(set5.getSegByIx(2).getID(), 1);
	BOOST_CHECK_EQUAL(set5.getSegIx(1), 2);
	BOOST_CHECK_EQUAL(set5.getSegByIx(3).getID(), 3);
	BOOST_CHECK_EQUAL(set5.getSegIx(3), 3);
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Arcset_AppendSet)

/**
 *  @brief Try appending sets; all these cases should work.
 */
BOOST_AUTO_TEST_CASE(Arcset_Append_ForForBegin){
	// Append (+) time set to beginning of (+) time set
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

	// Append end of forSet2 to beginning of forSet1 with tof = 1.3
	int segID = forSet1.appendSetAtNode(&forSet2, 0, 2, 1.3);

	// forSet1.print();
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), 1.3);
	BOOST_CHECK_CLOSE(forSet1.getNode(0).getEpoch(), 0, 1e-4);
	
	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}

}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForEnd){
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

	// Append beginning of forSet2 to end of forSet1 with tof = 1.3
	int segID = forSet1.appendSetAtNode(&forSet2, 2, 0, 1.3);
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), 1.3);
	BOOST_CHECK_CLOSE(forSet1.getNode(3).getEpoch(), 4.6, 1e-4);

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForBegin_ZeroTOF){
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

	// Append end of forSet2 to beginning of forSet1 with tof = 0
	int segID = forSet1.appendSetAtNode(&forSet2, 0, 2, 0);
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), 2.2);
	BOOST_CHECK_CLOSE(forSet1.getNode(0).getEpoch(), 0, 1e-4);

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForEnd_ZeroTOF){
	// Append (+) time set to end of (+) time set, TOF = 0
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

	// Append beginning of forSet2 to end of forSet1 with tof = 0
	int segID = forSet1.appendSetAtNode(&forSet2, 2, 0, 0);
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), 2.2);
	BOOST_CHECK_CLOSE(forSet1.getNode(3).getEpoch(), 3.3, 1e-4);

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegBegin){
	// Append (-) time set to beginning of (-) time set
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

	// Append beginning of revSet2 to the end of revSet1 with tof = -1.3
	int segID = revSet1.appendSetAtNode(&revSet2, 2, 0, -1.3);
	std::vector<ArcPiece> chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(revSet1.getSeg(segID).getTOF(), -1.3);
	BOOST_CHECK_CLOSE(revSet1.getNode(3).getEpoch(), -4.6, 1e-4);

	for(unsigned int s = 0; s < revSet1.getNumSegs(); s++){
		const Segment& seg = revSet1.getSegRefByIx_const(s);

		double termEpoch = revSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = revSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegEnd){
	// Append (-) time set to end of (-) time set
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

	// Append end of revSet2 to beginning of revSet1 with tof = -1.3
	int segID = revSet1.appendSetAtNode(&revSet2, 0, 2, -1.3);
	std::vector<ArcPiece> chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(revSet1.getSeg(segID).getTOF(), -1.3);
	BOOST_CHECK_CLOSE(revSet1.getNode(5).getEpoch(), 1.3, 1e-4);

	for(unsigned int s = 0; s < revSet1.getNumSegs(); s++){
		const Segment& seg = revSet1.getSegRefByIx_const(s);

		double termEpoch = revSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = revSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegBegin_ZeroTOF){
	// Append (-) time set to beginning of (-) time set, TOF = 0
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

	// append beginning of revSet2 to end of revSet1 with tof = 0
	int segID = revSet1.appendSetAtNode(&revSet2, 2, 0, 0);
	std::vector<ArcPiece> chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(revSet1.getSeg(segID).getTOF(), -2.2);
	BOOST_CHECK_CLOSE(revSet1.getNode(3).getEpoch(), -3.3, 1e-4);

	for(unsigned int s = 0; s < revSet1.getNumSegs(); s++){
		const Segment& seg = revSet1.getSegRefByIx_const(s);

		double termEpoch = revSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = revSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegEnd_ZeroTOF){
	// Append (-) time set to end of (-) time set, TOF = 0
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

	// Append the end of revSet2 to the beginning of revSet1
	int segID = revSet1.appendSetAtNode(&revSet2, 0, 2, 0);
	std::vector<ArcPiece> chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(revSet1.getSeg(segID).getTOF(), -2.2);
	BOOST_CHECK_CLOSE(revSet1.getNode(4).getEpoch(), 2.2, 1e-4);

	for(unsigned int s = 0; s < revSet1.getNumSegs(); s++){
		const Segment& seg = revSet1.getSegRefByIx_const(s);

		double termEpoch = revSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = revSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_PosTOF){
	// Append (-) time set to beginning of (+) time set, TOF > 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

	// Append beginning of revSet1 to beginning of forSet1 with TOF = 1.3
	int segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 1.3);

	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), 1.3);
	BOOST_CHECK_CLOSE(forSet1.getNode(3).getEpoch(), -1.3, 1e-4);

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_NegTOF){
	// Append (-) time set to beginning of (+) time set, TOF < 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

	// Append beginning of revSet1 to beginning of forSet1 with tof = -1.3
	int segID = forSet1.appendSetAtNode(&revSet1, 0, 0, -1.3);

	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), -1.3);

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_ZeroTOF){
	// Append (-) time set to beginning of (+) time set, TOF = 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

	// Append beginning of revSet1 to beginning of forSet1 with tof = 0
	int segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 0);
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK_EQUAL(forSet1.getSeg(segID).getTOF(), revSet1.getSegByIx(0).getTOF());

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_ZeroTOF_Short){
	// Append (-) time set to beginning of (+) time set, TOF = 0
	Arcset forSet1(&sys);
	forSet1.addNode(Node(state1, 6, 0));
	forSet1.addNode(Node(state2, 6, 1.1));
	Segment seg1(0, 1, 1.1);
	std::vector<double> t1 {0, 1.1};
	std::vector<double> q1(state1+0, state1+6);
	q1.insert(q1.end(), state2+0, state2+6);
	seg1.setTimeVector(t1);
	seg1.setStateVector(q1);
	forSet1.addSeg(seg1);

	Arcset revSet1(&sys);
	revSet1.addNode(Node(state1, 6, 0));
	revSet1.addNode(Node(state2, 6, -1.1));
	Segment seg2(0, 1, -1.1);
	std::vector<double> t2 {0, -1.1};
	std::vector<double> q2(state1+0, state1+6);
	q2.insert(q2.end(), state2+0, state2+6);
	seg2.setTimeVector(t2);
	seg2.setStateVector(q2);
	revSet1.addSeg(seg2);

	int segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 0);

	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == revSet1.getSegByIx(0).getTOF());

	for(unsigned int s = 0; s < forSet1.getNumSegs(); s++){
		const Segment& seg = forSet1.getSegRefByIx_const(s);

		double termEpoch = forSet1.getNode(seg.getTerminus()).getEpoch();
		double origEpoch = forSet1.getNode(seg.getOrigin()).getEpoch();
		BOOST_CHECK_SMALL(termEpoch - seg.getTOF() - origEpoch, 1e-8);

		// Check to make sure the segment times were updated with the node epochs
		BOOST_CHECK_SMALL(origEpoch - seg.getTimeByIx(0), 1e-8);
		BOOST_CHECK_SMALL(termEpoch - seg.getTimeByIx(-1), 1e-8);
	}
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Parallel01){
	// Append two (+) time sets at node 0
	initForwardSet();
	Arcset set1 = forwardSet;
	Arcset set2 = forwardSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 0, 0, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Parallel02){
	// Append two (+) time sets at end node
	initForwardSet();
	Arcset set1 = forwardSet;
	Arcset set2 = forwardSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 2, 2, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Parallel03){
	//Append two (+) time sets in the middle
	initForwardSet();
	Arcset set1 = forwardSet;
	Arcset set2 = forwardSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 2, 1, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision01){
	// Create time collision with (+) and (+) time sets
	initForwardSet();
	Arcset set1 = forwardSet;
	Arcset set2 = forwardSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 2, 0, -1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision02){
	// Create time collision with (+) and (+) time sets, again
	initForwardSet();
	Arcset set1 = forwardSet;
	Arcset set2 = forwardSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 0, 2, -1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision03){
	// Create time collision with (-) and (-) time sets
	initRevSet();
	Arcset set1 = revSet;
	Arcset set2 = revSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 2, 0, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision04){
	// Create time collision with (-) and (-) time sets, again
	initRevSet();
	Arcset set1 = revSet;
	Arcset set2 = revSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 0, 2, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision05){
	// Create time collision with (+) and (-) time sets
	initForwardSet();
	initRevSet();
	Arcset set1 = forwardSet;
	Arcset set2 = revSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 2, 0, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_AppendSet_Err_Collision06){
	// Create time collision with (+) and (-) time sets again
	initForwardSet();
	initRevSet();
	Arcset set1 = forwardSet;
	Arcset set2 = revSet;
	BOOST_CHECK_THROW(set1.appendSetAtNode(&set2, 0, 2, 1.1), Exception);
}//=================================================

BOOST_AUTO_TEST_SUITE_END()





