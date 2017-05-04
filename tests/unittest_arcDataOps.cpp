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
	forwardSet.addNode(Node(state3, 6, 2.2));
	forwardSet.addSeg(Segment(0, 1, 1.1));
	forwardSet.addSeg(Segment(1, 2, 1.1));
}//====================================================

void initRevSet(){
	revSet = Arcset(&sys);
	revSet.addNode(Node(state1, 6, 0));
	revSet.addNode(Node(state2, 6, -1.1));
	revSet.addNode(Node(state3, 6, -2.2));
	revSet.addSeg(Segment(0, 1, -1.1));
	revSet.addSeg(Segment(1, 2, -1.1));
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

	BOOST_CHECK(set.getNumSegs() == 0);
	BOOST_CHECK(set.getNode(n1ID).getLink(0) == ivID);
	BOOST_CHECK(set.getNode(n2ID).getLink(0) == ivID);
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
	BOOST_CHECK(set.getSeg(s1ID).getOrigin() == ivID);
	BOOST_CHECK(set.getSeg(s1ID).getTerminus() == n2ID);
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
	BOOST_CHECK(set.getSeg(s3ID).getOrigin() == n3ID);
	BOOST_CHECK(set.getSeg(s3ID).getTerminus() == ivID);
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
	BOOST_CHECK(seg.getOrigin() == 1);
	BOOST_CHECK(seg.getTerminus() == 3);
	BOOST_CHECK(std::abs(seg.getTOF() - 2.2) < 1e-4);

	// Correct, updated nodeIDMap
	BOOST_CHECK(set.getNodeIx(1) == 0);
	BOOST_CHECK(set.getNodeIx(2) == 1);
	BOOST_CHECK(set.getNodeIx(3) == 2);
	BOOST_CHECK(set.getNodeIx(4) == 3);

	// Make sure the deleted node is, in fact, deleted
	BOOST_CHECK_THROW(set.getNodeIx(0), Exception);

	// Correct, updated segIDMap
	BOOST_CHECK(set.getSegIx(1) == 0);
	BOOST_CHECK(set.getSegIx(3) == 1);
	BOOST_CHECK(set.getSegIx(4) == 2);

	// Check deleted segments
	BOOST_CHECK_THROW(set.getSegIx(0), Exception);
	BOOST_CHECK_THROW(set.getSegIx(2), Exception);
}//==================================================


/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the forward time segment has no terminal point
 */
BOOST_AUTO_TEST_CASE(deleteMiddleNOde_doubleSource2){
	Arcset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, -1.1));
	set.addNode(Node(state3, 6, -2.2));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	int sID_last = set.addSeg(Segment(0, ivID, 1.1));

	set.deleteNode(0);
	Segment seg = set.getSeg(sID_last+1);

	BOOST_CHECK(seg.getOrigin() == 1);
	BOOST_CHECK(seg.getTerminus() == ivID);
	BOOST_CHECK(std::abs(seg.getTOF() - 2.2) < 1e-4);
}//==================================================

/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the reverse time segment has no terminal point
 */
BOOST_AUTO_TEST_CASE(deleteMiddleNOde_doubleSource3){
	Arcset set(&sys);
	set.addNode(Node(state1, 6, 0));
	set.addNode(Node(state2, 6, 1.1));
	set.addNode(Node(state3, 6, 2.2));
	set.addSeg(Segment(0, ivID, -1.1));
	set.addSeg(Segment(0, 1, 1.1));
	int sID_last = set.addSeg(Segment(1, 2, 1.1));

	set.deleteNode(0);
	Segment seg = set.getSeg(sID_last+1);

	BOOST_CHECK(seg.getOrigin() == 1);
	BOOST_CHECK(seg.getTerminus() == ivID);
	BOOST_CHECK(std::abs(seg.getTOF() + 2.2) < 1e-4);
}//==================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* Arcset Creation
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

BOOST_AUTO_TEST_SUITE(Arcset_putInChronoOrder)

BOOST_AUTO_TEST_CASE(ForwardTime){
	Arcset set1(&sys);
	set1.addNode(Node(state1, 6, 0));
	set1.addNode(Node(state2, 6, 1.1));
	set1.addNode(Node(state3, 6, 2.2));
	set1.addSeg(Segment(0, 1, 1.1));
	set1.addSeg(Segment(1, 2, 1.1));

	set1.putInChronoOrder();
	
	// Node Order
	BOOST_CHECK(set1.getNodeByIx(0).getID() == 0); 
	BOOST_CHECK(set1.getNodeByIx(1).getID() == 1); 
	BOOST_CHECK(set1.getNodeByIx(2).getID() == 2);

	// Segment Order
	BOOST_CHECK(set1.getSegByIx(0).getID() == 0);
	BOOST_CHECK(set1.getSegByIx(1).getID() == 1);
}//====================================================

BOOST_AUTO_TEST_CASE(ReverseTime){
	Arcset set2(&sys);
	set2.addNode(Node(state1, 6, 0));
	set2.addNode(Node(state2, 6, -1.1));
	set2.addNode(Node(state3, 6, -2.2));
	set2.addSeg(Segment(0, 1, -1.1));
	set2.addSeg(Segment(1, 2, -1.1));

	set2.putInChronoOrder();
	BOOST_CHECK(set2.getNodeByIx(0).getID() == 2); 
	BOOST_CHECK(set2.getNodeByIx(1).getID() == 1);
	BOOST_CHECK(set2.getNodeByIx(2).getID() == 0);
	
	BOOST_CHECK(set2.getSegByIx(0).getID() == 1);
	BOOST_CHECK(set2.getSegByIx(1).getID() == 0);
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledForwardTime){
	Arcset set3(&sys);
	set3.addNode(Node(state2, 6, 1.1));
	set3.addNode(Node(state1, 6, 0));
	set3.addNode(Node(state3, 6, 2.2));
	set3.addSeg(Segment(1, 0, 1.1));
	set3.addSeg(Segment(0, 2, 1.1));
	
	set3.putInChronoOrder();
	
	BOOST_CHECK(set3.getNodeByIx(0).getID() == 1);
	BOOST_CHECK(set3.getNodeByIx(1).getID() == 0);
	BOOST_CHECK(set3.getNodeByIx(2).getID() == 2);
	
	BOOST_CHECK(set3.getSegByIx(0).getID() == 0);
	BOOST_CHECK(set3.getSegByIx(1).getID() == 1);
}//====================================================

BOOST_AUTO_TEST_CASE(ShuffledReverseTime){
	Arcset set4(&sys);
	set4.addNode(Node(state2, 6, -1.1));
	set4.addNode(Node(state1, 6, 0));
	set4.addNode(Node(state3, 6, -2.2));
	set4.addSeg(Segment(1, 0, -1.1));
	set4.addSeg(Segment(0, 2, -1.1));
	
	set4.putInChronoOrder();
	
	BOOST_CHECK(set4.getNodeByIx(0).getID() == 2);
	BOOST_CHECK(set4.getNodeByIx(1).getID() == 0);
	BOOST_CHECK(set4.getNodeByIx(2).getID() == 1);
	
	BOOST_CHECK(set4.getSegByIx(0).getID() == 1);
	BOOST_CHECK(set4.getSegByIx(1).getID() == 0);
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
	
	set5.putInChronoOrder();
	BOOST_CHECK(set5.getNodeByIx(0).getID() == 3);
	BOOST_CHECK(set5.getNodeByIx(1).getID() == 4);
	BOOST_CHECK(set5.getNodeByIx(2).getID() == 1);
	BOOST_CHECK(set5.getNodeByIx(3).getID() == 2);
	BOOST_CHECK(set5.getNodeByIx(4).getID() == 0);
	
	BOOST_CHECK(set5.getSegByIx(0).getID() == 0);
	BOOST_CHECK(set5.getSegByIx(1).getID() == 2);
	BOOST_CHECK(set5.getSegByIx(2).getID() == 1);
	BOOST_CHECK(set5.getSegByIx(3).getID() == 3);
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
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == 1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForEnd){
	// Append (+) time set to end of (+) time set
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

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
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == 1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForBegin_ZeroTOF){
	// Append (+) time set to beginning of (+) time set, TOF = 0
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

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
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == 1.1);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_ForForEnd_ZeroTOF){
	// Append (+) time set to end of (+) time set, TOF = 0
	initForwardSet();
	Arcset forSet1 = forwardSet;
	Arcset forSet2 = forwardSet;

	int segID = forSet1.appendSetAtNode(&forSet2, 2, 0, 0);
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == 1.1);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegBegin){
	// Append (-) time set to beginning of (-) time set
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

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
	BOOST_CHECK(revSet1.getSeg(segID).getTOF() == -1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegEnd){
	// Append (-) time set to end of (-) time set
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

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
	BOOST_CHECK(revSet1.getSeg(segID).getTOF() == -1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegBegin_ZeroTOF){
	// Append (-) time set to beginning of (-) time set, TOF = 0
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

	int segID = revSet1.appendSetAtNode(&revSet2, 2, 0, 0);
	std::vector<ArcPiece> chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
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
	BOOST_CHECK(revSet1.getSeg(segID).getTOF() == -1.1);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegNegEnd_ZeroTOF){
	// Append (-) time set to end of (-) time set, TOF = 0
	initRevSet();
	Arcset revSet1 = revSet;
	Arcset revSet2 = revSet;

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
	BOOST_CHECK(revSet1.getSeg(segID).getTOF() == -1.1);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_PosTOF){
	// Append (-) time set to beginning of (+) time set, TOF > 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

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
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == 1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_NegTOF){
	// Append (-) time set to beginning of (+) time set, TOF < 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

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
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == -1.3);
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_ZeroTOF){
	// Append (-) time set to beginning of (+) time set, TOF = 0
	initForwardSet();
	initRevSet();
	Arcset forSet1 = forwardSet;
	Arcset revSet1 = revSet;

	int segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 0);

	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, 2));

	BOOST_CHECK(pieceVecsAreEqual(chrono, chrono_ans));
	BOOST_CHECK(forSet1.getSeg(segID).getTOF() == revSet1.getSegByIx(1).getTOF());
}//=================================================

BOOST_AUTO_TEST_CASE(Arcset_Append_NegForBegin_ZeroTOF_Short){
	// Append (-) time set to beginning of (+) time set, TOF = 0
	Arcset forSet1(&sys);
	forSet1.addNode(Node(state1, 6, 0));
	forSet1.addNode(Node(state2, 6, 1.1));
	forSet1.addSeg(Segment(0, 1, 1.1));

	Arcset revSet1(&sys);
	revSet1.addNode(Node(state1, 6, 0));
	revSet1.addNode(Node(state2, 6, -1.1));
	revSet1.addSeg(Segment(0, 1, -1.1));

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





