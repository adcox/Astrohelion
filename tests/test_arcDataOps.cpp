#include "Exceptions.hpp"
#include "Linkable.hpp"
#include "Node.hpp"
#include "Nodeset.hpp"
#include "Segment.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

#include <iostream>

using namespace astrohelion;

#define RESET   "\033[0m"					/**< Reset ASCII text to default values */
#define BOLDRED     "\033[1m\033[31m"		/**< Make ASCII text bold and red */
#define BOLDGREEN   "\033[1m\033[32m"		/**< Make ASCII text bold and green */

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

SysData_cr3bp sys("earth", "moon");
double state1[] = {1, 0, 0, 0, 0, 0};
double state2[] = {2, 2, 0, 0, 0, 0};
double state3[] = {3, 0, 3, 0, 0, 0};
double state4[] = {4, 0, 0, 4, 0, 0};
double state5[] = {5, 0, 0, 0, 5, 0};

/**
 *  @brief Create a basic Linkable object (in this case, a node) and test the functions
 *  to ensure everything works properly
 */
void testCreateLinkable(){
	int ivID = Linkable::INVALID_ID;

	Node n(state1, 10);

	std::cout << "Testing Linkable Creation:" << std::endl;
	std::cout << "  Initial links are emtpy: " << (n.getLink(0) == ivID && n.getLink(1) == ivID ? PASS : FAIL) << std::endl;

	// Add one link; second slot should still have INVALID_ID
	n.addLink(3);
	std::cout << "  One link added: " << (n.getLink(0) == 3 && n.getLink(1) == ivID ? PASS : FAIL) << std::endl;

	// Add a second link
	n.addLink(7);
	std::cout << "  Two links added: " << (n.getLink(0) == 3 && n.getLink(1) == 7 ? PASS : FAIL) << std::endl;
	
	// Check the isLinkedTo() function
	std::cout << "  isLinkedTo(): " << (n.isLinkedTo(3) ? PASS : FAIL) << std::endl;
	std::cout << "  isLinkedTo() 2: " << (n.isLinkedTo(4) ? FAIL : PASS) << std::endl;
	std::cout << "  isLinkedTo() 3: " << (n.isLinkedTo(7) ? PASS : FAIL) << std::endl;

	// Adding three links should result in an error
	std::cout << "  Three links added: ";
	try{
		n.addLink(8);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	// Test to make sure links are cleared
	n.clearLinks();
	std::cout << "  clearLinks(): " << (n.getLink(0) == ivID && n.getLink(1) == ivID ? PASS : FAIL) << std::endl;

	// Test to make sure attempting to add a duplicate link results in an error
	std::cout << "  Duplicate links added: ";
	try{
		n.addLink(3);
		n.addLink(3);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	// Test to make sure links are removed correctly
	n.clearLinks();
	n.addLink(3);
	n.addLink(7);
	n.removeLink(3);
	std::cout << "  removeLink(): " << (n.getLink(0) == ivID && n.getLink(1) == 7 ? PASS : FAIL) << std::endl;

	n.removeLink(3);	// Nothing should change
	std::cout << "  Remove invalid link: " << (n.getLink(0) == ivID && n.getLink(1) == 7 ? PASS : FAIL) << std::endl;
}//====================================================

/**
 *  @brief Create an arcset out of a set of nodes and segments; test basic_arcset functionality
 */
void testCreateArcset(){
	int ivID = Linkable::INVALID_ID;

	Nodeset set(&sys);

	Node n1(state1, 10);
	Node n2(state2, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	Segment s_bad1(n1ID, 999, 12);
	Segment s_bad2(n2ID, 987, 12);

	int sID = set.addSeg(s);

	Node node1 = set.getNode(n1ID);
	Node node2 = set.getNode(n2ID);

	std::cout << "Testing Arcset Creation:" << std::endl;
	std::cout << "  basic node initialization: " << (n1.getEpoch() == 10 && n1.getLink(0) == ivID && n1.getLink(1) == ivID ? PASS : FAIL) << std::endl;
	std::cout << "  basic segment initialization: " << (s.getTOF() == 15 && s.getLink(0) == n1ID && s.getLink(1) == n2ID ? PASS : FAIL) << std::endl;
	std::cout << "  nodes added correctly: " << (n1ID != ivID && n2ID != ivID ? PASS : FAIL) << std::endl;
	std::cout << "  segment added correctly: " << (sID != ivID ? PASS : FAIL) << std::endl;
	std::cout << "  linked node 1: " << (sID != ivID && node1.getLink(0) == sID && node1.getLink(1) == ivID ? PASS : FAIL) << std::endl;
	std::cout << "  linked node 2: " << (sID != ivID && node2.getLink(0) == sID && node2.getLink(1) == ivID ? PASS : FAIL) << std::endl;

	std::cout << "  bad segment 1: ";
	try{
		set.addSeg(s_bad1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  bad segment 2: ";
	try{
		set.addSeg(s_bad2);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}
}//==================================================

/**
 *  @brief Try to delete a segment from a simiple 2-node set
 */
void tryDeleteSeg(){
	int ivID = Linkable::INVALID_ID;

	Nodeset set(&sys);
	Node n1(state1, 10);
	Node n2(state2, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);

	int sID = set.addSeg(s);
	set.deleteSeg(sID);
	
	std::cout << "Delete segment: " << (set.getNumSegs() == 0 &&
		set.getNode(n1ID).getLink(0) == ivID && set.getNode(n2ID).getLink(0) == ivID ? PASS : FAIL) << std::endl;
}//==================================================

/**
 *  @brief Try to delete the first node in a set; this should throw an exception
 * 
 *  @return wether or not the test was successful
 */
bool tryDeleteFirstNode(){
	Nodeset set(&sys);
	Node n1(state1, 10);
	Node n2(state2, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	set.addSeg(s);

	try{
		set.deleteNode(0);
		return true;
	}catch(Exception &e){
		return false;
	}
}//==================================================

/**
 *  @brief Try to delete the last node in a set; this should throw an exception
 * 
 *  @return wether or not the test was successful
 */
bool tryDeleteLastNode(){
	Nodeset set(&sys);
	Node n1(state1, 10);
	Node n2(state2, 25);

	int n1ID = set.addNode(n1);
	int n2ID = set.addNode(n2);

	Segment s(n1ID, n2ID, 15);
	set.addSeg(s);
	try{
		set.deleteNode(set.getNumNodes()-1);
		return true;
	}catch(Exception &e){
		return false;
	}
}//==================================================

/**
 *  @brief Attempt to delete a node from the middle of a linear-time arcset
 *  @details Make sure exceptions are thrown when they should be and that the arcset
 *  "heals" itself properly
 */
void tryDeleteMiddleNode(){
	Nodeset set(&sys);
	set.addNode(Node(state1, 0));
	set.addNode(Node(state2, 1.1));
	set.addNode(Node(state3, 2.2));
	set.addNode(Node(state4, 3.3));
	set.addNode(Node(state5, 4.4));
	set.addSeg(Segment(0, 1, 1.1));
	set.addSeg(Segment(1, 2, 1.1));
	set.addSeg(Segment(2, 3, 1.1));
	int sID_last = set.addSeg(Segment(3, 4, 1.1));

	std::cout << "Delete middle node (LINEAR, FORWARD TIME): ";
	try{
		set.deleteNode(1);
		Segment seg = set.getSeg(sID_last+1);	// should retrieve the new segment
		if(seg.getOrigin() == 0 && seg.getTerminus() == 2 && std::abs(seg.getTOF() - 2.2) < 1e-4)
			std::cout << PASS << std::endl;
		else
			std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "  Correct, updated nodeIDMap: " << (set.getNodeIx(0) == 0 && 
		set.getNodeIx(2) == 1 &&
		set.getNodeIx(3) == 2 && set.getNodeIx(4) == 3 ? PASS : FAIL) << std::endl;
	
	std::cout << "  Checking deleted node... ";
	try{
		set.getNodeIx(1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Correct, updated segIDMap: " << (set.getSegIx(2) == 0 && set.getSegIx(3) == 1 ? PASS : FAIL) << std::endl;
	std::cout << "  Checking deleted segments... ";
	try{
		set.getSegIx(0);
		set.getSegIx(1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Delete another middle node: ";
	try{
		set.deleteNode(2);
		Segment seg = set.getSeg(sID_last+2);
		if(seg.getOrigin() == 0 && seg.getTerminus() == 3 && std::abs(seg.getTOF() - 3.3) < 1e-4)
			std::cout << PASS << std::endl;
		else{
			std::cout << FAIL << std::endl;
		}
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}
	std::cout << "  Correct, updated nodeIDMap: " << (set.getNodeIx(0) == 0 && 
		set.getNodeIx(3) == 1 && set.getNodeIx(4) == 2 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted node... ";
	try{
		set.getNodeIx(1);
		set.getNodeIx(2);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Correct, updated segIDMap: " << (set.getSegIx(3) == 0 && set.getSegIx(5) == 1 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted segments ... ";
	try{
		set.getSegIx(0);
		set.getSegIx(1);
		set.getSegIx(2);
		set.getSegIx(4);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}
}//==================================================

/**
 *  @brief Same procedure as tryDeleteMiddleNode(), but using an arcset that
 *  progresses linearly in reverse time
 */
void tryDeleteMiddleNode_revTime(){
	Nodeset set(&sys);
	set.addNode(Node(state1, 0));
	set.addNode(Node(state2, -1.1));
	set.addNode(Node(state3, -2.2));
	set.addNode(Node(state4, -3.3));
	set.addNode(Node(state5, -4.4));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	set.addSeg(Segment(2, 3, -1.1));
	int sID_last = set.addSeg(Segment(3, 4, -1.1));

	std::cout << "Delete middle node (LINEAR, REVERSE TIME): ";
	try{
		set.deleteNode(1);
		// set.print();
		Segment seg = set.getSeg(sID_last+1);	// should retrieve the new segment
		if(seg.getOrigin() == 0 && seg.getTerminus() == 2 && std::abs(seg.getTOF() + 2.2) < 1e-4)
			std::cout << PASS << std::endl;
		else
			std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "  Correct, updated nodeIDMap: " << (set.getNodeIx(0) == 0 && 
		set.getNodeIx(2) == 1 &&
		set.getNodeIx(3) == 2 && set.getNodeIx(4) == 3 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted node... ";
	try{
		set.getNodeIx(1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Correct, updated segIDMap: " << (set.getSegIx(2) == 0 &&
		set.getSegIx(3) == 1 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted segments... ";
	try{
		set.getSegIx(0);
		set.getSegIx(1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Delete another middle node: ";
	try{
		set.deleteNode(2);
		// set.print();
		Segment seg = set.getSeg(sID_last+2);
		if(seg.getOrigin() == 0 && seg.getTerminus() == 3 && std::abs(seg.getTOF() + 3.3) < 1e-4)
			std::cout << PASS << std::endl;
		else{
			std::cout << FAIL << std::endl;
		}
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "  Correct, updated nodeIDMap: " << (set.getNodeIx(0) == 0 && 
		set.getNodeIx(3) == 1 && set.getNodeIx(4) == 2 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted node... ";
	try{
		set.getNodeIx(1);
		set.getNodeIx(2);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Correct, updated segIDMap: " << (
		set.getSegIx(3) == 0 && set.getSegIx(5) == 1 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted segments... ";
	try{
		set.getSegIx(0);
		set.getSegIx(1);
		set.getSegIx(2);
		set.getSegIx(4);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}
}//====================================================

/**
 *  @brief Try to delete a node that is the origin of two arcs.
 *  @details Make sure the arcset is healed properly.
 */
void tryDeleteMiddleNode_doubleSource1(){
	Nodeset set(&sys);
	set.addNode(Node(state1, 0));
	set.addNode(Node(state2, -1.1));
	set.addNode(Node(state3, -2.2));
	set.addNode(Node(state4, 1.1));
	set.addNode(Node(state5, 2.2));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	set.addSeg(Segment(0, 3, 1.1));
	int sID_last = set.addSeg(Segment(3, 4, 1.1));

	std::cout << "Delete middle node (DOULBE DIRECTION 1): ";
	try{
		set.deleteNode(0);
		// set.print();
		Segment seg = set.getSeg(sID_last+1);
		if(seg.getOrigin() == 1 && seg.getTerminus() == 3 && std::abs(seg.getTOF() - 2.2) < 1e-4)
			std::cout << PASS << std::endl;
		else
			std::cout << FAIL << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "  Correct, updated nodeIDMap: " << (set.getNodeIx(1) == 0 &&
		set.getNodeIx(2) == 1 &&
		set.getNodeIx(3) == 2 && set.getNodeIx(4) == 3 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted node... ";
	try{
		set.getNodeIx(0);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}

	std::cout << "  Correct, updated segIDMap: " << (set.getSegIx(1) == 0 &&
		set.getSegIx(3) == 1 && set.getSegIx(4) == 2 ? PASS : FAIL) << std::endl;

	std::cout << "  Checking deleted segments... ";
	try{
		set.getSegIx(0);
		set.getSegIx(2);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}
}//==================================================

/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the forward time segment has no terminal point
 */
void tryDeleteMiddleNode_doubleSource2(){
	Nodeset set(&sys);
	set.addNode(Node(state1, 0));
	set.addNode(Node(state2, -1.1));
	set.addNode(Node(state3, -2.2));
	set.addSeg(Segment(0, 1, -1.1));
	set.addSeg(Segment(1, 2, -1.1));
	int sID_last = set.addSeg(Segment(0, Linkable::INVALID_ID, 1.1));

	std::cout << "Delete middle node (DOULBE DIRECTION 2): ";
	try{
		set.deleteNode(0);
		// set.print();
		Segment seg = set.getSeg(sID_last+1);
		if(seg.getOrigin() == 1 && seg.getTerminus() == Linkable::INVALID_ID && std::abs(seg.getTOF() - 2.2) < 1e-4)
			std::cout << PASS << std::endl;
		else
			std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}
}//==================================================

/**
 *  @brief Try to delete a node that is the origin of two arcs where
 *  the reverse time segment has no terminal point
 */
void tryDeleteMiddleNode_doubleSource3(){
	Nodeset set(&sys);
	set.addNode(Node(state1, 0));
	set.addNode(Node(state2, 1.1));
	set.addNode(Node(state3, 2.2));
	set.addSeg(Segment(0, Linkable::INVALID_ID, -1.1));
	set.addSeg(Segment(0, 1, 1.1));
	int sID_last = set.addSeg(Segment(1, 2, 1.1));

	std::cout << "Delete middle node (DOULBE DIRECTION 3): ";
	try{
		set.deleteNode(0);
		// set.print();
		Segment seg = set.getSeg(sID_last+1);
		if(seg.getOrigin() == 1 && seg.getTerminus() == Linkable::INVALID_ID && std::abs(seg.getTOF() + 2.2) < 1e-4)
			std::cout << PASS << std::endl;
		else
			std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}
}//==================================================

bool pieceVecsAreEqual(std::vector<ArcPiece> v1, std::vector<ArcPiece> v2){
	if(v1.size() != v2.size())
		return false;

	for(unsigned int i = 0; i < v1.size(); i++){
		if(v1[i] != v2[i])
			return false;
	}

	return true;
}//====================================================

void testPutInChrono(){
	std::cout << "Testing getChronoOrder() function:" << std::endl;
	// Forward time set
	Nodeset set1(&sys);
	set1.addNode(Node(state1, 0));
	set1.addNode(Node(state2, 1.1));
	set1.addNode(Node(state3, 2.2));
	set1.addSeg(Segment(0, 1, 1.1));
	set1.addSeg(Segment(1, 2, 1.1));

	std::vector<ArcPiece> set1_ans;
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set1_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::vector<ArcPiece> set1_pieces = set1.getChronoOrder();
	set1.printInChrono();

	std::cout << "  Forward Time Set: " << (pieceVecsAreEqual(set1_pieces, set1_ans) ? PASS : FAIL) << std::endl;

	// Reverse time set
	Nodeset set2(&sys);
	set2.addNode(Node(state1, 0));
	set2.addNode(Node(state2, -1.1));
	set2.addNode(Node(state3, -2.2));
	set2.addSeg(Segment(0, 1, -1.1));
	set2.addSeg(Segment(1, 2, -1.1));

	std::vector<ArcPiece> set2_ans;
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set2_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));

	std::vector<ArcPiece> set2_pieces = set2.getChronoOrder();
	set2.printInChrono();

	std::cout << "  Reverse Time Set: " << (pieceVecsAreEqual(set2_pieces, set2_ans) ? PASS : FAIL) << std::endl;

	// Shuffled forward time
	Nodeset set3(&sys);
	set3.addNode(Node(state2, 1.1));
	set3.addNode(Node(state1, 0));
	set3.addNode(Node(state3, 2.2));
	set3.addSeg(Segment(1, 0, 1.1));
	set3.addSeg(Segment(0, 2, 1.1));
	
	std::vector<ArcPiece> set3_ans;
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set3_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::vector<ArcPiece> set3_pieces = set3.getChronoOrder();
	set3.printInChrono();

	std::cout << "  Shuffled forward Time Set: " << (pieceVecsAreEqual(set3_pieces, set3_ans) ? PASS : FAIL) << std::endl;

	// Shuffled reverse time
	Nodeset set4(&sys);
	set4.addNode(Node(state2, -1.1));
	set4.addNode(Node(state1, 0));
	set4.addNode(Node(state3, -2.2));
	set4.addSeg(Segment(1, 0, -1.1));
	set4.addSeg(Segment(0, 2, -1.1));
	
	std::vector<ArcPiece> set4_ans;
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set4_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));

	std::vector<ArcPiece> set4_pieces = set4.getChronoOrder();
	set4.printInChrono();

	std::cout << "  Shuffled reverse Time Set: " << (pieceVecsAreEqual(set4_pieces, set4_ans) ? PASS : FAIL) << std::endl;

	// Mixed time set
	Nodeset set5(&sys);
	set5.addNode(Node(state1, 2.2));
	set5.addNode(Node(state2, 0));
	set5.addNode(Node(state3, 1.1));
	set5.addNode(Node(state4, -2.2));
	set5.addNode(Node(state5, -1.1));
	set5.addSeg(Segment(4, 3, -1.1));
	set5.addSeg(Segment(1, 2, 1.1));
	set5.addSeg(Segment(1, 4, -1.1));
	set5.addSeg(Segment(2, 0, 1.1));
	
	std::vector<ArcPiece> set5_ans;
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	set5_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));

	std::vector<ArcPiece> set5_pieces = set5.getChronoOrder();
	set5.printInChrono();
	std::cout << "  Shuffled Mixed Time Set: " << (pieceVecsAreEqual(set5_pieces, set5_ans) ? PASS : FAIL) << std::endl;

	// Mixed Time Set, Linked Segments
	Nodeset set6(&sys);
	set6.addNode(Node(state1, 0));
	set6.addNode(Node(state2, 1.1));
	// set6.addNode(Node(state3, 2.2));
	// set6.addNode(Node(state4, -1.1));
	set6.addNode(Node(state5, -2.2));
	set6.addSeg(Segment(0, 1, 1.1));
	set6.addSeg(Segment(1, -1, 1.1));
	set6.addSeg(Segment(2, -1, -1.1));

	double segLinkData[] = {2, 2, 2, NAN, NAN, NAN};
	Constraint segLinkCon(Constraint_tp::SEG_CONT_PV, 1, segLinkData, 6);
	set6.addConstraint(segLinkCon);
	// set6.print();

	std::vector<ArcPiece> set6_ans;
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	set6_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::vector<ArcPiece> set6_pieces = set6.getChronoOrder();
	set6.printInChrono();
	std::cout << "  Linked Seg Set: " << (pieceVecsAreEqual(set6_pieces, set6_ans) ? PASS : FAIL) << std::endl;

	// Mixed Time Set, Linked Segments, reversed constraint IDs
	Nodeset set7(&sys);
	set7.addNode(Node(state1, 0));
	set7.addNode(Node(state2, 1.1));
	set7.addNode(Node(state5, -2.2));
	set7.addSeg(Segment(0, 1, 1.1));
	set7.addSeg(Segment(1, -1, 1.1));
	set7.addSeg(Segment(2, -1, -1.1));

	double segLinkData2[] = {1, 1, 1, NAN, NAN, NAN};
	Constraint segLinkCon2(Constraint_tp::SEG_CONT_PV, 2, segLinkData2, 6);
	set7.addConstraint(segLinkCon2);

	std::vector<ArcPiece> set7_ans;
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	set7_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::vector<ArcPiece> set7_pieces = set7.getChronoOrder();
	set7.printInChrono();
	std::cout << "  Linked Seg Set #2: " << (pieceVecsAreEqual(set7_pieces, set7_ans) ? PASS : FAIL) << std::endl;
}//====================================================

void testPutInChrono2(){
	std::cout << "Testing putInChronoOrder() function:" << std::endl;

	// Forward time set
	Nodeset set1(&sys);
	set1.addNode(Node(state1, 0));
	set1.addNode(Node(state2, 1.1));
	set1.addNode(Node(state3, 2.2));
	set1.addSeg(Segment(0, 1, 1.1));
	set1.addSeg(Segment(1, 2, 1.1));

	set1.putInChronoOrder();
	bool nodeOrder = set1.getNodeByIx(0).getID() == 0 && 
					set1.getNodeByIx(1).getID() == 1 && 
					set1.getNodeByIx(2).getID() == 2;
	bool segOrder = set1.getSegByIx(0).getID() == 0 &&
					set1.getSegByIx(1).getID() == 1;

	std::cout << "  Forward Time Set: " << (nodeOrder && segOrder ? PASS : FAIL) << std::endl;

	// Reverse time set
	Nodeset set2(&sys);
	set2.addNode(Node(state1, 0));
	set2.addNode(Node(state2, -1.1));
	set2.addNode(Node(state3, -2.2));
	set2.addSeg(Segment(0, 1, -1.1));
	set2.addSeg(Segment(1, 2, -1.1));

	set2.putInChronoOrder();
	nodeOrder = set2.getNodeByIx(0).getID() == 2 && 
				set2.getNodeByIx(1).getID() == 1 && 
				set2.getNodeByIx(2).getID() == 0;
	segOrder = set2.getSegByIx(0).getID() == 1 &&
				set2.getSegByIx(1).getID() == 0;

	std::cout << "  Reverse Time Set: " << (nodeOrder && segOrder ? PASS : FAIL) << std::endl;

	// Shuffled forward time
	Nodeset set3(&sys);
	set3.addNode(Node(state2, 1.1));
	set3.addNode(Node(state1, 0));
	set3.addNode(Node(state3, 2.2));
	set3.addSeg(Segment(1, 0, 1.1));
	set3.addSeg(Segment(0, 2, 1.1));
	
	set3.putInChronoOrder();
	nodeOrder = set3.getNodeByIx(0).getID() == 1 && 
				set3.getNodeByIx(1).getID() == 0 && 
				set3.getNodeByIx(2).getID() == 2;
	segOrder = set3.getSegByIx(0).getID() == 0 &&
				set3.getSegByIx(1).getID() == 1;

	std::cout << "  Shuffled forward Time Set: " << (nodeOrder && segOrder ? PASS : FAIL) << std::endl;

	// Shuffled reverse time
	Nodeset set4(&sys);
	set4.addNode(Node(state2, -1.1));
	set4.addNode(Node(state1, 0));
	set4.addNode(Node(state3, -2.2));
	set4.addSeg(Segment(1, 0, -1.1));
	set4.addSeg(Segment(0, 2, -1.1));
	
	set4.putInChronoOrder();
	nodeOrder = set4.getNodeByIx(0).getID() == 2 && 
				set4.getNodeByIx(1).getID() == 0 && 
				set4.getNodeByIx(2).getID() == 1;
	segOrder = set4.getSegByIx(0).getID() == 1 &&
				set4.getSegByIx(1).getID() == 0;

	std::cout << "  Shuffled reverse Time Set: " << (nodeOrder && segOrder ? PASS : FAIL) << std::endl;

	// Mixed time set
	Nodeset set5(&sys);
	set5.addNode(Node(state1, 2.2));
	set5.addNode(Node(state2, 0));
	set5.addNode(Node(state3, 1.1));
	set5.addNode(Node(state4, -2.2));
	set5.addNode(Node(state5, -1.1));
	set5.addSeg(Segment(4, 3, -1.1));
	set5.addSeg(Segment(1, 2, 1.1));
	set5.addSeg(Segment(1, 4, -1.1));
	set5.addSeg(Segment(2, 0, 1.1));
	
	set5.putInChronoOrder();
	nodeOrder = set5.getNodeByIx(0).getID() == 3 && 
				set5.getNodeByIx(1).getID() == 4 && 
				set5.getNodeByIx(2).getID() == 1 &&
				set5.getNodeByIx(3).getID() == 2 && 
				set5.getNodeByIx(4).getID() == 0;
	segOrder = set5.getSegByIx(0).getID() == 0 &&
				set5.getSegByIx(1).getID() == 2 && 
				set5.getSegByIx(2).getID() == 1 && 
				set5.getSegByIx(3).getID() == 3;

	std::cout << "  Shuffled Mixed Time Set: " << (nodeOrder && segOrder ? PASS : FAIL) << std::endl;
}//====================================================

/**
 *  @brief Try appending sets; all these cases should work.
 */
void tryAppendSet(){
	Nodeset forwardSet(&sys);
	forwardSet.addNode(Node(state1, 0));
	forwardSet.addNode(Node(state2, 1.1));
	forwardSet.addNode(Node(state3, 2.2));
	forwardSet.addSeg(Segment(0, 1, 1.1));
	forwardSet.addSeg(Segment(1, 2, 1.1));

	Nodeset revSet(&sys);
	revSet.addNode(Node(state1, 0));
	revSet.addNode(Node(state2, -1.1));
	revSet.addNode(Node(state3, -2.2));
	revSet.addSeg(Segment(0, 1, -1.1));
	revSet.addSeg(Segment(1, 2, -1.1));

	//****************************************
	Nodeset forSet1 = forwardSet;
	Nodeset forSet2 = forwardSet;

	int segID = forSet1.appendSetAtNode(&forSet2, 0, 2, 1.3);
	// forSet1.print();
	std::vector<ArcPiece> chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	std::vector<ArcPiece> chrono_ans;
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::cout << "Append (+) time set to beginning of (+) time set: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == 1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	forSet2 = forwardSet;

	segID = forSet1.appendSetAtNode(&forSet2, 2, 0, 1.3);
	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));

	std::cout << "Append (+) time set to end of (+) time set: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == 1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	forSet2 = forwardSet;

	segID = forSet1.appendSetAtNode(&forSet2, 0, 2, 0);
	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::cout << "Append (+) time set to beginning of (+) time set, TOF = 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == 1.1 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	forSet2 = forwardSet;

	segID = forSet1.appendSetAtNode(&forSet2, 2, 0, 0);
	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));

	std::cout << "Append (+) time set to end of (+) time set, TOF = 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == 1.1 ? PASS : FAIL) << std::endl;

	//****************************************
	Nodeset revSet1 = revSet;
	Nodeset revSet2 = revSet;

	segID = revSet1.appendSetAtNode(&revSet2, 2, 0, -1.3);
	chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));

	std::cout << "Append (-) time set to beginning of (-) time set: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (revSet1.getSeg(segID).getTOF() == -1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	revSet1 = revSet;
	revSet2 = revSet;

	segID = revSet1.appendSetAtNode(&revSet2, 0, 2, -1.3);
	chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));

	std::cout << "Append (-) time set to end of (-) time set: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (revSet1.getSeg(segID).getTOF() == -1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	revSet1 = revSet;
	revSet2 = revSet;

	segID = revSet1.appendSetAtNode(&revSet2, 2, 0, 0);
	chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));

	std::cout << "Append (-) time set to beginning of (-) time set, TOF = 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (revSet1.getSeg(segID).getTOF() == -1.1 ? PASS : FAIL) << std::endl;

	//****************************************
	revSet1 = revSet;
	revSet2 = revSet;

	segID = revSet1.appendSetAtNode(&revSet2, 0, 2, 0);
	chrono = revSet1.getChronoOrder();
	// revSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, segID));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));

	std::cout << "Append (-) time set to end of (-) time set, TOF = 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (revSet1.getSeg(segID).getTOF() == -1.1 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	revSet1 = revSet;

	segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 1.3);

	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::cout << "Append (-) time set to beginning of (+) time set, TOF > 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == 1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	revSet1 = revSet;

	segID = forSet1.appendSetAtNode(&revSet1, 0, 0, -1.3);

	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 5));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::cout << "Append (-) time set to beginning of (+) time set, TOF < 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == -1.3 ? PASS : FAIL) << std::endl;

	//****************************************
	forSet1 = forwardSet;
	revSet1 = revSet;

	segID = forSet1.appendSetAtNode(&revSet1, 0, 0, 0);

	chrono = forSet1.getChronoOrder();
	// forSet1.printInChrono();

	chrono_ans.clear();
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 4));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 2));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 3));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 0));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::SEG, 1));
	chrono_ans.push_back(ArcPiece(ArcPiece::Piece_Tp::NODE, 2));

	std::cout << "Append (-) time set to beginning of (+) time set, TOF = 0: " << (pieceVecsAreEqual(chrono, chrono_ans) ? PASS : FAIL) << std::endl;
	std::cout << "  >> New segment has correct TOF: " << (forSet1.getSeg(segID).getTOF() == revSet1.getSegByIx(1).getTOF() ? PASS : FAIL) << std::endl;
}//=================================================

/**
 *  @brief Try appending sets; all these cases should throw errors
 */
void tryAppendSet_errs(){
	Nodeset forwardSet(&sys);
	forwardSet.addNode(Node(state1, 0));
	forwardSet.addNode(Node(state2, 1.1));
	forwardSet.addNode(Node(state3, 2.2));
	forwardSet.addSeg(Segment(0, 1, 1.1));
	forwardSet.addSeg(Segment(1, 2, 1.1));

	Nodeset revSet(&sys);
	revSet.addNode(Node(state1, 0));
	revSet.addNode(Node(state2, -1.1));
	revSet.addNode(Node(state3, -2.2));
	revSet.addSeg(Segment(0, 1, -1.1));
	revSet.addSeg(Segment(1, 2, -1.1));

	std::cout << "Append two (+) time sets at node 0: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = forwardSet;

		set1.appendSetAtNode(&set2, 0, 0, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Append two (+) time sets at end node: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = forwardSet;

		set1.appendSetAtNode(&set2, 2, 2, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Append two (+) time sets in the middle: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = forwardSet;

		set1.appendSetAtNode(&set2, 2, 1, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (+) and (+) time sets: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = forwardSet;

		set1.appendSetAtNode(&set2, 2, 0, -1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (+) and (+) time sets, again: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = forwardSet;

		set1.appendSetAtNode(&set2, 0, 2, -1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (-) and (-) time sets: ";
	try{
		Nodeset set1 = revSet;
		Nodeset set2 = revSet;

		set1.appendSetAtNode(&set2, 2, 0, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (-) and (-) time sets, again: ";
	try{
		Nodeset set1 = revSet;
		Nodeset set2 = revSet;

		set1.appendSetAtNode(&set2, 0, 2, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (+) and (-) time sets: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = revSet;

		set1.appendSetAtNode(&set2, 2, 0, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}

	std::cout << "Create time collision with (+) and (-) time sets again: ";
	try{
		Nodeset set1 = forwardSet;
		Nodeset set2 = revSet;

		set1.appendSetAtNode(&set2, 0, 2, 1.1);
		std::cout << FAIL << std::endl;
	}catch(Exception &e){
		std::cout << PASS << std::endl;
	}catch(std::exception &e){
		std::cout << FAIL << std::endl;
		printf("  %s\n", e.what());
	}
}//=================================================

int main(){
	testCreateLinkable();
	testCreateArcset();
	tryDeleteSeg();
	std::cout << "Delete first node: " << (tryDeleteFirstNode() ? PASS : FAIL) << std::endl;
	std::cout << "Delete last node: " << (tryDeleteLastNode() ? PASS : FAIL) << std::endl;
	tryDeleteMiddleNode();
	tryDeleteMiddleNode_revTime();
	tryDeleteMiddleNode_doubleSource1();
	tryDeleteMiddleNode_doubleSource2();
	tryDeleteMiddleNode_doubleSource3();

	testPutInChrono();
	testPutInChrono2();
	tryAppendSet_errs();
	tryAppendSet();
}//====================================================


//