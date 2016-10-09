#define BOOST_TEST_MODULE ArcDataOps

#include <boost/test/unit_test.hpp>

#include "Exceptions.hpp"
#include "Linkable.hpp"
#include "Node.hpp"
// #include "Nodeset.hpp"
// #include "Segment.hpp"
// #include "SysData_cr3bp.hpp"
// #include "Utilities.hpp"

#include <iostream>

using namespace astrohelion;

// SysData_cr3bp sys("earth", "moon");
double state1[] = {1, 0, 0, 0, 0, 0};
// double state2[] = {2, 2, 0, 0, 0, 0};
// double state3[] = {3, 0, 3, 0, 0, 0};
// double state4[] = {4, 0, 0, 4, 0, 0};
// double state5[] = {5, 0, 0, 0, 5, 0};

int ivID = Linkable::INVALID_ID;

BOOST_AUTO_TEST_SUITE(Linkable)

BOOST_AUTO_TEST_CASE(Creation){
	Node n(state1, 10);
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
	Node n(state1, 10);
	n.addLink(3);
	n.addLink(7);

	BOOST_CHECK(n.isLinkedTo(3));
	BOOST_CHECK(n.isLinkedTo(7));
	BOOST_CHECK(!n.isLinkedTo(4));
}//====================================================

BOOST_AUTO_TEST_CASE(func_clearLinks){
	Node n(state1, 10);
	n.addLink(3);
	n.addLink(7);
	n.clearLinks();

	BOOST_CHECK(n.getLink(0) == ivID);
	BOOST_CHECK(n.getLink(1) == ivID);
}//====================================================

BOOST_AUTO_TEST_CASE(duplicateLinks){
	Node n(state1, 10);
	n.addLink(3);
	BOOST_CHECK_THROW(n.addLink(3), Exception);
}//====================================================

BOOST_AUTO_TEST_CASE(func_removeLink){
	// Test to make sure links are removed correctly
	Node n(state1, 10);
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