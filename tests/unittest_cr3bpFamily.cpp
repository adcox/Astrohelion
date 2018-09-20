#define BOOST_TEST_MODULE CR3BP_Family

#include <boost/test/unit_test.hpp>

/**
 *	Test out CR3BP Families of Orbits
 */
#include "Family_PO_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "SysData_cr3bp.hpp"

using namespace astrohelion;

BOOST_AUTO_TEST_CASE(FAMILY_OPERATIONS){
	// Load the family
	SysData_cr3bp sys("../../data/families/cr3bp_earth-moon/L1_Lyap.mat");
	Family_PO_cr3bp fam(&sys);
	std::vector<ControlLaw*> laws;
	fam.readFromMat("../../data/families/cr3bp_earth-moon/L1_Lyap.mat", laws);

	// fam.sortEigs();
	// std::vector<unsigned int> bifs = fam.findBifurcations();
	// if(bifs.size() > 0){
	// 	printf("Found bifurcations at:\n");
	// 	for(unsigned int i = 0; i < bifs.size(); i++){
	// 		printf("  Ix = %04d\n", bifs[i]);
	// 	}
	// }
	// fam.saveToMat("data/LoadedButterflyFam.mat");	// Check to see if data was re-loaded correctly

	printf("Checing Match State: X\n");
	double matchX = 0.9;
	std::vector<Arcset_periodic> matches = fam.getMemberByState(matchX, 0);
	for(unsigned int i = 0; i < matches.size(); i++){
		// printf("   %03u: x = %f ", i, matches[i].getStateByIx(0)[0]);
		// std::cout << (std::abs(matches[i].getIC()[0] - matchX) < 1e-9 ? PASS : FAIL) << std::endl;
		BOOST_CHECK_SMALL(matches[i].getStateByIx(0)[0] - matchX, 1e-9);
	}

	printf("Checking Match Jacobi:\n");
	double matchJC = 3;
	matches.clear();
	matches = fam.getMemberByJacobi(matchJC);
	// printf("  Found %zu Potential members\n", matches.size());
	for(unsigned int i = 0; i < matches.size(); i++){
		// printf("   %03u: JC = %f ", i, matches[i].getJacobi());
		// std::cout << (std::abs(matches[i].getJacobi() - matchJC) < 1e-9 ? PASS : FAIL) << std::endl;
		Arcset_cr3bp temp(matches[i]);
		BOOST_CHECK_SMALL(temp.getJacobiByIx(0) - matchJC, 1e-9);
	}

	printf("Checking Match TOF:\n");
	double matchTOF = 4.5;
	matches.clear();
	matches = fam.getMemberByTOF(matchTOF);
	// printf("  Found %zu Potential members\n", matches.size());
	for(unsigned int i = 0; i < matches.size(); i++){
		// printf("   %03u: TOF = %f ", i, matches[i].getTOF());
		// std::cout << (std::abs(matches[i].getTOF() - matchTOF) < 1e-9 ? PASS : FAIL) << std::endl;
		BOOST_CHECK_SMALL(matches[i].getTotalTOF() - matchTOF, 1e-9);
	}

	for(unsigned int i = 0; i < laws.size(); i++){
		if(laws[i]){
			delete laws[i];
			laws[i] = nullptr;
		}
	}
}