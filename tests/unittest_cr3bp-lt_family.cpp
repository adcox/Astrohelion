#define BOOST_TEST_MODULE CR3BP_LT_FAMILY

#include <boost/test/unit_test.hpp>
#include <vector>

#include "Arcset_periodic.hpp"
#include "Arcset_cr3bp_lt.hpp"
#include "Family_PO_cr3bp_lt.hpp"
#include "SysData_cr3bp_lt.hpp"

using namespace astrohelion;

void freeLaws(std::vector<ControlLaw*>&);

void freeLaws(std::vector<ControlLaw*>& laws){
	for(auto p : laws){
		delete p;
		p = nullptr;
	}
}//====================================================

/*	
 *	Find orbits with the specified alpha angle
 */
BOOST_AUTO_TEST_CASE(FIND_2D_ANGLE){
	char filename[] = "../../data/families/cr3bp-lt_earth-moon/"
		"E1_Lyap_f7.0e-02_Hlt-1.558_law2112.mat";
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt fam(&sys);
	std::vector<ControlLaw*> laws;
	fam.readFromMat(filename, laws);

	double alpha = PI/5;
	std::vector<Arcset_periodic> matches = fam.getMemberBy2DThrustAngle(alpha);
	for(unsigned int i = 0; i < matches.size(); i++){
		std::vector<double> ctrl0 = matches[i].getNodeRefByIx(0).\
			getExtraParamVec(PARAMKEY_CTRL);

		BOOST_CHECK_SMALL(ctrl0[0] - alpha, 1e-9);
	}

	freeLaws(laws);
}//====================================================

/*
 *	Find orbits with the specified low-thrust Hamiltonian value
 */
BOOST_AUTO_TEST_CASE(FIND_HLT){
	char filename[] = "../../data/families/cr3bp-lt_earth-moon/"
		"E1_Lyap_f7.0e-02_alph000.00_law2112.mat";
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt fam(&sys);
	std::vector<ControlLaw*> laws;
	fam.readFromMat(filename, laws);

	double H = -1.4;
	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(H);
	for(unsigned int i = 0; i < matches.size(); i++){
		Arcset_cr3bp_lt temp(matches[i]);
		BOOST_CHECK_SMALL(temp.getHltByIx(0) - H, 1e-9);
	}

	freeLaws(laws);
}//====================================================