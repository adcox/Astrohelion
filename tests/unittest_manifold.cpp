#define BOOST_TEST_MODULE ManifoldEngine

#include <vector>

#include <boost/test/unit_test.hpp>

#include "Calculations.hpp"
#include "ManifoldEngine.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

SysData_cr3bp emSys("earth", "moon");

/* Check the eigenvector / eigenvalue function to make sure it returns the
 * correct number and type of eigenvalue/eigenvector pairs
 */
BOOST_AUTO_TEST_CASE(EIG_VEC_VAL){
	// Small Lyapunov
	std::vector<double> ic {0.839367079288131, 0, 0, 0, -0.0201563145344225, 0};
	double period = 2.7;
	Arcset_cr3bp perOrbit = cr3bp_getPeriodic(&emSys, ic, period, Mirror_tp::MIRROR_XZ);
	
	ManifoldEngine engine;
	std::vector<cdouble> eigVals;

	// **************************************
	// Stable Eigenvector
	// **************************************
	MatrixXRd eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S_P, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) < 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S_M, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) < 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) < 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	// **************************************
	// Unstable Eigenvector
	// **************************************
	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U_P, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) > 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U_M, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) > 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK(std::abs(eigVals[0]) > 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	// **************************************
	// Stable and Unstable Eigenvector
	// **************************************
	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_ALL, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 2);
	BOOST_CHECK(std::abs(eigVals[0]) < 1);	// Stable first
	BOOST_CHECK(std::abs(eigVals[1]) > 1);	// Unstable second
	BOOST_CHECK(eigVecs.cols() == 2);
}//====================================================

BOOST_AUTO_TEST_CASE(Single_Manifold){
	// Small Lyapunov
	std::vector<double> ic {0.839367079288131, 0, 0, 0, -0.0201563145344225, 0};
	double period = 2.7;
	Arcset_cr3bp perOrbit = cr3bp_getPeriodic(&emSys, ic, period, Mirror_tp::MIRROR_XZ);
	perOrbit.saveToMat("data/lyap.mat");
	ManifoldEngine engine;

	// **************************************
	// Stable and Unstable Manifolds
	// **************************************
	std::vector<Arcset_cr3bp> manifolds = engine.computeSingleFromPeriodic(Manifold_tp::MAN_ALL,
		&perOrbit, period, 1, Manifold_StepOff_tp::STEP_MATCH_JC);

	BOOST_CHECK(manifolds.size() == 4);

	unsigned int numRev = 0, numFor = 0;
	for(unsigned int i = 0; i < manifolds.size(); i++){
		numRev += manifolds[i].getTimeByIx(-1) < 0;
		numFor += manifolds[i].getTimeByIx(-1) > 0;
		BOOST_CHECK(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1) < 1e-8);
		BOOST_CHECK(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)) < 1e-8);
	}

	BOOST_CHECK(numRev == 2);
	BOOST_CHECK(numFor == 2);

	// **************************************
	// Stable Manifolds
	// **************************************
	manifolds = engine.computeSingleFromPeriodic(Manifold_tp::MAN_S, &perOrbit, period, 1, Manifold_StepOff_tp::STEP_MATCH_JC);

	BOOST_CHECK(manifolds.size() == 2);
	numRev = 0;
	numFor = 0;
	for(unsigned int i = 0; i < manifolds.size(); i++){
		numRev += manifolds[i].getTimeByIx(-1) < 0;
		numFor += manifolds[i].getTimeByIx(-1) > 0;
		BOOST_CHECK(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1) < 1e-8);
		BOOST_CHECK(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)) < 1e-8);
	}

	BOOST_CHECK(numRev == 2);
	BOOST_CHECK(numFor == 0);

	// **************************************
	// Unstable Manifolds
	// **************************************
	manifolds = engine.computeSingleFromPeriodic(Manifold_tp::MAN_U, &perOrbit, period, 1, Manifold_StepOff_tp::STEP_MATCH_JC);

	BOOST_CHECK(manifolds.size() == 2);
	numRev = 0;
	numFor = 0;
	for(unsigned int i = 0; i < manifolds.size(); i++){
		numRev += manifolds[i].getTimeByIx(-1) < 0;
		numFor += manifolds[i].getTimeByIx(-1) > 0;
		BOOST_CHECK(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1) < 1e-8);
		BOOST_CHECK(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)) < 1e-8);
	}

	BOOST_CHECK(numRev == 0);
	BOOST_CHECK(numFor == 2);
}//====================================================




