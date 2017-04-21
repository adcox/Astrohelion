#define BOOST_TEST_MODULE ManifoldEngine

#include <vector>

#include <boost/test/unit_test.hpp>

#include "Calculations.hpp"
#include "ManifoldEngine.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

SysData_cr3bp emSys("earth", "moon");

/* Check the eigenvector / eigenvalue function to make sure it returns the
 * correct number and type of eigenvalue/eigenvector pairs
 */
BOOST_AUTO_TEST_CASE(EIG_VEC_VAL){
	std::vector<double> ic {0.839367079288131, 0, 0, 0, -0.0201563145344225, 0};
	double period = 2.7;
	Traj_cr3bp perOrbit = cr3bp_getPeriodic(&emSys, ic, period, Mirror_tp::MIRROR_XZ);
	
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
}