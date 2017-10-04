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
	engine.setVerbosity(Verbosity_tp::NO_MSG);
	std::vector<cdouble> eigVals;

	// **************************************
	// Stable Eigenvector
	// **************************************
	MatrixXRd eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S_RIGHT, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 1);
	BOOST_CHECK_LT(std::abs(eigVals[0]), 1);
	BOOST_CHECK_EQUAL(eigVecs.cols(), 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S_LEFT, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 1);
	BOOST_CHECK_LT(std::abs(eigVals[0]), 1);
	BOOST_CHECK_EQUAL(eigVecs.cols(), 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_S, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 1);
	BOOST_CHECK_LT(std::abs(eigVals[0]), 1);
	BOOST_CHECK_EQUAL(eigVecs.cols(), 1);

	// **************************************
	// Unstable Eigenvector
	// **************************************
	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U_RIGHT, &perOrbit, &eigVals);

	BOOST_CHECK(eigVals.size() == 1);
	BOOST_CHECK_GT(std::abs(eigVals[0]), 1);
	BOOST_CHECK(eigVecs.cols() == 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U_LEFT, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 1);
	BOOST_CHECK_GT(std::abs(eigVals[0]), 1);
	BOOST_CHECK_EQUAL(eigVecs.cols(), 1);

	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_U, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 1);
	BOOST_CHECK_GT(std::abs(eigVals[0]), 1);
	BOOST_CHECK_EQUAL(eigVecs.cols(), 1);

	// **************************************
	// Stable and Unstable Eigenvector
	// **************************************
	eigVals.clear();
	eigVecs = engine.eigVecValFromPeriodic(Manifold_tp::MAN_ALL, &perOrbit, &eigVals);

	BOOST_CHECK_EQUAL(eigVals.size(), 2);
	BOOST_CHECK_LT(std::abs(eigVals[0]), 1);	// Stable first
	BOOST_CHECK_GT(std::abs(eigVals[1]),  1);	// Unstable second
	BOOST_CHECK_EQUAL(eigVecs.cols(), 2);
}//====================================================

BOOST_AUTO_TEST_CASE(Single_Manifold){
	// Small Lyapunov
	std::vector<double> ic {0.839367079288131, 0, 0, 0, -0.0201563145344225, 0};
	double period = 2.7;
	Arcset_cr3bp perOrbit = cr3bp_getPeriodic(&emSys, ic, period, Mirror_tp::MIRROR_XZ);
	perOrbit.saveToMat("data/lyap.mat");
	ManifoldEngine engine;
	engine.setVerbosity(Verbosity_tp::NO_MSG);

	// **************************************
	// Stable and Unstable Manifolds
	// **************************************
	std::vector<Arcset_cr3bp> manifolds = engine.computeSingleFromPeriodic(Manifold_tp::MAN_ALL,
		&perOrbit, period, 1, Manifold_StepOff_tp::STEP_MATCH_JC);

	BOOST_CHECK_EQUAL(manifolds.size(), 4);

	unsigned int numRev = 0, numFor = 0;
	for(unsigned int i = 0; i < manifolds.size(); i++){
		numRev += manifolds[i].getTimeByIx(-1) < 0;
		numFor += manifolds[i].getTimeByIx(-1) > 0;
		BOOST_CHECK_SMALL(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1), 1e-8);
		BOOST_CHECK_SMALL(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)), 1e-8);
	}

	BOOST_CHECK_EQUAL(numRev, 2);
	BOOST_CHECK_EQUAL(numFor, 2);

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
		BOOST_CHECK_SMALL(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1), 1e-8);
		BOOST_CHECK_SMALL(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)), 1e-8);
	}

	BOOST_CHECK_EQUAL(numRev, 2);
	BOOST_CHECK_EQUAL(numFor, 0);

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
		BOOST_CHECK_SMALL(std::abs(std::abs(manifolds[i].getTimeByIx(-1)) - 1), 1e-8);
		BOOST_CHECK_SMALL(std::abs(manifolds[i].getJacobiByIx(0) - perOrbit.getJacobiByIx(0)), 1e-8);
	}

	BOOST_CHECK_EQUAL(numRev, 0);
	BOOST_CHECK_EQUAL(numFor, 2);
}//====================================================

BOOST_AUTO_TEST_CASE(Many_Manifold){
	// Small Lyapunov
	std::vector<double> ic {0.839367079288131, 0, 0, 0, -0.0201563145344225, 0};
	double period = 2.7;
	Arcset_cr3bp perOrbit = cr3bp_getPeriodic(&emSys, ic, period, Mirror_tp::MIRROR_XZ);
	unsigned int numPts = 100;
	
	Arcset_cr3bp discretizedOrbit(&emSys);
	SimEngine sim;
	sim.setVarStepSize(false);
	sim.setNumSteps(numPts);
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(perOrbit.getStateByIx(0), perOrbit.getTotalTOF(), numPts+1, &discretizedOrbit);
	// discretizedOrbit.saveToMat("data/discretizedOrbit.mat");

	ManifoldEngine engine;
	std::vector<Arcset_cr3bp> manifolds = engine.computeSetFromPeriodic(Manifold_tp::MAN_ALL,
		&discretizedOrbit, numPts, 0.1, Manifold_StepOff_tp::STEP_MATCH_JC);

	BOOST_CHECK_EQUAL(manifolds.size(), 4*numPts);

	std::vector<double> allICs;
	unsigned int numRev = 0, numFor = 0;
	for(unsigned int i = 0; i < numPts; i++){
		for(unsigned int m = 0; m < 4; m++){
			numRev += manifolds[4*i+m].getTimeByIx(-1) < 0;
			numFor += manifolds[4*i+m].getTimeByIx(-1) > 0;

			// Make sure the manifold initial condition is near the corresponding fixed point
			std::vector<double> manifoldIC = manifolds[4*i+m].getStateByIx(0);
			std::vector<double> orbitState = discretizedOrbit.getStateByIx(i);
			double dist = sqrt(pow(manifoldIC[0] - orbitState[0], 2) + pow(manifoldIC[1] - orbitState[1], 2) + 
				pow(manifoldIC[2] - orbitState[2], 2));

			allICs.insert(allICs.end(), manifoldIC.begin(), manifoldIC.end());
			BOOST_CHECK_LT(std::abs(dist*emSys.getCharL() - engine.getStepOffDist()), 0.1);
		}
	}

	// mat_t *matfp = Mat_CreateVer("data/returns.mat", NULL, MAT_FT_DEFAULT);
	// saveMatrixToFile(matfp, "ics", allICs, allICs.size()/6, 6);
	// Mat_Close(matfp);

	BOOST_CHECK_EQUAL(numRev, 2*numPts);
	BOOST_CHECK_EQUAL(numFor, 2*numPts);	
}//====================================================




