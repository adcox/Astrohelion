#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

int main(int argc, char** argv){
	// Check the A-matrix for the CR3BP-LT with GEN_INERT
	SysData_cr3bp_lt sys("earth", "moon", 100);
	unsigned int lawID = ltlaw::GEN_INERT | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	std::vector<double> params {1.234, 1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawID, params);

	double q0[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.123, 0.321, 
		1, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 1};
	double t = 0.4;

	EOM_ParamStruct eomParams(&sys, &law);

	MatrixXRd A_analytic(9,9), A_numeric(9,9);

	// Get the A matrix analytically
	double dq[90] = {0};
	DynamicsModel_cr3bp_lt::fullEOMs(t, q0, dq, &eomParams);
	A_analytic = Eigen::Map<MatrixXRd>(dq+9, 9, 9);

	double err = 1e-7;
	for(unsigned int i = 0; i < 9; i++){
		double qp[90] = {0}, qdot[90] = {0}, qdot2[90] = {0};

		// Perturb the state vector in forward direction
		std::copy(q0, q0+90, qp);
		qp[i] += err;

		// Compute the state derivative
		DynamicsModel_cr3bp_lt::fullEOMs(t, qp, qdot, &eomParams);
		Eigen::VectorXd Af = Eigen::Map<Eigen::VectorXd>(qdot, 9, 1);

		// Perturb the state vector in the backward direction
		std::copy(q0, q0+90, qp);
		qp[i] -= err;
		
		// Compute state derivatives
		DynamicsModel_cr3bp_lt::fullEOMs(t, qp, qdot2, &eomParams);
		Eigen::VectorXd Ab = Eigen::Map<Eigen::VectorXd>(qdot2, 9, 1);

		A_numeric.col(i) = (Af - Ab)/(2*err);
	}

	MatrixXRd diff = A_analytic - A_numeric;

	std::cout << "Analytic:\n" << A_analytic << std::endl;
	std::cout << "Numeric:\n" << A_numeric << std::endl;
	std::cout << "Difference:\n" << diff.cwiseAbs() << std::endl;
	return EXIT_SUCCESS;
}