/**
 *	Generate a Lissajous trajectory using linearization, discretize it, and then correct in the CR3BP
 *	
 *	Compile by running:
 *		
 *		>> make generateLis
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_linear_motion_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"

#include <cstdio>

int main(){

	tpat_sys_data_cr3bp sys("sun", "earth");
	
	// Define physical characteristics of Lissajous Orbit
	double Ay = 1000000/sys.getCharL();		// In-plane amplitude
	double Az = 150000/sys.getCharL();		// Out-of-plane amplitude
	double phi = 0, psi = 0;				// Phase angles

	// Create a linear motion engine
	tpat_linear_motion_engine linEngine;
	linEngine.setNumRevs(35);

	// Create a linear approximation to a Lissajous
	tpat_traj_cr3bp SE_L1_Liss_linear = linEngine.getCR3BPLiss(1, Ay, false, phi, Az, psi, &sys);

	// Discretize the trajectory w/o using integration
	tpat_nodeset_cr3bp SE_L1_Liss_linear_nodes = SE_L1_Liss_linear.discretize(250);

	// Save to mat files to check with Matlab
	SE_L1_Liss_linear.saveToMat("data/SE_L1_Liss_linear.mat");
	SE_L1_Liss_linear_nodes.saveToMat("data/SE_L1_Liss_linear_nodes.mat");

	// Compute mean Jacobi for linear trajector
	double sumC = 0;
	for(int i = 0; i < SE_L1_Liss_linear.getLength(); i++){
		sumC += SE_L1_Liss_linear.getJacobi(i);
	}
	double meanC = sumC/SE_L1_Liss_linear.getLength();

	// Constraint Jacobi of non-linear Liss to be the average JC of the linear Liss
	// tpat_constraint fixC(tpat_constraint::JC, std::floor(SE_L1_Liss_linear_nodes.getNumNodes()/2), &meanC, 1);
	// SE_L1_Liss_linear_nodes.addConstraint(fixC);

	// Attempt to correct Liss into a continuous arc
	tpat_correction_engine corrector;
	corrector.setVarTime(true);
	corrector.setEqualArcTime(true);
	corrector.setMaxIts(50);

	corrector.multShoot(&SE_L1_Liss_linear_nodes);

	tpat_nodeset_cr3bp SE_L1_Liss_nonlin_nodes = corrector.getCR3BP_Output();
	tpat_traj_cr3bp SE_L1_Liss_nonlin = tpat_traj_cr3bp::fromNodeset(SE_L1_Liss_nonlin_nodes);

	SE_L1_Liss_nonlin_nodes.saveToMat("data/SE_L1_Liss_nonlin_nodes.mat");
	SE_L1_Liss_nonlin.saveToMat("data/SE_L1_liss_nonlin.mat");

	return EXIT_SUCCESS;
}