/**
 *	Generate a Lissajous trajectory "Tom's Method"
 *	
 *	Compile by running:
 *		
 *		>> make generateLisTom
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
 
#include <cmath>
#include <cstdio>

int main(){

	tpat_sys_data_cr3bp sys("earth", "moon");	

	// Load Vertical family
	tpat_family_cr3bp vertFam("../share/families/EM_L1_Vert.mat");

	// Get a vertical orbit
	tpat_family_member_cr3bp vertMember = vertFam.getMember(5);

	// Compute Nodeset from family member
	tpat_nodeset_cr3bp vertNodes(vertMember.getIC(), vertFam.getSysDataPtr(), vertMember.getTOF(), 7);

	// Stack the nodes on top of each other
	int numStack = 30;
	tpat_nodeset_cr3bp stackedNodes = vertNodes;
	for(int n = 0; n < numStack; n++){
		tpat_nodeset_cr3bp temp = stackedNodes + vertNodes;
		stackedNodes = temp;
	}

	// Constraint first node to have same x as original vertical, step away in y, fix z to 0
	double data[] = {vertMember.getIC().at(0), vertMember.getIC().at(1) + 0.001, 0, NAN, NAN, NAN};
	tpat_constraint initNodeCon(tpat_constraint::STATE, 0, data, 6);

	stackedNodes.addConstraint(initNodeCon);

	// Correct the entire trajectory
	tpat_correction_engine corrector;
	corrector.correct(&stackedNodes);

	tpat_nodeset_cr3bp correctedNodes = corrector.getCR3BP_Output();
	correctedNodes.saveToMat("data/Liss_TomMethod_Nodes.mat");
	
	tpat_traj_cr3bp lissTraj = tpat_traj_cr3bp::fromNodeset(correctedNodes);
	lissTraj.saveToMat("data/Liss_TomMethod_Traj.mat");
	
	return EXIT_SUCCESS;
}