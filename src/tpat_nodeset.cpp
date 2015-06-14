/**
 *	@file tpat_nodeset.cpp
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "tpat.hpp"

#include "tpat_nodeset.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_trajectory.hpp"
#include "tpat_utilities.hpp"
 
#include <cstdio>
#include <iostream>
#include <cmath>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a nodeset with space for one node and TOF
 *	@param n the size of a node (i.e. number of states)
 */
tpat_nodeset::tpat_nodeset(const int n) : nodeSize(n){
	// Reserve space for one node and TOF
	nodes.reserve(n);
	tofs.reserve(1);
}//========================================================

/**
 *	@brief Copy constructor
 *	@param n a nodeset
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset& n) : nodeSize(n.nodeSize){
	if(nodeSize == n.nodeSize){
		nodeDistro = n.nodeDistro;
		nodes = n.nodes;
		tofs = n.tofs;
		velConNodes = n.velConNodes;
		constraints = n.constraints;
		velConSet = n.velConSet;
	}else{
		throw tpat_exception("Nodesets must have same nodeSize");
	}
}//========================================================

/**
 *	@brief Destructor
 */
tpat_nodeset::~tpat_nodeset(){
	nodes.clear();
	tofs.clear();
	constraints.clear();
	velConNodes.clear();
}//======================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param n a different nodeset
 *	@return this nodeset, modified to be equal to n
 */
tpat_nodeset& tpat_nodeset::operator =(const tpat_nodeset &n){
	if(nodeSize == n.nodeSize){
		nodeDistro = n.nodeDistro;
		nodes = n.nodes;
		tofs = n.tofs;
		velConNodes = n.velConNodes;
		constraints = n.constraints;
		velConSet = n.velConSet;
		return *this;
	}else{
		throw tpat_exception("Nodesets must have same size");
	}
}//======================================================

/**
 *	@brief Concatenate two nodesets
 *	@param lhs the "first" nodeset
 *	@param rhs the "second" nodeset
 *	@param output the output nodeset; will be made equal to lhs + rhs, or [lhs, rhs]
 */
void tpat_nodeset::basicConcat(const tpat_nodeset &lhs, const tpat_nodeset &rhs, tpat_nodeset *output){
	// how much the index of the nodes in RHS will change when concatenated
	int rhsNodeShift = lhs.getNumNodes();
	// If the last node of LHS equals the first node of RHS, skip it when concatenating
	if(lhs.getNode(-1) == rhs.getNode(0)){
		rhsNodeShift--;
		vector<double> lhsNodes(lhs.nodes.begin(), lhs.nodes.end()-1*lhs.nodeSize);
		output->appendNode(lhsNodes);
	}else{
		output->appendNode(lhs.nodes);
	}
	output->appendNode(rhs.nodes);

	// Concatenate TOF vectors
	output->tofs.insert(output->tofs.end(), lhs.tofs.begin(), lhs.tofs.end());
	output->tofs.insert(output->tofs.end(), rhs.tofs.begin(), rhs.tofs.end());

	// Concatenate constraint vectors, adjusting node numbers on any constraints in RHS
	// to account for the shift in node index
	for(int c = 0; c < lhs.getNumCons(); c++){
		output->addConstraint(lhs.constraints[c]);
	}
	for(int c = 0; c < rhs.getNumCons(); c++){
		tpat_constraint temp(rhs.getConstraint(c));
		temp.setNode(temp.getNode() + rhsNodeShift);
		output->addConstraint(temp);
	}

	// Concatenate velocity continuity nodes, adjusting node numbers for RHS
	// velConNodes to account for the shift in node index
	vector<int> newVelCon(lhs.velConNodes.begin(), lhs.velConNodes.end());
	for(int n = 0; n < ((int)rhs.velConNodes.size()); n++){
		newVelCon.push_back(rhs.velConNodes[n] + rhsNodeShift);
	}
	output->setVelConNodes(newVelCon);

	output->setNodeDistro( lhs.nodeDistro == rhs.nodeDistro ? lhs.nodeDistro : NONE);
}//=======================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return a pointer to the beginning of the node array; useful for in-place editing
 */
std::vector<double>* tpat_nodeset::getNodes(){ return &nodes; }

/**
 *	@return a pointer to the beginning of the TOF array; useful for in-place editing
 */
std::vector<double>* tpat_nodeset::getTOFs(){ return &tofs; }

/**
 *	@param i the index of the node (begins at zero). Negative values count
 *	backwards from the end, i.e. -1 gives the final value, -2 gives the second
 *	to last value, etc.
 *	@return the node specified by <tt>i</tt>
 */
std::vector<double> tpat_nodeset::getNode(int i) const {
	if(i < 0)
		i += getNumNodes();

	vector<double> temp(nodes.begin()+i*nodeSize, nodes.begin()+(i+1)*nodeSize);
	return temp;
}

/**
 *	@param i the index of the TOF; note that if there are n nodes, there will be n-1 
 *	TOFs. Negative values count backwards from the end, i.e. -1 gives the final value,
 *	-2 gives the second to last value, etc.
 *	@return the TOF between nodes <tt>i</tt> and <tt>i+1</tt>
 */
double tpat_nodeset::getTOF(int i) const {
	if(i < 0)
		i += tofs.size();
	return tofs.at(i);
}

/**
 *	@return the number of nodes contained in this node set
 */
int tpat_nodeset::getNumNodes() const { return (int)(nodes.size()/nodeSize); }

/**
 *	@return the number of states contained in one node
 */
int tpat_nodeset::getNodeSize() const { return nodeSize; }

/**
 *	@return the type of node distribution employed by this nodeset
 */
tpat_nodeset::node_distro_t tpat_nodeset::getNodeDistro() const { return nodeDistro; }

/**
 *	@return the indices of nodes that are continuous in velocity with arcs beforehand.
 */
std::vector<int> tpat_nodeset::getVelConNodes() {
	// If the user hasn't set which nodes are continuous, automatically
	// initialize the list so ALL nodes are continuous
	if(!velConSet){
		velConNodes.reserve(getNumNodes());
		// The first node isn't continuous since it has no preceding arc
		for(int n = 1; n < getNumNodes(); n++){
			velConNodes.push_back(n);
		}
		velConSet = true;
	}
	return velConNodes; 
}//=========================================

/**
 *	@brief Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
tpat_constraint tpat_nodeset::getConstraint(int i) const{
	tpat_constraint temp(constraints.at(i)); 
	return temp;
}

/**
 *	@return the number of constraints stored in this nodeset
 */
int tpat_nodeset::getNumCons() const{ return constraints.size(); }

/**
 *	@brief Add a constraint to the nodeset
 *	@param c a constraint to add
 */
void tpat_nodeset::addConstraint(tpat_constraint c){ constraints.push_back(c); }

/**
 *	@brief Append a new node to the end of the nodes vector
 *	@param node a new node; should have <tt>nodeSize</tt> non-dimensional states
 */
void tpat_nodeset::appendNode(std::vector<double> node){
	nodes.insert(nodes.end(), node.begin(), node.end());
}

/**
 *	@brief Append a new node to the end of the nodes vector
 *	@param node an array that MUST have <tt>nodeSize</tt> number
 *	of elements. Fewer elements will result in reading past the
 *	end of the array, more will simply be ignored.
 */	
void tpat_nodeset::appendNode(double *node){
	nodes.insert(nodes.end(), node, node+nodeSize);
}

/**
 *	@brief Append a new TOF to the end of the TOF vector
 *	@param tof a new TOF, non-dimensional units
 */
void tpat_nodeset::appendTOF(double tof){
	tofs.push_back(tof);
}

/**
 *	@brief Set the node distribution type; this will not perform any calculations, it is for 
 *	informational purposes only.
 *	@param type the type of node distribution
 */
void tpat_nodeset::setNodeDistro(node_distro_t type){ nodeDistro = type; }

/**
 *	@brief Set the specified nodes to be continuous in velocity
 *	@param n indices of the nodes that are continuous in velocity
 */
void tpat_nodeset::setVelConNodes(std::vector<int> n){
	velConNodes = n;
	velConSet = true;
}

/**
 *	@brief Set all nodes to be continuous in velocity except for those specified
 *	@param notCont the indices of the nodes that are NOT continuous in velocity
 */
void tpat_nodeset::setVelConNodes_allBut(std::vector<int> notCont){
	velConNodes.clear();
	// Don't ever include node 0 (there is no segment with which it should be continuous)
	for(int n = 1; n < ((int)nodes.size()/nodeSize); n++){
		if(std::find(notCont.begin(), notCont.end(), n) == notCont.end()){
			velConNodes.push_back(n);
		}
	}
	velConSet = true;
}//================================================

/**
 *	Clear all constraints from the constraint vector; useful to remove default constraints
 */
void tpat_nodeset::clearConstraints(){
	constraints.clear();
}//=====================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param sysData a pointer to a system data object describing the system the nodeset will exist in
 *	@param t0 time that corresponds to IC, non-dimensional
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
void tpat_nodeset::initSetFromICs(double IC[6], tpat_sys_data *sysData, double t0, double tof, int numNodes, 
		node_distro_t type){
	// Set up the simulation engine
	tpat_simulation_engine engine(sysData);
	engine.setVerbose(true);
	engine.clearEvents();	// Don't use default crash events to avoid infinite loop

	switch(type){
		default:
		case tpat_nodeset::NONE:
			printWarn("Nodeset type is NONE or not specified, using TIME\n");
		case tpat_nodeset::TIME:
			engine.setVarStepSize(false);
			engine.setNumSteps(numNodes);
			nodeDistro = tpat_nodeset::TIME;
			break;
		case tpat_nodeset::ARCLENGTH:
			engine.setVarStepSize(true);
			engine.setNumSteps(abs(tof)*500);
			nodeDistro = tpat_nodeset::ARCLENGTH;
			break;
	}
	engine.setRevTime(tof < 0);
	printColor(BLUE, "    Set RevTime to %s\n", engine.usesRevTime() ? "ON" : "OFF");

	// Run the simulation and get the trajectory
	engine.runSim(IC, t0, tof);
	tpat_trajectory traj = engine.getTraj();

	printColor(BLUE, "    Retrieved trajectory\n");

	// Reserve space in the nodes and tof vectors
	nodes.reserve(numNodes*nodeSize);
	tofs.reserve(numNodes);

	// Save the nodes and TOFs
	vector<double> *trajState = traj.getState();
	vector<double> *trajTime = traj.getTime();

	double sumArclen = 0;
	double desiredArclen = 0;
	const int k = tpat_trajectory::STATE_WIDTH;

	printColor(BLUE, "    Moving nodes from sim into nodeset\n");
	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == tpat_nodeset::ARCLENGTH){
			if(n > 0){
				// Compute the total length of the trajectory (approx.)
				double dx = trajState->at(n*k) - trajState->at((n-1)*k);
				double dy = trajState->at(n*k+1) - trajState->at((n-1)*k+1);
				double dz = trajState->at(n*k+2) - trajState->at((n-1)*k+2);
				sumArclen += sqrt(dx*dx + dy*dy + dz*dz);
			}
		}else{
			// if TIME is the type, every state on the trajectory is a node, so just copy them over
			nodes.insert(nodes.end(), trajState->begin()+n*k, trajState->begin()+n*k+6);
		}
	}

	desiredArclen = sumArclen/(numNodes-1);
	sumArclen = 0;

	printColor(BLUE, "    Computing TOFs\n");
	int lastNode = 0;
	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == tpat_nodeset::ARCLENGTH){
			if(n == 0){
				nodes.insert(nodes.end(), trajState->begin(), trajState->begin()+6);
			}else{
				// Use linear approximation to compute distance between points
				double dx = trajState->at(n*k) - trajState->at((n-1)*k);
				double dy = trajState->at(n*k+1) - trajState->at((n-1)*k+1);
				double dz = trajState->at(n*k+2) - trajState->at((n-1)*k+2);
				sumArclen += sqrt(dx*dx + dy*dy + dz*dz);

				if(sumArclen > desiredArclen){
					// Save this state as a node
					nodes.insert(nodes.end(), trajState->begin()+n*k, trajState->begin()+ n*k + 6);
					tofs.push_back(trajTime->at(n) - trajTime->at(lastNode));
					lastNode = n;	// update
					sumArclen = 0;	// reset
				}
			}
		}else{
			// Compute the times-of-flight between nodes
			if(n > 0)
				tofs.push_back(trajTime->at(n) - trajTime->at(n-1));
		}
	}

	// Add the last state if it hasn't been already for ARCLENGTH type
	if(nodeDistro == tpat_nodeset::ARCLENGTH && lastNode < traj.getLength()-1){
		int n = traj.getLength()-1;
		nodes.insert(nodes.end(), trajState->begin()+n*k, trajState->begin()+n*k + 6);
		tofs.push_back(trajTime->at(n) - trajTime->at(lastNode));
	}
}//==========================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_nodeset::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file @version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		saveNodes(matfp);
		saveTOFs(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the nodes to a file
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_nodeset::saveNodes(mat_t *matFile){

	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	vector<double> posVel(nodes.size());
	int numNodes = getNumNodes();
	for(int r = 0; r < numNodes; r++){
		for(int c = 0; c < nodeSize; c++){
			posVel[c*numNodes + r] = nodes[r*nodeSize + c];
		}
	}

	// Next, create a matlab variable for the nodes and save it to the file
	size_t dims[2] = {static_cast<size_t>(numNodes), 6};
	matvar_t *matvar = Mat_VarCreate("Nodes", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(posVel[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Nodes", MAT_COMPRESSION_NONE);
}//==================================================

/**
 *	@brief Save the tof values to a file
 *	@param matFile a pointer to the destination matlab file
 */
void tpat_nodeset::saveTOFs(mat_t *matFile){
	size_t dims[2] = {tofs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("TOFs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(tofs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "TOFs", MAT_COMPRESSION_NONE);
}
