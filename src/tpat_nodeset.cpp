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
#include "tpat_traj.hpp"
#include "tpat_utilities.hpp"
 
#include <algorithm>
#include <cstdio>
#include <cmath>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Default constructor
 */
tpat_nodeset::tpat_nodeset(){}

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset &orig, int first, int last){
	nodeDistro = orig.nodeDistro;

	nodes.insert(nodes.end(), orig.nodes.begin()+first, orig.nodes.begin()+last);

	for(int c = 0; c < ((int)orig.constraints.size()); c++){
		if(orig.constraints[c].getNode() <= last && orig.constraints[c].getNode() >= first){
			constraints.push_back(orig.constraints[c]);
			constraints.back().setNode(orig.constraints[c].getNode() - first);
		}
	}
}//==================================================

/**
 *	@brief Copy constructor
 *	@param n a nodeset
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset &n){
	copyMe(n);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_nodeset::~tpat_nodeset(){
	nodes.clear();
	constraints.clear();
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
	copyMe(n);
	return *this;
}//======================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Get the total TOF for the whole nodeset
 *	@return the total TOF for the whole nodeset (non-dim)
 */
double tpat_nodeset::getTotalTOF() const { 
	double total = 0;
	for(size_t n = 0; n < nodes.size(); n++){
		total += nodes[n].getTOF();
	}
	return total;
}//===================================================

/**
 *	@param i the index of the node (begins at zero). Negative values count
 *	backwards from the end, i.e. -1 gives the final value, -2 gives the second
 *	to last value, etc.
 *	@return the node specified by <tt>i</tt>
 */
tpat_node tpat_nodeset::getNode(int i) const {
	if(i < 0)
		i += nodes.size();

	return nodes.at(i);
}//=========================================

/**
 *	@param i the index of the TOF; note that if there are n nodes, there will be n-1 
 *	TOFs. Negative values count backwards from the end, i.e. -1 gives the final value,
 *	-2 gives the second to last value, etc.
 *	@return the TOF between nodes <tt>i</tt> and <tt>i+1</tt>
 */
double tpat_nodeset::getTOF(int i) const {
	if(i < 0)
		i += nodes.size();
	return nodes[i].getTOF();
}//=========================================

/**
 *	@return the number of nodes contained in this node set
 */
int tpat_nodeset::getNumNodes() const { return nodes.size(); }

/**
 *	@return the type of node distribution employed by this nodeset
 */
tpat_nodeset::node_distro_t tpat_nodeset::getNodeDistro() const { return nodeDistro; }

/**
 *	@brief Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
tpat_constraint tpat_nodeset::getConstraint(int i) const{
	tpat_constraint temp(constraints.at(i)); 
	return temp;
}//=========================================

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
void tpat_nodeset::appendNode(tpat_node node){
	nodes.push_back(node);
}

/**
 *	@brief Remove a node from the nodeset
 *	@param ix the index of the node. If negative, this index will count backwards
 *	from the end of the set
 */
void tpat_nodeset::deleteNode(int ix){
	if(ix < 0)
		ix += nodes.size();

	nodes.erase(nodes.begin()+ix);
}//============================================

/**
 *	@brief Insert a node at the specified position
 *	@param idx the index the new node should have (begins at 0). If the index is negative,
 *	it will count backwards from the end of the vector.
 *	@param newNode a new node to be inserted
 */
void tpat_nodeset::insertNode(int idx, tpat_node newNode){
	if(idx < 0)
		idx += getNumNodes();

	nodes.insert(nodes.begin()+idx, newNode);
}//=========================================

/**
 *	@brief Set the node distribution type; this will not perform any calculations, it is for 
 *	informational purposes only.
 *	@param type the type of node distribution
 */
void tpat_nodeset::setNodeDistro(node_distro_t type){ nodeDistro = type; }

/**
 *	@brief Set all nodes to be continuous in velocity except for those specified
 *	@param notCont the indices of the nodes that are NOT continuous in velocity
 */
void tpat_nodeset::setVelConNodes_allBut(std::vector<int> notCont){
	
	// Don't ever include node 0 (there is no segment with which it should be continuous)
	for(size_t n = 1; n < nodes.size(); n++){
		if(std::find(notCont.begin(), notCont.end(), n) == notCont.end())
			nodes[n].setVel_AllCon();
	}
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
 *	integrated trajectory into pieces (nodes)
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
	engine.setVerbose(false);
	engine.clearEvents();	// Don't use default crash events to avoid infinite loop

	switch(type){
		default:
		case tpat_nodeset::DISTRO_NONE:
			printWarn("Nodeset type is NONE or not specified, using DISTRO_TIME\n");
		case tpat_nodeset::DISTRO_TIME:
			engine.setVarStepSize(false);
			engine.setNumSteps(numNodes);
			nodeDistro = tpat_nodeset::DISTRO_TIME;
			break;
		case tpat_nodeset::DISTRO_ARCLENGTH:
			engine.setVarStepSize(true);
			engine.setNumSteps(std::abs(tof)*500);
			nodeDistro = tpat_nodeset::DISTRO_ARCLENGTH;
			break;
	}
	engine.setRevTime(tof < 0);

	// Run the simulation and get the trajectory
	engine.runSim(IC, t0, tof);
	tpat_traj traj = engine.getTraj();

	// Save the nodes and TOFs
	std::vector<double> *trajState = traj.getState();
	std::vector<double> *trajTime = traj.getTime();

	double sumArclen = 0;
	double desiredArclen = 0;
	const int k = tpat_traj::STATE_SIZE;

	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == tpat_nodeset::DISTRO_ARCLENGTH){
			if(n > 0){
				// Compute the total length of the trajectory (approx.)
				double dx = trajState->at(n*k) - trajState->at((n-1)*k);
				double dy = trajState->at(n*k+1) - trajState->at((n-1)*k+1);
				double dz = trajState->at(n*k+2) - trajState->at((n-1)*k+2);
				sumArclen += sqrt(dx*dx + dy*dy + dz*dz);
			}
		}
	}

	desiredArclen = sumArclen/(numNodes-1);
	sumArclen = 0;

	std::vector<double> nodeStates;
	std::vector<double> allTOFs;

	int lastNode = 0;
	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == tpat_nodeset::DISTRO_ARCLENGTH){
			if(n == 0){
				nodeStates.insert(nodeStates.end(), trajState->begin(), trajState->begin()+6);
			}else{
				// Use linear approximation to compute distance between points
				double dx = trajState->at(n*k) - trajState->at((n-1)*k);
				double dy = trajState->at(n*k+1) - trajState->at((n-1)*k+1);
				double dz = trajState->at(n*k+2) - trajState->at((n-1)*k+2);
				sumArclen += sqrt(dx*dx + dy*dy + dz*dz);

				if(sumArclen > desiredArclen){
					// Save this state as a node
					nodeStates.insert(nodeStates.end(), trajState->begin()+n*k, trajState->begin()+ n*k + 6);
					allTOFs.push_back(trajTime->at(n) - trajTime->at(lastNode));
					lastNode = n;	// update
					sumArclen = 0;	// reset
				}
			}
		}else{
			nodeStates.insert(nodeStates.end(), trajState->begin()+n*k, trajState->begin()+n*k+6);
			// Compute the times-of-flight between nodes
			if(n > 0)
				allTOFs.push_back(trajTime->at(n) - trajTime->at(n-1));
		}
	}

	// Add the last state if it hasn't been already for ARCLENGTH type
	if(nodeDistro == tpat_nodeset::DISTRO_ARCLENGTH && lastNode < traj.getLength()-1){
		int n = traj.getLength()-1;
		nodeStates.insert(nodeStates.end(), trajState->begin()+n*k, trajState->begin()+n*k + 6);
		allTOFs.push_back(trajTime->at(n) - trajTime->at(lastNode));
	}

	// Reserve space in the nodes and tof vectors
	nodes.reserve(numNodes);

	// Create node objects from saved states and tofs
	for(size_t i = 0; i < nodeStates.size()/6; i++){
		std::vector<double> state;
		state.insert(state.end(), nodeStates.begin()+i*6, nodeStates.begin()+(i+1)*6);
		double tof = 0;
		if(i < allTOFs.size())
			tof = allTOFs[i];

		tpat_node node(state, tof);
		nodes.push_back(node);
	}
}//==========================================================

/**
 *	@brief Split a trajectory into nodes using the specified distribution type
 *	@param traj a trajectory to make into a nodeset
 *	@param sysData a pointer to the system data object used to create traj (cannot extract from trajectory base class)
 *	@param numNodes the number of nodes to create, including IC
 *	@param type the node distribution type
 */
void tpat_nodeset::initSetFromTraj(tpat_traj traj, tpat_sys_data *sysData, int numNodes, node_distro_t type){
	/* Could I code this more intelligently? Probably. Am I too lazy? Definitely */ 
	double ic[] = {0,0,0,0,0,0};
	std::vector<double> trajIC = traj.getState(0);
	std::copy(trajIC.begin(), trajIC.begin()+6, ic);
	
	initSetFromICs(ic, sysData, traj.getTime(0), traj.getTime(-1) - traj.getTime(0), numNodes, type);
}//==============================================

/**
 *	@brief Reverse the order of the nodes in this nodeset
 *
 *	The cosntraints are automatically adjusted so they still 
 *	constraint the same states even though the node indices 
 *	have changed.
 */
void tpat_nodeset::reverseOrder(){
	// Re-order nodes and TOFs
	for(int n = 0; n < floor(getNumNodes()/2); n++){
		std::swap(nodes[n], nodes[nodes.size()-n-1]);
	}

	// Reverse sign on TOFs
	for(size_t n = 0; n < nodes.size(); n++){
		nodes[n].setTOF(-1*nodes[n].getTOF());
	}

	// Update node indices for constraints and velConNodes
	for(int c = 0; c < ((int)constraints.size()); c++){
		constraints[c].setNode(nodes.size()-1-constraints[c].getNode());
	}

	// Make the first node discontinuous in velocity (no preceding node)
	nodes[0].setVel_AllDiscon();
	// Make last node continuous (no previous information about its continuity)
	nodes[nodes.size()-1].setVel_AllCon();
}//=====================================

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
		// Append all but the last node
		for(size_t n = 0; n < lhs.nodes.size() - 1; n++){
			output->appendNode(lhs.nodes[n]);
		}
	}else{
		// Append all nodes
		for(size_t n = 0; n < lhs.nodes.size(); n++){
			output->appendNode(lhs.nodes[n]);
		}
	}

	// Append all the RHS nodes
	for(size_t n = 0; n < rhs.nodes.size(); n++){
		output->appendNode(rhs.nodes[n]);
	}

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

	output->setNodeDistro( lhs.nodeDistro == rhs.nodeDistro ? lhs.nodeDistro : DISTRO_NONE);
}//=======================================

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
	std::vector<double> posVel(nodes.size()*6);
	int numNodes = getNumNodes();
	for(int r = 0; r < numNodes; r++){
		std::vector<double> nodeState = nodes[r].getPosVelState();
		for(int c = 0; c < 6; c++){
			posVel[c*numNodes + r] = nodeState[c];
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
	std::vector<double> tofs;
	for(size_t n = 0; n < nodes.size()-1; n++){
		tofs.push_back(nodes[n].getTOF());
	}
	size_t dims[2] = {tofs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("TOFs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(tofs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "TOFs", MAT_COMPRESSION_NONE);
}//====================================================

/**
 *	@brief Copy the nodeset
 *	@param n a nodeset reference
 */
void tpat_nodeset::copyMe(const tpat_nodeset &n){
	nodeDistro = n.nodeDistro;
	nodes = n.nodes;
	constraints = n.constraints;
}//====================================================




