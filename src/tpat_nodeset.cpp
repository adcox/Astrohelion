/**
 *  @file tpat_nodeset.cpp
 *	@brief Contains a set of nodes
 *
 *	@author Andrew Cox
 *	@version August 30, 2015
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_exceptions.hpp"
#include "tpat_node.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_traj.hpp"
#include "tpat_utilities.hpp"

#include <cmath>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Construct a nodeset for the specified system
 *	@param sys a pointer to a system data object
 */
tpat_nodeset::tpat_nodeset(tpat_sys_data *sys) : tpat_arc_data(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a nodeset from another nodeset
 *	@param n a nodeset reference
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset &n) : tpat_arc_data (n){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a nodeset from its base object
 *	@param a an arc data object
 */
tpat_nodeset::tpat_nodeset(const tpat_arc_data &a) : tpat_arc_data (a){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param n Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset. If
 *	last is the same index as first, only one node (with index = first = last)
 *	will be put in the new nodeset
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset &n, int first, int last) : tpat_arc_data(n){
	steps.clear();

	if(first < 0 || last > (int)(n.steps.size()))
		throw tpat_exception("tpat_nodeset::tpat_node: node index out of bounds");
	
	if(last > first)	// Insert a range
		steps.insert(steps.end(), n.steps.begin()+first, n.steps.begin()+last);
	else	// first = last, so just insert the specified node
		steps.insert(steps.end(), n.steps[first]);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the constraints for a specific node
 *	@param nodeIx the index of the node. If less than 0, it will
 *	count backward from the end of the nodeset
 *	@return a vector of all constraints applied to the specified node
 */
std::vector<tpat_constraint> tpat_nodeset::getNodeCons(int nodeIx) const{
	if(nodeIx < 0)
		nodeIx += steps.size();

	return steps[nodeIx].getConstraints();
}//====================================================

/**
 *	@brief Retrieve a specific node
 *	@param ix the index of the node. If less than 0, the index
 *	will count backwards from the end of the nodeset
 *	@return the requested node
 */
tpat_node tpat_nodeset::getNode(int ix) const {
	if(ix < 0)
		ix += steps.size();

	return tpat_node(steps[ix]);
}//====================================================

/**
 *	@brief Retrieve the number of constraints for the entire nodeset
 *	@return the number of constraints for the entire nodeset
 */
int tpat_nodeset::getNumCons() const { 
	int count = 0;
	for(size_t i = 0; i < steps.size(); i++){
		count += steps[i].getConstraints().size();
	}
	return count;
}//====================================================

/**
 *	@brief Retrieve the number of nodes in the nodeset
 *	@return the number of nodes in the nodeset
 */
int tpat_nodeset::getNumNodes() const { return steps.size(); }

/**
 *	@brief Get the time-of-flight for a specific node
 *	@param ix node index; if less than 0, the index counts
 *	backwards from the end of the nodeset
 *	@return non-dimensional time-of-flight between nodes ix and ix+1
 */
double tpat_nodeset::getTOF(int ix) const {
	if(ix < 0)
		ix += steps.size();

	tpat_node node(steps[ix]);
	return node.getTOF();
}//====================================================

/**
 *	@brief Get the total time-of-flight along this nodeset
 *	@return the total time-of-flight along this nodeset (non-dimensional)
 */
double tpat_nodeset::getTotalTOF() const {
	double total = 0;
	for(size_t ix = 0; ix < steps.size(); ix++){
		tpat_node node(steps[ix]);
		
		// If there are any NAN values, don't add them
		if(node.getTOF() == node.getTOF())
			total += node.getTOF();
	}

	return total;
}//====================================================

/**
 *	@brief Add a constraint to the set
 *	
 *	The constraint object specifies the index of the node it will
 *	be added to.
 *
 *	@param con the constraint to add
 */
void tpat_nodeset::addConstraint(tpat_constraint con){
	if(con.getNode() >= 0 && con.getNode() < (int)(steps.size())){
		steps[con.getNode()].addConstraint(con);
	}else{
		throw tpat_exception("tpat_nodeset::addConstraint: constraint node out of range");
	}
}//====================================================

/**
 *	@brief Append a node to the end of the nodeset
 *	@param node a new node
 */
void tpat_nodeset::appendNode(tpat_node node){
	steps.push_back(node);
}//====================================================

/**
 *	@brief Delete the specified node from the nodeset
 *	@param ix node index; if < 0, index counts backwards from
 *	end of nodeset
 */
void tpat_nodeset::deleteNode(int ix){
	if(ix < 0)
		ix += steps.size();
	steps.erase(steps.begin()+ix);
	updateCons();
}//====================================================

/**
 *	@brief Insert a node at the desired index
 *	@param ix node index; if < 0, will count backwards from end of nodeset
 *	@param node new node
 */
void tpat_nodeset::insertNode(int ix, tpat_node node) {
	if(ix < 0)
		ix += steps.size();

	steps.insert(steps.begin() + ix, node);
	updateCons();
}//====================================================

/**
 *	@brief Set all nodes to be continuous in velocity except for those
 *	specified
 *	@param ix a vector of node indices to make discontinuous in 
 *	velocity
 */
void tpat_nodeset::setVelConNodes_allBut(std::vector<int> ix) {
	for(size_t i = 0; i < steps.size()-1; i++){
		tpat_node *np = static_cast<tpat_node*>(&(steps[i]));
		// Check to see if the node should have continuous velocity
		if(std::find(ix.begin(), ix.end(), i) == ix.end()){
			np->setVel_AllCon();
		}else{
			np->setVel_AllDiscon();
		}
	}
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Remove all constraints from all nodes
 */
void tpat_nodeset::clearConstraints() {
	for(size_t i = 0; i < steps.size(); i++)
		steps[i].clearConstraints();
}//====================================================

/**
 *	@brief Display a textual representation of this object in the standard output
 */
void tpat_nodeset::print() const{
	printf("%s Nodeset:\n Nodes: %zu\n", sysData->getTypeStr().c_str(), steps.size());
	for (size_t n = 0; n < steps.size(); n++){
		std::vector<double> node = steps[n].getPosVelState();
		printf("  %02lu: %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f", n+1,
			node.at(0), node.at(1), node.at(2), node.at(3), node.at(4), node.at(5));
		if(n < steps.size()-1){
			printf("   TOF = %.8f\n", getTOF(n));
		}else{
			printf("\n");
		}
	}
	printf(" Constraints:\n");
	for(size_t n = 0; n < steps.size(); n++){
		std::vector<tpat_constraint> nodeCons = getNodeCons(n);
		for(size_t c = 0; c < nodeCons.size(); c++){
			nodeCons[c].print();
		}
	}

	printf(" Velocity Discontinuities allowed at nodes ");
	char velEl[] = {'x', 'y', 'z'};
	bool anyDiscon = false;
	for(size_t n = 0; n < steps.size(); n++){
		tpat_node node = static_cast<tpat_node>(steps[n]);
		std::vector<bool> velCon = node.getVelCon();
		for(size_t i = 0; i < velCon.size(); i++){
			if(!velCon[i]){
				printf("%zuv_%c, ", n, velEl[i]);
				anyDiscon = true;
			}
		}
	}
	if(!anyDiscon)
		printf("None\n");
	else
		printf("\n");
}//====================================================

/**
 *	@brief Reverse the order of the nodes in this nodeset
 *
 *	The constraints are automatically adjusted so they still 
 *	constraint the same states even though the node indices 
 *	have changed.
 */
void tpat_nodeset::reverseOrder() {
	for(int n = 0; n < std::floor(steps.size()/2); n++){
		std::swap(steps[n], steps[steps.size()-n-1]);
	}

	// Shift TOF back one node, change sign
	for(size_t n = 0; n < steps.size()-1; n++){
		tpat_node *node = static_cast<tpat_node*>(&(steps[n]));
		tpat_node *nextNode = static_cast<tpat_node*>(&(steps[n+1]));
		node->setTOF(-1*(nextNode->getTOF()));
	}
	
	// Set TOF on final node to zero
	tpat_node *endNode = static_cast<tpat_node*>(&(steps[steps.size()-1]));
	endNode->setTOF(0);

	// Make first node discontinuous in velocity (no preceding node)
	// and the last node continuous (no previous info about its continuity)
	tpat_node *first = static_cast<tpat_node*>(&steps[0]);
	tpat_node *last = static_cast<tpat_node*>(&steps[steps.size()-1]);
	first->setVel_AllDiscon();
	last->setVel_AllCon();
}//====================================================

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
		saveState(matfp, "Nodes");
		saveTOFs(matfp);
		sysData->saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//====================================================

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
	engine.setVerbose(SOME_MSG);
	engine.clearEvents();	// Don't use default crash events to avoid infinite loop

	node_distro_t nodeDistro = tpat_nodeset::DISTRO_NONE;

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

	// Compute the total arc length using a linear approximation
	double sumArclen = 0;
	double desiredArclen = 0;
	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == tpat_nodeset::DISTRO_ARCLENGTH){
			if(n > 0){
				std::vector<double> state = traj.getState(n);
				std::vector<double> prevState = traj.getState(n-1);

				// Compute the total length of the trajectory (approx.)
				double dx = state[0] - prevState[0];
				double dy = state[1] - prevState[1];
				double dz = state[2] - prevState[2];
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
				std::vector<double> state = traj.getState(n);
				nodeStates.insert(nodeStates.end(), state.begin(), state.end());
			}else{
				// Use linear approximation to compute distance between points
				std::vector<double> state = traj.getState(n);
				std::vector<double> prevState = traj.getState(n-1);

				double dx = state[0] - prevState[0];
				double dy = state[1] - prevState[1];
				double dz = state[2] - prevState[2];
				sumArclen += sqrt(dx*dx + dy*dy + dz*dz);

				if(sumArclen > desiredArclen){
					// Save this state as a node
					nodeStates.insert(nodeStates.end(), state.begin(), state.end());
					allTOFs.push_back(traj.getTime(n) - traj.getTime(lastNode));
					lastNode = n;	// update
					sumArclen = sumArclen - desiredArclen;	// reset, remember remainder
				}
			}
		}else{
			std::vector<double> state = traj.getState(n);
			nodeStates.insert(nodeStates.end(), state.begin(), state.end());

			// Compute the times-of-flight between nodes
			if(n > 0)
				allTOFs.push_back(traj.getTime(n) - traj.getTime(n-1));
		}
	}

	// Add the last state if it hasn't been already for ARCLENGTH type
	if(nodeDistro == tpat_nodeset::DISTRO_ARCLENGTH && lastNode < traj.getLength()-1){
		int n = traj.getLength()-1;
		std::vector<double> state = traj.getState(n);
		nodeStates.insert(nodeStates.end(), state.begin(), state.end());
		allTOFs.push_back(traj.getTime(n) - traj.getTime(lastNode));
	}

	// Reserve space in the nodes and tof vectors
	steps.reserve(numNodes);

	// Create node objects from saved states and tofs
	for(size_t i = 0; i < nodeStates.size()/6; i++){
		std::vector<double> state;
		state.insert(state.end(), nodeStates.begin()+i*6, nodeStates.begin()+(i+1)*6);
		double tof = 0;
		if(i < allTOFs.size())
			tof = allTOFs[i];

		tpat_node node(state, tof);
		steps.push_back(node);
	}
}//==========================================================

/**
 *	@brief Split a trajectory into nodes using the specified distribution type
 *
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory.
 *
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
 *	@brief Save time-of-flight values to a mat file
 *	@param matFile a pointer to the open mat file
 */
void tpat_nodeset::saveTOFs(mat_t *matFile){
	saveExtraParam(matFile, 0, "TOFs");
}//==============================================

/**
 *	@brief Initialize the extraParam vector to hold nodeset-specific data
 */
void tpat_nodeset::initExtraParam(){
	// Time of Flight for each node/step
	numExtraParam = 1;
	extraParamRowSize.push_back(1);
}//==============================================

