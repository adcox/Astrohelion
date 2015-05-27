/**
 *	@file adtk_nodeset.cpp
 */

#include "adtk_nodeset.hpp"

#include "adtk_simulation_engine.hpp"
#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"

#include <cstdio>
#include <iostream>
#include <cmath>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_nodeset::adtk_nodeset(const int n) : nodeSize(n){
	// Reserve space for one node and TOF
	nodes.reserve(n);
	tofs.reserve(1);
}//========================================================

/**
 *	Copy constructor
 */
adtk_nodeset::adtk_nodeset(const adtk_nodeset& n) : nodeSize(n.nodeSize){
	if(nodeSize == n.nodeSize){
		nodes = n.nodes;
		tofs = n.tofs;
		velConNodes = n.velConNodes;
	}else{
		fprintf(stderr, "Cannot create a nodeset with %d-size nodes from one with %d-size nodes\n",
			nodeSize, n.nodeSize);
		throw;
	}
}//========================================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

adtk_nodeset& adtk_nodeset::operator =(const adtk_nodeset &n){
	if(nodeSize == n.nodeSize){
		nodes = n.nodes;
		tofs = n.tofs;
		velConNodes = n.velConNodes;
		return *this;
	}else{
		fprintf(stderr, "Cannot create a nodeset with %d-size nodes from one with %d-size nodes\n",
			nodeSize, n.nodeSize);
		throw;
	}
}//======================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

std::vector<double>* adtk_nodeset::getNodes(){ return &nodes; }
std::vector<double>* adtk_nodeset::getTOFs(){ return &tofs; }

/**
 *	@param i the index of the node (begins at zero)
 *	@return the node specified by <tt>i</tt>
 */
std::vector<double> adtk_nodeset::getNode(int i) const {
	vector<double> temp(nodes.begin()+i*nodeSize, nodes.begin()+(i+1)*nodeSize);
	return temp;
}

/**
 *	@oaram i the index of the TOF; note that if there are n nodes, there will be n-1 TOFs.
 *	@return the TOF between nodes <tt>i</tt> and <tt>i+1</tt>
 */
double adtk_nodeset::getTOF(int i) const { return tofs.at(i); }

/**
 *	@return the number of nodes contained in this node set
 */
int adtk_nodeset::getNumNodes() const { return (int)(nodes.size()/nodeSize); }

/**
 *	@return the number of states contained in one node
 */
int adtk_nodeset::getNodeSize() const { return nodeSize; }

/**
 *	@return the type of node distribution employed by this nodeset
 */
adtk_nodeset::node_distro_t adtk_nodeset::getNodeDistro() const { return nodeDistro; }

/**
 *	@return the indices of nodes that are continuous in velocity with arcs beforehand.
 */
std::vector<int> adtk_nodeset::getVelConNodes() {
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
 *	Append a new node to the end of the nodes vector
 *	@param node a new node; should have <tt>nodeSize</tt> non-dimensional states
 */
void adtk_nodeset::appendNode(std::vector<double> node){
	nodes.insert(nodes.end(), node.begin(), node.end());
}

/**
 *	Append a new TOF to the end of the TOF vector
 *	@param tof a new TOF, non-dimensional units
 */
void adtk_nodeset::appendTOF(double tof){
	tofs.push_back(tof);
}

/**
 *	Set the node distribution type; this will not perform any calculations, it is for 
 *	informational purposes only.
 *	@param type the type of node distribution
 */
void adtk_nodeset::setNodeDistro(node_distro_t type){ nodeDistro = type; }

/**
 *	Set the specified nodes to be continuous in velocity
 *	@param n indices of the nodes that are continuous in velocity
 */
void adtk_nodeset::setVelConNodes(std::vector<int> n){
	velConNodes = n;
	velConSet = true;
}

/**
 *	Set all nodes to be continuous in velocity except for those specified
 *	@param notCont the indices of the nodes that are NOT continuous in velocity
 */
void adtk_nodeset::setVelConNodes_allBut(std::vector<int> notCont){
	velConNodes.clear();
	// Don't ever include node 0 (there is no segment with which it should be continuous)
	for(int n = 1; n < ((int)nodes.size()/nodeSize); n++){
		if(std::find(notCont.begin(), notCont.end(), n) == notCont.end()){
			velConNodes.push_back(n);
		}
	}
	velConSet = true;
}//================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param t0 time that corresponds to IC, non-dimensional
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
void adtk_nodeset::initSetFromICs(double IC[6], adtk_sys_data *sysData, double t0, double tof, int numNodes, 
		node_distro_t type){
	// Set up the simulation engine
	adtk_simulation_engine engine(sysData);
	engine.setVerbose(true);
	switch(type){
		default:
		case adtk_nodeset::NONE:
			cout << "NONE or not specified, using TIME" << endl;
		case adtk_nodeset::TIME:
			engine.setVarStepSize(false);
			engine.setNumSteps(numNodes);
			nodeDistro = adtk_nodeset::TIME;
			break;
		case adtk_nodeset::ARCLENGTH:
			engine.setVarStepSize(true);
			engine.setNumSteps(abs(tof)*500);
			nodeDistro = adtk_nodeset::ARCLENGTH;
			break;
	}
	engine.setRevTime(tof < 0);

	// Run the simulation and get the trajectory
	engine.runSim(IC, t0, tof);
	adtk_trajectory traj = engine.getTraj();

	// Reserve space in the nodes and tof vectors
	nodes.reserve(numNodes*nodeSize);
	tofs.reserve(numNodes);

	// Save the nodes and TOFs
	vector<double> *trajState = traj.getState();
	vector<double> *trajTime = traj.getTime();

	double sumArclen = 0;
	double desiredArclen = 0;
	const int k = adtk_trajectory::STATE_WIDTH;

	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == adtk_nodeset::ARCLENGTH){
			// Compute the total length of the trajectory (approx.)
			double dx = trajState->at(n*k) - trajState->at((n-1)*k);
			double dy = trajState->at(n*k+1) - trajState->at((n-1)*k+1);
			double dz = trajState->at(n*k+2) - trajState->at((n-1)*k+2);
			sumArclen += sqrt(dx*dx + dy*dy + dz*dz);
		}else{
			// if TIME is the type, every state on the trajectory is a node, so just copy them over
			nodes.insert(nodes.end(), trajState->begin()+n*k, trajState->begin()+n*k+6);
		}
	}

	desiredArclen = sumArclen/(numNodes-1);
	sumArclen = 0;

	int lastNode = 0;
	for (int n = 0; n < traj.getLength(); n++){
		if(nodeDistro == adtk_nodeset::ARCLENGTH){
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
}//==========================================================
