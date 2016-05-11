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

#include "tpat_nodeset.hpp"

#include "tpat_event.hpp"
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
tpat_nodeset::tpat_nodeset(const tpat_sys_data *sys) : tpat_base_arcset(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a nodeset from another nodeset
 *	@param n a nodeset reference
 */
tpat_nodeset::tpat_nodeset(const tpat_nodeset &n) : tpat_base_arcset (n){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a nodeset from its base object
 *	@param a an arc data object
 */
tpat_nodeset::tpat_nodeset(const tpat_base_arcset &a) : tpat_base_arcset (a){
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
tpat_nodeset::tpat_nodeset(const tpat_nodeset &n, int first, int last) : tpat_base_arcset(n){
	(void) n;
	(void) first;
	(void) last;

	throw tpat_exception("tpat_nodeset::tpat_nodeset: Not yet implemented!\n");
	// steps.clear();

	// if(first < 0 || last > (int)(n.steps.size()))
	// 	throw tpat_exception("tpat_nodeset::tpat_node: node index out of bounds");
	
	// if(last > first)	// Insert a range
	// 	steps.insert(steps.end(), n.steps.begin()+first, n.steps.begin()+last);
	// else	// first = last, so just insert the specified node
	// 	steps.insert(steps.end(), n.steps[first]);
}//====================================================

tpat_nodeset::~tpat_nodeset(){}

/**
 *  @brief Create a new nodeset object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object
 *  @return a pointer to the newly created nodeset
 */
baseArcsetPtr tpat_nodeset::create( const tpat_sys_data *sys) const{
	return baseArcsetPtr(new tpat_nodeset(sys));
}//====================================================

/**
 *  @brief Create a new nodeset object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned nodeset
 */
 baseArcsetPtr tpat_nodeset::clone() const{
	return baseArcsetPtr(new tpat_nodeset(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

tpat_nodeset operator +(const tpat_nodeset &lhs, const tpat_nodeset &rhs){
	const tpat_nodeset lhs_cpy(lhs);
	const tpat_nodeset rhs_cpy(rhs);
	tpat_nodeset result(lhs.sysData);

	tpat_base_arcset::sum(&lhs, &rhs, &result);

	return result;
}//====================================================

tpat_nodeset& tpat_nodeset::operator +=(const tpat_nodeset &rhs){
	tpat_nodeset temp = *this + rhs;
	copyMe(temp);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Insert a node after the specified node at any locations where the
 *  specified event occurs
 *  @details This function <i>does</i> adjust the prior node to ensure that 
 *  times of flights and other parameters will lead to a nearly continuous 
 *  integrated path. The minimum acceptable time between events is assumed to 
 *  be 1e-2 nondimensional units.
 * 
 *  @param priorNodeIx Index of the node prior to this one; new nodes will be 
 *  inserted after this node.
 *  @param evt the event that identifies the node locations. If multiple occurences
 *  are located, multiple nodes will be inserted. If the event does not occur,
 *  no nodes are inserted
 * 
 *  @return the number of nodes created and inserted into the nodeset.
 */
int tpat_nodeset::createNodesAtEvent(int priorNodeIx, tpat_event evt){
	std::vector<tpat_event> events(1, evt);
	return createNodesAtEvents(priorNodeIx, events);
}//====================================================

/**
 *  @brief Insert a node after the specified node at any locations where the
 *  specified events occur
 *  @details This function <i>does</i> adjust the prior node to ensure that 
 *  times of flights and other parameters will lead to a nearly continuous 
 *  integrated path. The minimum acceptable time between events is assumed to 
 *  be 1e-2 nondimensional units.
 * 
 *  @param priorNodeIx Index of the node prior to this one; new nodes will be 
 *  inserted after this node.
 *  @param evts a vector of events that identify the node locations. If multiple occurences
 *  are located, multiple nodes will be inserted. If nond of the events occur,
 *  no nodes are inserted
 
 *  @return the number of nodes created and inserted into the nodeset.
 *  @see createNodesAtEvent
 *  @see createNodesAtEvents(int, std::vector<tpat_event>, double)
 */
int tpat_nodeset::createNodesAtEvents(int priorNodeIx, std::vector<tpat_event> evts){
	return createNodesAtEvents(priorNodeIx, evts, 1e-2);
}//====================================================

/**
 *  @brief Insert nodes on the specified segment at locations where the
 *  specified events occur
 *  
 *  @param segID the ID of the segment to add nodes to; this segment will be deleted and
 *  replaced by a series of smaller segments (and nodes) if one or more of the input
 *  events occurs
 *  @param evts a vector of events that identify the node locations. If multiple occurences
 *  are located, multiple nodes will be inserted. If nond of the events occur,
 *  no nodes are inserted
 *  @param minTimeDiff Minimum time (nondimensional) between nodes; all segments *must* have
 *  times-of-flight greater than or equal to this amount
 *  
 *  @return the number of nodes created and inserted into the nodeset.
 *  @throws tpat_exception if <tt>segID</tt> is out of bounds
 */
int tpat_nodeset::createNodesAtEvents(int segID, std::vector<tpat_event> evts, double minTimeDiff){
	if(segID < 0 || segID >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_nodeset::createNodesAtEvents: Segment ID is out of bounds");

	// Get a copy of the segment we are replacing
	tpat_segment seg = segs[segIDMap[segID]];
	tpat_node origin = nodes[nodeIDMap[seg.getOrigin()]];
	tpat_node terminus = nodes[nodeIDMap[seg.getTerminus()]];

	// Create a simulation engine and add the events to it
	tpat_simulation_engine engine(sysData);
	engine.setRevTime(seg.getTOF() < 0);
	engine.clearEvents();
	for(size_t i = 0; i < evts.size(); i++){
		evts[i].setStopOnEvent(false);		// Ignore stopping conditions that other processes may have imposed
		engine.addEvent(evts[i]);
	}

	/*	Get an arc that spans the entire segement less a small amount of time at the end
	 * 	to ensure that no new nodes are created with TOF less than minTimeDiff
	 */
	double segTOF = seg.getTOF() < 0 ? seg.getTOF() + std::abs(minTimeDiff) : seg.getTOF() - std::abs(minTimeDiff);
	engine.runSim(origin.getState(), origin.getEpoch(), segTOF);
	tpat_traj traj = engine.getTraj();

	double T0 = traj.getEpochByIx(0);
	std::vector<tpat_event> events = engine.getEvents();
	std::vector<eventRecord> evtRecs = engine.getEventRecords();
	int evtCount = 0, prevNodeID = origin.getID();
	double tof = 0;
	for(size_t e = 0; e < evtRecs.size(); e++){
		for(size_t i = 0; i < evts.size(); i++){

			// If the event occurred, find the corresponding trajectory state and add that to the nodeset
			if(events[evtRecs[e].eventIx] == evts[i]){
				
				// If at least one event is found, we need to delete the segment that is being replaced
				if(evtCount == 0)
					deleteSeg(segID);

				int stepIx = evtRecs[e].stepIx;
				tof = traj.getEpochByIx(stepIx) - T0;

				if(tof > minTimeDiff){
					int newID = addNode(tpat_node(traj.getStateByIx(stepIx), traj.getEpochByIx(stepIx)));
					addSeg(tpat_segment(prevNodeID, newID, tof));

					prevNodeID = newID;
					T0 += tof;
					evtCount++;
				}
			}
		}
	}

	if(evtCount > 0){
		// Add a final segment connecting the last node to the original terminus
		tof = terminus.getEpoch() - nodes[nodeIDMap[prevNodeID]].getEpoch();
		addSeg(tpat_segment(prevNodeID, terminus.getID(), tof));
	}

	return evtCount;
}//====================================================

/**
 *	@brief Allow velocity discontinuities (i.e., delta-Vs) at the specified nodes
 *	@param ix a vector of node indices that can have velocity discontinuities
 */
void tpat_nodeset::allowDV_at(std::vector<int> ix) {
	for(size_t i = 0; i < nodes.size(); i++){
		// Check to see if the node should have continuous velocity
		if(std::find(ix.begin(), ix.end(), i) == ix.end()){
			nodes[i].setVel_AllCon();
		}else{
			nodes[i].setVel_AllDiscon();
		}
	}
}//====================================================

/**
 *  @brief Allow velocity discontinuities (i.e., delta-Vs) at all nodes
 */
void tpat_nodeset::allowDV_all(){
	for(size_t i = 0; i < nodes.size()-1; i++){
		nodes[i].setVel_AllDiscon();
	}
}//====================================================

/**
 *  @brief Allow velocity discontinuities (i.e., delta-Vs) at none of the nodes
 */
void tpat_nodeset::allowDV_none(){
	for(size_t i = 0; i < nodes.size()-1; i++){
		nodes[i].setVel_AllCon();
	}
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Display a textual representation of this object in the standard output
 */
void tpat_nodeset::print() const{
	printf("%s Nodeset:\n Nodes: %zu\n Segments: %zu\n", sysData->getTypeStr().c_str(),
		nodes.size(), segs.size());
	printf("List of Nodes:\n");
	for(size_t n = 0; n < nodeIDMap.size(); n++){
		printf("  %02lu (ix %02d):", n, nodeIDMap[n]);
		if(nodeIDMap[n] != tpat_linkable::INVALID_ID){
			std::vector<double> state = nodes[nodeIDMap[n]].getState();
			printf(" @ %13.8f -- {%13.8f, %13.8f, %13.8f, %13.8f, %13.8f, %13.8f}\n",
				nodes[nodeIDMap[n]].getEpoch(), state[0], state[1], state[2], state[3],
				state[4], state[5]);
		}else{
			printf(" [N/A]\n");
		}
	}

	printf("List of Segments:\n");
	for (size_t s = 0; s < segIDMap.size(); s++){
		printf("  %02lu (ix %02d):", s, segIDMap[s]);
		if(segIDMap[s] != tpat_linkable::INVALID_ID && segIDMap[s] < (int)(segs.size())){
			printf(" origin @ %02d, terminus @ %02d, TOF = %13.8f\n", segs[segIDMap[s]].getOrigin(),
				segs[segIDMap[s]].getTerminus(), segs[segIDMap[s]].getTOF());
		}else{
			printf(" [N/A]\n");
		}
	}

	printf(" Constraints:\n");
	for(size_t n = 0; n < nodes.size(); n++){
		std::vector<tpat_constraint> nodeCons = nodes[n].getConstraints();
		for(size_t c = 0; c < nodeCons.size(); c++){
			nodeCons[c].print();
		}
	}
	for(size_t s = 0; s < segs.size(); s++){
		std::vector<tpat_constraint> segCons = segs[s].getConstraints();
		for(size_t c = 0; c < segCons.size(); c++){
			segCons[c].print();
		}
	}
	for(size_t c = 0; c < cons.size(); c++){
		cons[c].print();
	}

	printf(" Velocity Discontinuities allowed at nodes ");
	char velEl[] = {'x', 'y', 'z'};
	bool anyDiscon = false;
	for(size_t n = 0; n < nodes.size(); n++){
		std::vector<bool> velCon = nodes[n].getVelCon();
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
	for(size_t s = 0; s < segs.size(); s++){
		int o = segs[s].getOrigin();
		segs[s].setOrigin(segs[s].getTerminus());
		segs[s].setTerminus(o);
		segs[s].setTOF(segs[s].getTOF()*-1);
	}

	// Make first node discontinuous in velocity (no preceding node)
	nodes[0].setVel_AllDiscon();
	// and the last node continuous (no previous info about its continuity)
	nodes[nodes.size()-1].setVel_AllCon();
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_nodeset::saveToMat(const char* filename) const{
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
		saveEpoch(matfp, "Epochs");
		saveTOF(matfp, "TOFs");
		sysData->saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *  @brief Populate data in this nodeset from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 *  @throws tpat_exception if the file cannot be loaded
 */
void tpat_nodeset::readFromMat(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw tpat_exception("tpat_nodeset: Could not load data from file");
	}

	initNodesSegsFromMat(matfp, "Nodes");	// This function MUST be called before other data reading functions
	readStateFromMat(matfp, "Nodes");
	readEpochFromMat(matfp, "Epochs");
	readTOFFromMat(matfp, "TOFs");

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions
 *	@param IC a set of initial conditions, non-dimensional units associated with the 
 *	@param t0 time that corresponds to IC, non-dimensional
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC (must be at least 2)
 *	@param distroType node distribution type
 *	@throws tpat_exception if <tt>numNodes</tt> is less than two
 */
void tpat_nodeset::initFromICs(const double IC[6], double t0, double tof, int numNodes, tpat_nodeDistro_tp distroType){

	if(numNodes < 2){
		throw tpat_exception("tpat_nodeset::initFromICs: Nodeset must have at least two nodes!");
	}

	// Prepare to add nodes
	nodes.reserve(numNodes);
	segs.reserve(numNodes-1);

	switch(distroType){
		default:
		case tpat_nodeset::DISTRO_NONE:
			printWarn("Nodeset type is NONE or not specified, using DISTRO_TIME\n");
		case tpat_nodeset::DISTRO_TIME:
			// Initialize using time
			initFromICs_time(IC, t0, tof, numNodes);
			break;
		case tpat_nodeset::DISTRO_ARCLENGTH:
			// Initialize using arclength
			initFromICs_arclength(IC, t0, tof, numNodes);
			break;
	}
}//==========================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions; discretize the arc such that
 *	each segment has approximately the same time-of-flight.
 *	@param IC a set of initial conditions, non-dimensional
 *	@param t0 time that corresponds to IC, non-dimensional
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
void tpat_nodeset::initFromICs_time(const double IC[6], double t0, double tof, int numNodes){
	tpat_simulation_engine engine(sysData);
	engine.setVerbose(SOME_MSG);
	engine.clearEvents();	// Don't use default crash events to avoid infinite loop
	engine.setRevTime(tof < 0);

	int id = addNode(tpat_node(IC, t0));

	double segTOF = tof/(numNodes-1);
	std::vector<double> ic(IC, IC+6);
	for(int n = 0; n < numNodes-1; n++){
		engine.runSim(ic, t0 + n*segTOF, segTOF);
		tpat_traj traj = engine.getTraj();

		id = addNode(tpat_node(traj.getStateByIx(-1), traj.getTimeByIx(-1)));
		addSeg(tpat_segment(id-1, id, segTOF));

		ic = traj.getStateByIx(-1);
	}
}//==========================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions; discretize the arc such that
 *	each segment has approximately the same arclength in distance.
 *	@param IC a set of initial conditions, non-dimensional units associated with the 
 *	@param t0 time that corresponds to IC, non-dimensional
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
void tpat_nodeset::initFromICs_arclength(const double IC[6], double t0, double tof, int numNodes){
	tpat_simulation_engine engine(sysData);
	engine.setVerbose(SOME_MSG);
	engine.clearEvents();	// Don't use default crash events to avoid infinite loop
	engine.setRevTime(tof < 0);

	// Run the simulation and get the trajectory
	engine.runSim(IC, t0, tof);
	tpat_traj traj = engine.getTraj();

	// Compute the total arc length using a linear approximation
	double sumArclen = 0;
	std::vector<double> allArcLen(traj.getNumSegs(), 0);
	std::vector<double> allTOF(traj.getNumSegs(), 0);

	for (int n = 1; n < traj.getNumNodes(); n++){
		std::vector<double> state = traj.getStateByIx(n);
		std::vector<double> prevState = traj.getStateByIx(n-1);

		// Compute the total length of the trajectory (approx.)
		double dx = state[0] - prevState[0];
		double dy = state[1] - prevState[1];
		double dz = state[2] - prevState[2];
		double d = sqrt(dx*dx + dy*dy + dz*dz);
		
		allTOF[n-1] = traj.getEpochByIx(n) - traj.getEpochByIx(n-1);
		allArcLen[n-1] = d;
		sumArclen += d;
	}
	double desiredArclen = sumArclen/(numNodes-1);

	// Add one node to the set
	int prevID = addNode(tpat_node(traj.getStateByIx(0), traj.getEpochByIx(0)));

	// Loop through again to create trajectory
	sumArclen = 0;
	double sumTOF = 0;
	for(size_t s = 0; s < allArcLen.size(); s++){
		
		// Keep adding arclength between steps until the desired length is reached
		if(sumArclen < desiredArclen){
			sumArclen += allArcLen[s];
			sumTOF += allTOF[s];
		}else{
			// reached desired length: save node
			int id = addNode(tpat_node(traj.getStateByIx(s), traj.getEpochByIx(s)));
			
			// Add a segment to link the previous node and this new one
			addSeg(tpat_segment(prevID, id, sumTOF));

			// Reset counters and index variables
			sumArclen = 0;
			sumTOF = 0;
			prevID = id;
		}
	}

	// Save the final state as the last node
	int id2 = addNode(tpat_node(traj.getStateByIx(-1), traj.getEpochByIx(-1)));
	addSeg(tpat_segment(prevID, id2, sumTOF));
}//==========================================================

/**
 *	@brief Split a trajectory into nodes using the specified distribution type
 *
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory. Additionally,
 *	if the input trajectory contains discontinuities in position, velocity, or time, this 
 *	algorithm will fail to capture these discontinuities and could result in a very different
 *	path.
 *
 *	@param traj a trajectory to make into a nodeset
 *	@param numNodes the number of nodes to create, including IC
 *	@param type the node distribution type
 */
void tpat_nodeset::initFromTraj(tpat_traj traj, int numNodes, tpat_nodeDistro_tp type){
	/* Could I code this more intelligently? Probably. Am I too lazy? Definitely */ 
	double ic[] = {0,0,0,0,0,0};
	std::vector<double> trajIC = traj.getStateByIx(0);
	std::copy(trajIC.begin(), trajIC.begin()+6, ic);
	
	initFromICs(ic, traj.getEpochByIx(0), traj.getEpochByIx(-1) - traj.getEpochByIx(0), numNodes, type);
}//==============================================

/**
 *	@brief Initialize the extraParam vector to hold nodeset-specific data
 */
void tpat_nodeset::initExtraParam(){
	// Nothing to do right now!
}//==============================================

