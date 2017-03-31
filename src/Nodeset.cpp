/**
 *  \file Nodeset.cpp
 *	\brief Contains a set of nodes
 *
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Nodeset.hpp"

#include "Event.hpp"
#include "Exceptions.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "Traj.hpp"
#include "Utilities.hpp"

#include <cmath>


namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief Construct a nodeset for the specified system
 *	\param sys a pointer to a system data object
 */
Nodeset::Nodeset(const SysData *sys) : BaseArcset(sys){}

/**
 *	\brief Create a nodeset from another nodeset
 *	\param n a nodeset reference
 */
Nodeset::Nodeset(const Nodeset &n) : BaseArcset (n){}

/**
 *	\brief Create a nodeset from its base object
 *	\param a an arc data object
 */
Nodeset::Nodeset(const BaseArcset &a) : BaseArcset (a){}

/**
 *	\brief Create a nodeset as a subset of another
 *	\param n Original nodeset
 *	\param first index of the first node to be included in the new nodeset
 *	\param last index of the last node to be included in the new nodeset. If
 *	last is the same index as first, only one node (with index = first = last)
 *	will be put in the new nodeset
 */
Nodeset::Nodeset(const Nodeset &n, int first, int last) : BaseArcset(n){
	(void) n;
	(void) first;
	(void) last;

	throw Exception("Nodeset::Nodeset: Not yet implemented!\n");
	// steps.clear();

	// if(first < 0 || last > (int)(n.steps.size()))
	// 	throw Exception("Nodeset::Node: node index out of bounds");
	
	// if(last > first)	// Insert a range
	// 	steps.insert(steps.end(), n.steps.begin()+first, n.steps.begin()+last);
	// else	// first = last, so just insert the specified node
	// 	steps.insert(steps.end(), n.steps[first]);
}//====================================================

/**
 *  \brief Default Destructor
 */
Nodeset::~Nodeset(){}

/**
 *  \brief Create a new nodeset object on the stack
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param sys pointer to a system data object
 *  \return a pointer to the newly created nodeset
 */
baseArcsetPtr Nodeset::create( const SysData *sys) const{
	return baseArcsetPtr(new Nodeset(sys));
}//====================================================

/**
 *  \brief Create a new nodeset object on the stack that is a 
 *  duplicate of this object
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \return a pointer to the newly cloned nodeset
 */
 baseArcsetPtr Nodeset::clone() const{
	return baseArcsetPtr(new Nodeset(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  \brief Combine two nodesets.
 *  \details This function concatenates two nodeset objects. It is assumed
 *  that the first state on <tt>rhs</tt> is identical to the final state on
 *  <tt>rhs</tt>. The <tt>rhs</tt> object is also assumed to occur after
 *  (chronologically) <tt>lhs</tt>
 * 
 *  \param lhs reference to a nodeset object
 *  \param rhs reference to a nodeset object
 * 
 *  \return the concatenation of lhs + rhs.
 */
 Nodeset operator +(const Nodeset &lhs, const Nodeset &rhs){
	const Nodeset lhs_cpy(lhs);
	const Nodeset rhs_cpy(rhs);
	Nodeset result(lhs.pSysData);

	BaseArcset::sum(&lhs, &rhs, &result);

	return result;
}//====================================================

/**
 *  \brief Concatenate this object with another nodeset
 * 
 *  \param rhs reference to a nodeset object
 *  \return the concatenation of this and <tt>rhs</tt>
 *  @see operator +()
 */
Nodeset& Nodeset::operator +=(const Nodeset &rhs){
	Nodeset temp = *this + rhs;
	copyMe(temp);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Insert a node after the specified node at any locations where the
 *  specified event occurs
 *  \details This function <i>does</i> adjust the prior node to ensure that 
 *  times of flights and other parameters will lead to a nearly continuous 
 *  integrated path. A numerical simulation is used to propagate the nonlinear
 *  solution between segments, so the events are defined in the dynamical system
 *  associated with the nodeset.
 * 
 *  \param segID the ID of the segment to add nodes to; this segment will be deleted and
 *  replaced by a series of smaller segments (and nodes) if one or more of the input
 *  events occurs
 *  \param evt the event that identifies the node locations. If multiple occurences
 *  are located, multiple nodes will be inserted. If the event does not occur,
 *  no nodes are inserted
 *  \param minTimeDiff Minimum time (nondimensional) between nodes; all segments *must* have
 *  times-of-flight greater than or equal to this amount (default is 1e-2)
 * 
 *  \return the number of nodes created and inserted into the nodeset.
 */
int Nodeset::createNodesAtEvent(int segID, Event evt, double minTimeDiff){
	std::vector<Event> events(1, evt);
	return createNodesAtEvents(segID, events, minTimeDiff);
}//====================================================

/**
 *  \brief Insert nodes on the specified segment at locations where the
 *  specified events occur.
 *  
 *  \details This function <i>does</i> adjust the prior node to ensure that 
 *  times of flights and other parameters will lead to a nearly continuous 
 *  integrated path. A numerical simulation is used to propagate the nonlinear
 *  solution between segments, so the events are defined in the dynamical system
 *  associated with the nodeset.
 *  
 *  \param segID the ID of the segment to add nodes to; this segment will be deleted and
 *  replaced by a series of smaller segments (and nodes) if one or more of the input
 *  events occurs
 *  \param evts a vector of events that identify the node locations. If multiple occurences
 *  are located, multiple nodes will be inserted. If none of the events occur,
 *  no nodes are inserted
 *  \param minTimeDiff Minimum time (nondimensional) between nodes; all segments *must* have
 *  times-of-flight greater than or equal to this amount (default is 1e-2)
 *  
 *  \return the number of nodes created and inserted into the nodeset.
 *  \throws Exception if <tt>segID</tt> is out of bounds
 */
int Nodeset::createNodesAtEvents(int segID, std::vector<Event> evts, double minTimeDiff){
	if(segIDMap.count(segID) == 0)
		throw Exception("Nodeset::createNodesAtEvents: Segment ID is out of bounds");

	// Get a copy of the segment we are replacing
	Segment seg = segs[segIDMap[segID]];
	Node origin = nodes[nodeIDMap[seg.getOrigin()]];
	Node terminus = nodes[nodeIDMap[seg.getTerminus()]];

	// Create a simulation engine and add the events to it
	SimEngine engine;
	engine.setRevTime(seg.getTOF() < 0);
	engine.setMakeDefaultEvents(false);
	engine.setCtrlLaw(seg.getCtrlLaw());
	for(unsigned int i = 0; i < evts.size(); i++){
		evts[i].setStopOnEvent(false);		// Ignore stopping conditions that other processes may have imposed
		engine.addEvent(evts[i]);
	}

	/*	Get an arc that spans the entire segement less a small amount of time at the end
	 * 	to ensure that no new nodes are created with TOF less than minTimeDiff
	 */
	double segTOF = seg.getTOF() < 0 ? seg.getTOF() + std::abs(minTimeDiff) : seg.getTOF() - std::abs(minTimeDiff);
	Traj traj(pSysData);
	engine.runSim(origin.getState(), origin.getEpoch(), segTOF, &traj);

	double T0 = traj.getEpochByIx(0);
	std::vector<Event> events = engine.getEvents();
	std::vector<SimEventRecord> evtRecs = engine.getEventRecords();
	int evtCount = 0, prevNodeID = origin.getID();
	double tof = 0;
	for(unsigned int e = 0; e < evtRecs.size(); e++){
		for(unsigned int i = 0; i < evts.size(); i++){

			// If the event occurred, find the corresponding trajectory state and add that to the nodeset
			if(events[evtRecs[e].eventIx] == evts[i]){
				
				// If at least one event is found, we need to delete the segment that is being replaced
				if(evtCount == 0)
					deleteSeg(segID);

				int stepIx = evtRecs[e].stepIx;
				tof = traj.getEpochByIx(stepIx) - T0;

				if(tof > minTimeDiff){
					int newID = addNode(Node(traj.getStateByIx(stepIx), traj.getEpochByIx(stepIx)));
					Segment newSeg = Segment(prevNodeID, newID, tof);
					newSeg.setSTM(traj.getSTMByIx(stepIx));
					newSeg.setCtrlLaw(traj.getCtrlLawByIx(0));
					addSeg(newSeg);

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
		Traj temp(pSysData);
		engine.clearEvents();
		engine.runSim(nodes[nodeIDMap[prevNodeID]].getState(), nodes[nodeIDMap[prevNodeID]].getEpoch(), tof, &temp);

		Segment newSeg = Segment(prevNodeID, terminus.getID(), tof);
		newSeg.setSTM(temp.getSTMByIx(-1));
		newSeg.setCtrlLaw(temp.getCtrlLawByIx(0));
		addSeg(newSeg);
	}

	return evtCount;
}//====================================================

/**
 *	\brief Allow velocity discontinuities (i.e., delta-Vs) at the specified segments
 *	\param id a vector of segment IDs that can have velocity discontinuities
 */
void Nodeset::allowDV_at(std::vector<int> id) {
	for(unsigned int i = 0; i < segs.size(); i++){
		// Check to see if the node should have continuous velocity
		if(std::find(id.begin(), id.end(), segs[i].getID()) == id.end()){
			segs[i].setVel_AllCon();
		}else{
			segs[i].setVel_AllDiscon();
		}
	}
}//====================================================

/**
 *  \brief Allow velocity discontinuities (i.e., delta-Vs) on all segments
 */
void Nodeset::allowDV_all(){
	for(unsigned int i = 0; i < segs.size(); i++){
		segs[i].setVel_AllDiscon();
	}
}//====================================================

/**
 *  \brief Allow velocity discontinuities (i.e., delta-Vs) on none of the segments
 */
void Nodeset::allowDV_none(){
	for(unsigned int i = 0; i < segs.size(); i++){
		segs[i].setVel_AllCon();
	}
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	\brief Display a textual representation of this object in the standard output
 */
void Nodeset::print() const{
	printf("%s Nodeset:\n Nodes: %zu\n Segments: %zu\n", pSysData->getTypeStr().c_str(),
		nodes.size(), segs.size());
	printf("List of Nodes:\n");
	for(const auto &index : nodeIDMap){
		printf("  %02d (ix %02d):", index.first, index.second);
		if(index.second != Linkable::INVALID_ID){
			std::vector<double> state = nodes[index.second].getState();
			printf(" @ %13.8f -- {%13.8f, %13.8f, %13.8f, %13.8f, %13.8f, %13.8f}\n",
				nodes[index.second].getEpoch(), state[0], state[1], state[2], state[3],
				state[4], state[5]);
		}else{
			printf(" [N/A]\n");
		}
	}

	printf("List of Segments:\n");
	for(const auto &index : segIDMap){
		printf("  %02d (ix %02d):", index.first, index.second);
		if(index.second != Linkable::INVALID_ID && index.second < static_cast<int>(segs.size())){
			printf(" origin @ %02d, terminus @ %02d, TOF = %13.8f, Ctrl Law ID = %u\n", segs[index.second].getOrigin(),
				segs[index.second].getTerminus(), segs[index.second].getTOF(), segs[index.second].getCtrlLaw());
		}else{
			printf(" [N/A]\n");
		}
	}

	printf(" Constraints:\n");
	for(unsigned int n = 0; n < nodes.size(); n++){
		std::vector<Constraint> nodeCons = nodes[n].getConstraints();
		for(unsigned int c = 0; c < nodeCons.size(); c++){
			nodeCons[c].print();
		}
	}
	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<Constraint> segCons = segs[s].getConstraints();
		for(unsigned int c = 0; c < segCons.size(); c++){
			segCons[c].print();
		}
	}
	for(unsigned int c = 0; c < cons.size(); c++){
		cons[c].print();
	}

	printf(" Velocity Discontinuities allowed on segments: ");
	char velEl[] = {'x', 'y', 'z'};
	bool anyDiscon = false;
	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<bool> velCon = segs[s].getVelCon();
		for(unsigned int i = 0; i < velCon.size(); i++){
			if(!velCon[i]){
				printf("%uv_%c, ", s, velEl[i]);
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
 *	\brief Reverse the order of the nodes in this nodeset
 *
 *	The constraints are automatically adjusted so they still 
 *	constraint the same states even though the node indices 
 *	have changed.
 */
void Nodeset::reverseOrder() {
	for(unsigned int s = 0; s < segs.size(); s++){
		int o = segs[s].getOrigin();
		segs[s].setOrigin(segs[s].getTerminus());
		segs[s].setTerminus(o);
		segs[s].setTOF(segs[s].getTOF()*-1);
	}
}//====================================================

/**
 *	\brief Save the trajectory to a file
 *	\param filename the name of the .mat file
 */
void Nodeset::saveToMat(const char* filename) const{
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file \version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		astrohelion::printErr("Error creating MAT file\n");
	}else{
		saveCmds(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *  \brief Execute commands to save data to a Matlab file
 *  \details This function is called from saveToMat() and should
 *  be overridden in derived classes as necessary.
 * 
 *  \param pMatFile pointer to an open Matlab file
 */
void Nodeset::saveCmds(mat_t *pMatFile) const{
	saveState(pMatFile, VARNAME_NODE);
	saveEpoch(pMatFile);
	saveTOF(pMatFile);
	saveSTMs(pMatFile);
	saveCtrlLaw(pMatFile);
	pSysData->saveToMat(pMatFile);
	// TODO: Add these functions:
	// saveCons(pMatFile);
	// saveVelCon(pMatFile);
}//====================================================

/**
 *  \brief Populate data in this nodeset from a matlab file
 * 
 *  \param filepath the path to the matlab data file
 *  \throws Exception if the file cannot be loaded
 */
void Nodeset::readFromMat(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Nodeset: Could not load data from file");
	}else{
		readCmds(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \details This function is called from readFromMat() and should
 *  be overridden in derived classes as necessary.
 * 
 *  \param pMatFile pointer to an open Matlab file
 */
void Nodeset::readCmds(mat_t *pMatFile){
	initNodesSegsFromMat(pMatFile, VARNAME_NODE);	// This function MUST be called before other data reading functions
	readStateFromMat(pMatFile, VARNAME_NODE);
	readEpochFromMat(pMatFile);
	readTOFFromMat(pMatFile);
	readSTMFromMat(pMatFile);
	readCtrlLawFromMat(pMatFile);
}//====================================================

/**
 *	\brief Compute a set of nodes by integrating from initial conditions
 *	\param IC a set of initial conditions, non-dimensional units
 *	\param t0 time that corresponds to IC, non-dimensional
 *	\param tof duration of the simulation, non-dimensional
 *	\param numNodes number of nodes to create, including IC (must be at least 2)
 *	\param distroType node distribution type
 *	\param ctrlLawID ID of the control law to use while propagating arcs
 *	
 *	\throws Exception if <tt>numNodes</tt> is less than two
 */
void Nodeset::initFromICs(std::vector<double> IC, double t0, double tof, int numNodes, NodeDistro_tp distroType, unsigned int ctrlLawID){

	if(numNodes < 2){
		throw Exception("Nodeset::initFromICs: Nodeset must have at least two nodes!");
	}

	// Prepare to add nodes
	nodes.reserve(numNodes);
	segs.reserve(numNodes-1);

	switch(distroType){
		default:
		case Nodeset::NONE:
			astrohelion::printWarn("Nodeset type is NONE or not specified, using TIME\n");
		case Nodeset::TIME:
			// Initialize using time
			initFromICs_time(IC, t0, tof, numNodes, ctrlLawID);
			break;
		case Nodeset::ARCLENGTH:
			// Initialize using arclength
			initFromICs_arclength(IC, t0, tof, numNodes, ctrlLawID);
			break;
	}
}//==========================================================

/**
 *	\brief Compute a set of nodes by integrating from initial conditions; discretize the arc such that
 *	each segment has approximately the same time-of-flight.
 *	\param IC a set of initial conditions, non-dimensional
 *	\param t0 time that corresponds to IC, non-dimensional
 *	\param tof duration of the simulation, non-dimensional
 *	\param numNodes number of nodes to create, including IC
 *	\param ctrlLawID ID of the control law to use while propagating arcs
 */
void Nodeset::initFromICs_time(std::vector<double> IC, double t0, double tof, int numNodes, unsigned int ctrlLawID){
	SimEngine engine;
	engine.setMakeDefaultEvents(false);	// Don't use default crash events to avoid infinite loop
	engine.setRevTime(tof < 0);
	engine.setCtrlLaw(ctrlLawID);

	double segTOF = tof/(numNodes-1);
	int id = addNode(Node(IC, t0));
	
	for(int n = 0; n < numNodes-1; n++){
		Traj traj(pSysData);
		engine.runSim(IC, t0 + n*segTOF, segTOF, &traj);

		id = addNode(Node(traj.getStateByIx(-1), traj.getTimeByIx(-1)));
		Segment newSeg = Segment(id-1, id, segTOF);
		newSeg.setSTM(traj.getSTMByIx(-1));
		newSeg.setCtrlLaw(traj.getCtrlLawByIx(0));
		addSeg(newSeg);

		IC = traj.getStateByIx(-1);
	}
}//==========================================================

/**
 *	\brief Compute a set of nodes by integrating from initial conditions; discretize the arc such that
 *	each segment has approximately the same arclength in distance.
 *	\param IC a set of initial conditions, non-dimensional units associated with the 
 *	\param t0 time that corresponds to IC, non-dimensional
 *	\param tof duration of the simulation, non-dimensional
 *	\param numNodes number of nodes to create, including IC
 *	\param ctrlLawID ID of the control law to use while propagating arcs
 */
void Nodeset::initFromICs_arclength(std::vector<double> IC, double t0, double tof, int numNodes, unsigned int ctrlLawID){
	SimEngine engine;
	engine.setMakeDefaultEvents(false);	// Don't use default crash events to avoid infinite loop
	engine.setRevTime(tof < 0);
	engine.setCtrlLaw(ctrlLawID);

	// Run the simulation and get the trajectory
	Traj traj(pSysData);
	engine.runSim(IC, t0, tof, &traj);

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
	int prevID = addNode(Node(traj.getStateByIx(0), traj.getEpochByIx(0)));

	// Loop through again to create trajectory
	sumArclen = 0;
	double sumTOF = 0;
	for(unsigned int s = 0; s < allArcLen.size(); s++){
		
		// Keep adding arclength between steps until the desired length is reached
		if(sumArclen < desiredArclen){
			sumArclen += allArcLen[s];
			sumTOF += allTOF[s];
		}else{
			// reached desired length: save node
			int id = addNode(Node(traj.getStateByIx(s), traj.getEpochByIx(s)));
			
			// Add a segment to link the previous node and this new one
			// I'm not setting the STM here; it will be initialized to 
			// Identity on every arc segment
			Segment newSeg = Segment(prevID, id, sumTOF);
			newSeg.setCtrlLaw(traj.getCtrlLawByIx(0));
			addSeg(newSeg);

			// Reset counters and index variables
			sumArclen = 0;
			sumTOF = 0;
			prevID = id;
		}
	}

	// Save the final state as the last node
	int id2 = addNode(Node(traj.getStateByIx(-1), traj.getEpochByIx(-1)));
	Segment newSeg = Segment(prevID, id2, sumTOF);
	newSeg.setCtrlLaw(traj.getCtrlLawByIx(0));
	addSeg(newSeg);
}//==========================================================

/**
 *	\brief Split a trajectory into nodes using the specified distribution type
 *
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory. Additionally,
 *	if the input trajectory contains discontinuities in position, velocity, or time, this 
 *	algorithm will fail to capture these discontinuities and could result in a very different
 *	path.
 *
 *	\param traj a trajectory to make into a nodeset
 *	\param numNodes the number of nodes to create, including IC
 *	\param type the node distribution type
 */
void Nodeset::initFromTraj(Traj traj, int numNodes, NodeDistro_tp type){
	/* Could I code this more intelligently? Probably. Am I too lazy? Definitely */ 
	initFromICs(traj.getStateByIx(0), traj.getEpochByIx(0), traj.getEpochByIx(-1) - traj.getEpochByIx(0), numNodes, type);
}//==============================================

}// END of Astrohelion namespace