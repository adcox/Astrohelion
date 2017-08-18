/**
 *  \file Arcset.cpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version April 28, 2017
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
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
#include <exception>

#include "Arcset.hpp"
#include "Exceptions.hpp"
#include "SysData.hpp"
#include "Utilities.hpp"

#include "SimEngine.hpp"
#include "Arcset.hpp"
#include "Arcset.hpp"

namespace astrohelion{


//-----------------------------------------------------
//      Constructors and Desctructor
//-----------------------------------------------------

/**
 *	\brief Create a arcset for a specific system
 *	\param data a pointer to a system data object
 */
Arcset::Arcset(const SysData *data) : BaseArcset(data) {}

/**
 *	\brief Create a arcset from another arcset
 *	\param t a arcset reference
 */
Arcset::Arcset(const Arcset &t) : BaseArcset(t) {}

/**
 *	\brief Create a arcset from its base class
 *	\param a an arc data reference
 */
Arcset::Arcset(const BaseArcset &a) : BaseArcset(a) {}

/**
 *  \brief Default destructor
 */
Arcset::~Arcset(){}

/**
 *  \brief Create a new arcset object on the stack
 *  \details The <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param sys pointer to a system data object
 *  \return a pointer to the newly created arcset
 */
baseArcsetPtr Arcset::create( const SysData *sys) const{
	return baseArcsetPtr(new Arcset(sys));
}//====================================================

/**
 *  \brief Create a new arcset object on the stack that is a 
 *  duplicate of this object
 *  \details The <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \return a pointer to the newly cloned arcset
 */
baseArcsetPtr Arcset::clone() const{
	return baseArcsetPtr(new Arcset(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  \brief Combine two arcsets.
 *  \details This function concatenates two arcset objects. It is assumed
 *  that the first state on <tt>rhs</tt> is identical to the final state on
 *  <tt>rhs</tt>. The <tt>rhs</tt> object is also assumed to occur after
 *  (chronologically) <tt>lhs</tt>
 * 
 *  \param lhs reference to a arcset object
 *  \param rhs reference to a arcset object
 * 
 *  \return the concatenation of lhs + rhs.
 */
Arcset operator +(const Arcset &lhs, const Arcset &rhs){
	const Arcset lhs_cpy(lhs);
	const Arcset rhs_cpy(rhs);
	Arcset result(lhs.pSysData);

	BaseArcset::sum(&lhs, &rhs, &result);

	return result;
}//====================================================

/**
 *  \brief Concatenate this object with another arcset
 * 
 *  \param rhs reference to a arcset object
 *  \return the concatenation of this and <tt>rhs</tt>
 *  @see operator +()
 */
Arcset& Arcset::operator +=(const Arcset &rhs){
	Arcset temp = *this + rhs;
	copyMe(temp);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get
//-----------------------------------------------------

/**
 *	\brief Allow velocity discontinuities (i.e., delta-Vs) at the specified segments
 *	\param id a vector of segment IDs that can have velocity discontinuities
 */
void Arcset::allowDV_at(std::vector<int> id) {
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
void Arcset::allowDV_all(){
	for(unsigned int i = 0; i < segs.size(); i++){
		segs[i].setVel_AllDiscon();
	}
}//====================================================

/**
 *  \brief Allow velocity discontinuities (i.e., delta-Vs) on none of the segments
 */
void Arcset::allowDV_none(){
	for(unsigned int i = 0; i < segs.size(); i++){
		segs[i].setVel_AllCon();
	}
}//====================================================

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
int Arcset::createNodesAtEvent(int segID, Event evt, double minTimeDiff){
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
int Arcset::createNodesAtEvents(int segID, std::vector<Event> evts, double minTimeDiff){
	if(segIDMap.count(segID) == 0)
		throw Exception("Arcset::createNodesAtEvents: Segment ID is out of bounds");

	// Get a copy of the segment we are replacing
	Segment seg = segs[segIDMap[segID]];
	Node origin = nodes[nodeIDMap[seg.getOrigin()]];
	Node terminus = nodes[nodeIDMap[seg.getTerminus()]];

	// Create a simulation engine and add the events to it
	SimEngine engine;
	engine.setVerbosity(Verbosity_tp::NO_MSG);
	engine.setRevTime(seg.getTOF() < 0);
	engine.setMakeDefaultEvents(false);
	for(unsigned int i = 0; i < evts.size(); i++){
		evts[i].setStopOnEvent(false);		// Ignore stopping conditions that other processes may have imposed
		engine.addEvent(evts[i]);
	}

	Arcset traj(pSysData);
	engine.runSim(origin.getState(), origin.getEpoch(), seg.getTOF(), &traj, seg.getCtrlLaw());

	int evtCount = 0;
	for(unsigned int e = 0; e < evts.size(); e++){
		for(unsigned int n = 0; n < traj.getNumNodes(); n++){
			if(traj.getNodeByIx(n).getTriggerEvent() == evts[e].getType()){
				evtCount++;
			}
		}
	}

	if(evtCount > 0){
		// Delete the old segment, replace it with new ones
		deleteSeg(segID);

		int prevNodeID = origin.getID();
		int nextNodeID = Linkable::INVALID_ID;
		for(unsigned int s = 0; s < traj.getNumSegs(); s++){
			// Get a copy of the segment from the newly propagated arc
			Segment newSeg = traj.getSegByIx(s);
			
			if(std::abs(newSeg.getTOF()) < std::abs(minTimeDiff)){
				// Too close to a previous node!
				throw Exception("Arcset::createNodesAtEvent: TOF between nodes is less than specified minimum; TODO - Implement a way to continue");
			}

			// Add the terminating node if this isn't the last segment
			if(s < traj.getNumSegs()-1){
				Node newNode = traj.getNode(newSeg.getTerminus());
				nextNodeID = addNode(newNode);
			}else{
				nextNodeID = terminus.getID();	// ID of the terminal node from the original segment we deleted
			}

			// Reset the link, update the origin, and add the segment to the arcset
			newSeg.clearLinks();
			newSeg.setOrigin(prevNodeID);
			newSeg.setTerminus(nextNodeID);
			addSeg(newSeg);
			
			// Shift these for the next roun
			prevNodeID = nextNodeID;
		}
	}
	
	return evtCount;
}//====================================================

/**
 *	\brief Retrieve the time along the arcset at a specific step
 *	\param ix node index; if < 0, it will count backwards from end of arcset
 *	\return the non-dimensional time along the arcset at the specified step
 *	\throws Exception if <tt>ix</tt> is out of bounds
 */
double Arcset::getTimeByIx(int ix) const {
	return getEpochByIx(ix);
}//====================================================


/**
 *  \brief Set the time associated with a node
 * 
 *	\param ix node index; if < 0, it will count backwards from end of arcset
 *  \param t time associated with the node
 *  \throws Exception if <tt>ix</tt> is out of bounds
 */
void Arcset::setTimeByIx(int ix, double t){
	setEpochByIx(ix, t);
}//====================================================

/**
 *  \brief Set each Segment STM to represent the relationship
 *  between the begining (chronologically) of the arcset and
 *  the the end (chronologicay) of each Segment.
 *  \details No checks are made to ensure the arcset is continuous.
 *  If the arcset is discontinuous the STMs will lose any 
 *  accuracy or meaning past the first (chronological) discontinuity.
 *  
 *  Only the single STM matrix is updated. The STM elements stored 
 *  in the Segment state vector are not changed; these elements represent
 *  the STM for each individual segment (begins at Identity for each 
 *  segment).
 *  
 *  \throws Exception if the multiplication of two consecutive STMs 
 *  is illegal (e.g., if the matrices have different sizes, which
 *  is the case if two segments leverage control laws with different 
 *  numbers of control states)
 */
void Arcset::setSTMs_parallel(){
	putInChronoOrder();	// returns immediately if already in chrono order

	MatrixXRd stmPrev, stmSeg;
	for(unsigned int s = 1; s < segs.size(); s++){
		if(segs[s].getCtrlLaw() != segs[s-1].getCtrlLaw())
			printWarn("Arcset::setSTMs_parallel: Segments leverage different control laws; multiplying STMs may yield non-useful / non-physical values\n");
		try{
			segs[s].setSTM(segs[s].getSTM() * segs[s-1].getSTM());
		}catch(std::exception &e){
			printErr(e.what());
			throw Exception("Arcset::setSTMs_parallel: Eigen error, cannot multiply two STMs; likely have a different size because of control laws");
		}
	}
}//====================================================

/**
 *  \brief Set each Segment STM to represent the evolution of the state
 *  vector between the beginning (chronologically) of the Segment and the
 *  end.
 *  \details The individual segment STMs elements are always stored in the
 *  Segment state vector
 */
void Arcset::setSTMs_sequence(){
	const unsigned int coreDim = pSysData->getDynamicsModel()->getCoreStateSize();

	unsigned int ctrlDim = 0;
	unsigned int stmDim = 0;
	ControlLaw *pLaw = nullptr;

	for(unsigned int s = 0; s < segs.size(); s++){
		pLaw = segs[s].getCtrlLaw();

		if(pLaw)
			ctrlDim = pLaw->getNumStates();
		else
			ctrlDim = 0;

		stmDim = (ctrlDim + coreDim)*(ctrlDim + coreDim);
		std::vector<double> statef = segs[s].getStateByRow(-1);
		std::vector<double> stm(statef.begin() + coreDim + ctrlDim, statef.begin() + coreDim + ctrlDim + stmDim);
		segs[s].setSTM(stm);
	}
}//====================================================

/**
 *  \brief Shift all time values by a constant amount
 *  \details This can be useful for use with the EM2SE and SE2EM functions
 * 
 *  \param amount a constant, non-dimensional time shift to apply to 
 *  all time values for points on this trajectory
 */
void Arcset::shiftAllTimes(double amount){
	for(unsigned int i = 0; i < nodes.size(); i++){
		nodes[i].setEpoch(nodes[i].getEpoch() + amount);
	}

	for(unsigned int s = 0; s < segs.size(); s++){
		segs[s].shiftAllTimes(amount);
	}
}//====================================================

//-----------------------------------------------------
//      Utility
//-----------------------------------------------------

/**
 *	\brief Print a useful message describing this arcset to the standard output
 */
void Arcset::print() const {

	printf("%s Arcset:\n Nodes: %zu\n Segments: %zu\n", pSysData->getTypeStr().c_str(),
		nodes.size(), segs.size());
	printf("List of Nodes:\n");
	for(const auto &index : nodeIDMap){
		printf("  %02d (ix %02d):", index.first, index.second);
		if(index.second != Linkable::INVALID_ID){
			std::vector<double> state = nodes[index.second].getState();
			printf(" @ %13.8f -- {", nodes[index.second].getEpoch());

			if(state.size() > 0){
				for(unsigned int i = 0; i < state.size()-1; i++){
					printf("%13.8f, ", state[i]);
				}
				printf("%13.8f}\n", state.back());
			}else{
				printf("}\n");
			}

			printf("\t> Link[0] = %02d,  Link[1] = %02d\n", nodes[index.second].getLink(0), 
				nodes[index.second].getLink(1));

			try{
				std::vector<double> ctrlState = nodes[index.second].getExtraParamVec(PARAMKEY_CTRL);
				printf("\t> Ctrl = {");
				
				if(ctrlState.size() > 0){
					for(unsigned int i = 0; i < ctrlState.size()-1; i++){
						printf("%13.8f, ", ctrlState[i]);
					}
					printf("%13.8f}\n", ctrlState.back());
				}else{
					printf("}\n");
				}
			}catch(Exception &e){}
			
		}else{
			printf(" [N/A]\n");
		}
	}

	printf("List of Segments:\n");
	for(const auto &index : segIDMap){
		printf("  %02d (ix %02d):", index.first, index.second);
		if(index.second != Linkable::INVALID_ID && index.second < static_cast<int>(segs.size())){
			printf(" origin @ %02d, terminus @ %02d, TOF = %13.8f\n", segs[index.second].getOrigin(),
				segs[index.second].getTerminus(), segs[index.second].getTOF());

			if(segs[index.second].getCtrlLaw())
				printf("\t> Ctrl Law: %s\n", segs[index.second].getCtrlLaw()->getLawString().c_str());
			else
				printf("\t> Ctrl Law: None\n");
			
		}else{
			printf(" [N/A]\n");
		}
	}

	printf("Constraints:\n-----------------\n");
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
 *	\brief Save the arcset to a file
 *	\param filename the name of the .mat file
 */
void Arcset::saveToMat(const char* filename, Save_tp saveTp) const{
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		astrohelion::printErr("Arcset::saveToMat: Error creating MAT file\n");
	}else{
		try{
			saveCmds(matfp, saveTp);
		}catch(Exception &E){
			Mat_Close(matfp);
			throw E;
		}catch(std::exception &e){
			Mat_Close(matfp);
			throw e;
		}
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
void Arcset::saveCmds(mat_t* pMatFile, Save_tp saveTp) const{
	saveLinkTable(pMatFile);

	saveNodeStates(pMatFile);
	saveNodeStateDeriv(pMatFile);
	saveNodeTimes(pMatFile);
	saveNodeCtrl(pMatFile);

	saveSegStates(pMatFile, saveTp);
	saveSegTimes(pMatFile, saveTp);
	saveSegTOF(pMatFile, saveTp);
	saveSegSTMs(pMatFile, saveTp);
	saveSegCtrlLaw(pMatFile, saveTp);

	pSysData->saveToMat(pMatFile);
}//====================================================

/**
 *  \brief Populate data in this trajectory from a matlab file
 * 
 *  \param filepath the path to the matlab data file
 *  \param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 *  
 *  \throws Exception if the data file cannot be opened
 */
void Arcset::readFromMat(const char *filepath, std::vector<ControlLaw*> &refLaws){

	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Arcset::loadFromMat Could not load data from file");
	}else{
		try{
			readCmds(matfp, refLaws);
		}catch(std::exception &e){
			Mat_Close(matfp);
			throw e;
		}
	}

	Mat_Close(matfp);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \details This function is called from readFromMat() and should
 *  be overridden in derived classes as necessary.
 * 
 *  \param pMatFile pointer to an open Matlab file
 *  \param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 *  
 *  \todo Remove backward compatibility code in future (today: May 3 2017)
 */
void Arcset::readCmds(mat_t *pMatFile, std::vector<ControlLaw*> &refLaws){
	Save_tp saveTp = Save_tp::SAVE_ALL;

	try{
		try{
			readLinkTable(pMatFile);
		}catch(Exception &e){
			printErr("Arcset::readCmds: Could not read link table\n\t%s\n", e.what());
			initNodesSegsFromMat(pMatFile);
		}

		readNodeStatesFromMat(pMatFile);
		readNodeStateDerivFromMat(pMatFile);
		readNodeTimesFromMat(pMatFile);
		readNodeCtrlFromMat(pMatFile);

		readSegCtrlLawFromMat(pMatFile, refLaws, saveTp);
		readSegSTMFromMat(pMatFile, saveTp);
		readSegStatesFromMat(pMatFile, saveTp);
		readSegTimesFromMat(pMatFile, saveTp);
		readSegTOFFromMat(pMatFile, saveTp);
		
	}catch(Exception &e){
		// if file was saved using older style, try slightly different read commands
		printErr("Arcset::readCmds: Encountered error:\n\t%s\n", e.what());
		
		try{
			// Old trajectory save
			initNodesSegsFromMat(pMatFile, VARNAME_DEP_STATE);
			readNodeStatesFromMat(pMatFile, VARNAME_DEP_STATE);
			readNodeTimesFromMat(pMatFile, VARNAME_DEP_TIME);
		}catch(Exception &e){
			// Old nodeset save
			printErr("Arcset::readCmds: Encountered error:\n\t%s\n", e.what());

			initNodesSegsFromMat(pMatFile, VARNAME_DEP_NODE);
			readNodeStatesFromMat(pMatFile, VARNAME_DEP_NODE);
			readNodeTimesFromMat(pMatFile, VARNAME_DEP_EPOCH);
		}

		// Old save for both trajectory and nodeset
		readNodeStateDerivFromMat(pMatFile);
		readSegCtrlLawFromMat(pMatFile, refLaws, Save_tp::SAVE_ALL);
		readSegSTMFromMat(pMatFile, Save_tp::SAVE_ALL);
		readSegTOFFromMat(pMatFile, Save_tp::SAVE_ALL);
	}
}//====================================================

}// End of astrohelion namespace