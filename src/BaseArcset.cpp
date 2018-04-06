/**
 *  @file BaseArcset.cpp
 *	@brief Data object that stores information about an integrated arc
 *	
 *	@author Andrew Cox
 *	@version Feb 24, 2017
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Astrohelion.
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

#include <algorithm>
#include <iostream>


#include "BaseArcset.hpp"

#include "AsciiOutput.hpp"
#include "DynamicsModel.hpp"
#include "EigenDefs.hpp"
#include "SysData.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Constructor (requires system data object)
 *	@param sys a pointer to a system data object that describes
 *	the system this trajectory is integrated in
 */
BaseArcset::BaseArcset(const SysData *sys) : pSysData(sys){}

/**
 *	@brief Copy constructor
 *	@param d an arcset reference
 */
BaseArcset::BaseArcset(const BaseArcset &d) : pSysData(d.pSysData){
	copyMe(d);
}//====================================================

/**
 *	@brief Destructor
 */
BaseArcset::~BaseArcset(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Set this object equal to another
 *	@param d an arcset reference
 *	@return a reference to this arcset object
 */
BaseArcset& BaseArcset::operator =(const BaseArcset &d){
	copyMe(d);

	pSysData = d.pSysData;
	return *this;
}//====================================================

/**
 *  @brief Sum two arcset objects.
 *  @details This function returns `result` = `lhs` + `rhs`; 
 *  Both `lhs` and `rhs` are copied and sorted into chronological 
 *  order. The `rhs` is then appended to the end of `lhs` with no
 *  time-of-flight between the end of `lhs` and the beginning of `rhs`.
 *  The first node of `rhs` is deleted if the final node on `lhs` is 
 *  an origin node, i.e., if the two progress in different time directions.
 * 
 *  @param lhs pointer to an arcset object
 *  @param rhs pointer to an arcset object
 *  @param result a pointer to the an arcset object that will store the result of
 *  the summing operation.
 */
void BaseArcset::sum(const BaseArcset *lhs, const BaseArcset *rhs, BaseArcset *result){
	baseArcsetPtr lhs_cpy = lhs->clone();
	baseArcsetPtr rhs_cpy = rhs->clone();

	if(!lhs_cpy->isInChronoOrder())
		lhs_cpy->putInChronoOrder();

	if(!rhs_cpy->isInChronoOrder())
		rhs_cpy->putInChronoOrder();

	int lhs_lastNodeID = lhs_cpy->getNodeByIx(-1).getID();
	int rhs_firstNodeID = rhs_cpy->getNodeByIx(0).getID();

	*result = *lhs_cpy;

	result->appendSetAtNode(rhs_cpy.get(), lhs_lastNodeID, rhs_firstNodeID, 0);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Add a constraint to the arcset object.
 *  @details The constraint application type determines whether the constraint
 *  is applied to a node, segment, or the arc as a whole.
 * 
 *  @param con the constraint
 *  @see ConstraintApp_tp
 *  @throws Exception if node or segment ID is out of bounds, or if the
 *  constraint application type is unknown
 */
void BaseArcset::addConstraint(Constraint con){
	int id = con.getID();
	switch(con.getAppType()){
		case ConstraintApp_tp::APP_TO_NODE:
			if(nodeIDMap.count(id) == 0)
				throw Exception("BaseArcset::addConstraint: Node ID out of bounds");

			nodes[nodeIDMap[id]].addConstraint(con);
			break;
		case ConstraintApp_tp::APP_TO_SEG:
			if(segIDMap.count(id) == 0)
				throw Exception("BaseArcset::addConstraint: Segment ID out of bounds");

			segs[segIDMap[id]].addConstraint(con);
			break;
		case ConstraintApp_tp::APP_TO_ARC:
			cons.push_back(con);
			break;
		default:
			throw Exception("BaseArcset::addConstraint: Constraint application type is unknown");
	}
}//=====================================================

/**
 *  @brief Add a node to this data object
 *  @details A unique key is assigned to the node when it is added.
 *  Any links the input node has are cleared.
 * 
 *  @param n the node to add
 *  @return the ID assigned to the node
 */
int BaseArcset::addNode(Node n){
	n.clearLinks();			// Cannot have any links coming in
	n.setID(nextNodeID);
	nodes.push_back(n);
	nodeIDMap[nextNodeID] = nodes.size()-1;

	bInChronoOrder = false;
	return nextNodeID++;
}//=====================================================

/**
 *  @brief Add a segment to the arcset object.
 *  @details Besides adding the segment object to the storage vector,
 *  this function updates the nodes associated with the origin and terminus
 *  of the segment so they contain the proper links. Additionally, the
 *  function checks to make sure time flows in one direction only and that
 *  parallel structures are not created (i.e., two branches extend from
 *  one node in same time direction).
 * 
 *  @param s The segment to add
 *  @return the ID of the segment
 *  @throws Exception if adding the segment will result in a time direction
 *  collision or a parallel structure. Exception is also thrown if the segment
 *  is linked to a non-existant node other than the placeholder 
 *  Linkable::INVALID_ID
 */
int BaseArcset::addSeg(Segment s){
	s.setID(nextSegID);
	
	if(s.getOrigin() == Linkable::INVALID_ID)
		throw Exception("BaseArcset::addSeg: Segment must have a valid origin node\n");
	
	bool foundValidLink = false;
	for(int i = 0; i < Linkable::NUM_LINKS; i++){
		// Get the ID of one of the nodes this segment is linked to
		int linkedNodeID = s.getLink(i);
		if(linkedNodeID != Linkable::INVALID_ID){

			if(nodeIDMap.count(linkedNodeID) == 0){
				char msg[128];
				sprintf(msg, "BaseArcset::addSeg: Segment (ID = %d) has link to an invalid node (ID = %d)\n",
					s.getID(), linkedNodeID);
				throw Exception(msg);
			}

			// Get the index of that node within the storage array
			int linkedNodeIx = nodeIDMap[linkedNodeID];
			if(linkedNodeIx != Linkable::INVALID_ID){

				// See if the node is linked to any other segments
				Node *theNode = &(nodes[linkedNodeIx]);
				int secondaryLinks = 0;
				for(int j = 0; j < Linkable::NUM_LINKS; j++){
					// Check to make sure the link is to a real segment
					if(theNode->getLink(j) != Linkable::INVALID_ID){
						int nearSegIx = segIDMap[theNode->getLink(j)];
						// Make sure the segment is real
						if(nearSegIx != Linkable::INVALID_ID){
							secondaryLinks++;

							// If the node is linked to another segment, get that segment and compare it to this one
							const Segment *nearSeg = &(segs[nearSegIx]);
							bool sameLinkType = nearSeg->getLink(i) == linkedNodeID;
							bool sameTimeDir = nearSeg->getTOF()*s.getTOF() > 0;

							if(sameLinkType && i == Segment::TERM_IX){
								// print();
								// printErr("Adding segment (ID %d) O: %d, T: %d\n", s.getID(), s.getOrigin(), s.getTerminus());
								// printErr("Conflict at node (ID %d): seg (ID %d) also terminates here!\n", linkedNodeID, nearSeg->getID());
								throw Exception("BaseArcset::addSeg: would create node with two terminating segments");
							}else if(sameLinkType && sameTimeDir){
								// either a time collision (both terminate) or parallel structure (both originate)
								// printf("Nearby segment w/ ID %d originates at node %d and terminates at node %d w/ TOF = %.4e\n", nearSeg->getID(), nearSeg->getOrigin(), nearSeg->getTerminus(), nearSeg->getTOF());
								// printf("The new seg (ID %d) originates at node %d and terminates at node %d w/ TOF = %.4e\n", s.getID(), s.getOrigin(), s.getTerminus(), s.getTOF());
								// print();
								// printInChrono();
								throw Exception("BaseArcset::addSeg: either time collision or parallel structure");
							}else if(!sameLinkType && !sameTimeDir){
								// parallel structure
								printErr("Parallel structure!\n");
								// print();
								printErr("Adding segment (ID %d) O: %d, T: %d, tof = %.4f\n", s.getID(), s.getOrigin(), s.getTerminus(), s.getTOF());
								printErr("Conflict at node (ID %d): seg (ID %d) has O: %d, T:%d, tof = %.4f\n", linkedNodeID, nearSeg->getID(),
									nearSeg->getOrigin(), nearSeg->getTerminus(), nearSeg->getTOF());
								throw Exception("BaseArcset::addSeg: parallel structure");
							}else{
								foundValidLink = true;
							}
						}
					}
				}

				// If neither of the nodes are linked to any segments, we should be fine!
				if(secondaryLinks == 0)
					foundValidLink = true;

				if(foundValidLink){
					theNode->addLink(nextSegID);		// OK, looks good from here; node will check for duplicate linkage
					// printf("Node (ID %d) is now linked to segment (ID %d)\n", theNode.getID(), nextSegID);
				}
			}else{
				// Linked node has valid ID but isn't part of this arcset
				throw Exception("BaseArcset::addSeg: Linked node is not part of this arcset object");
			}
		}
	}

	if(foundValidLink){
		segs.push_back(s);
		segIDMap[s.getID()] = segs.size()-1;
		bInChronoOrder = false;
		return nextSegID++;	
	}else{
		return Linkable::INVALID_ID;
	}
}//===========================================

/**
 *  @brief Append an arcset object (i.e., a set of nodes and segments) to this one
 * 
 *  @param pArcsetIn a pointer to the arcset object that will be appended to this object
 *  @param localNodeID the ID of the node in *this* arcset object that `pArcsetIn` will be linked to
 *  @param appendNodeID the ID of the node in `pArcsetIn` to link from
 *  @param tof time-of-flight between appendNodeID to localNodeID; when this value is nonzero, an artificial
 *  segment is constructed two link the two arcsets. If `tof` is zero, then one of the nodes is deleted.
 *  If the node at `localNodeID` is an origin node, then the node at `appendNodeID` is deleted. If
 *  the local node is not an origin, then it is deleted instead. The segment left without a terminus is then 
 *  connected to the remaining node.
 *  
 *  @return the ID of a new segment that links the old and new arcset objects
 *	@throws Exception if the two arcset objects have different system data objects
 *  @throws Exception if either ID is out of bounds
 *  @throws Exception if one or both of the identifies nodes does not have
 *  a free link slot
 */
int BaseArcset::appendSetAtNode(const BaseArcset *pArcsetIn, int localNodeID, int appendNodeID, double tof){
	if(pArcsetIn->pSysData != pSysData)
		throw Exception("BaseArcset::appendSetAtNode: Cannot concatenate two arcsets with different system data objects");

	// First, check to make sure the specified nodes are valid
	if(nodeIDMap.count(localNodeID) == 0)
		throw Exception("BaseArcset::appendSetAtNode: localNodeID is out of bounds");

	// Create a copy so we don't affect the original
	baseArcsetPtr pAppendArc = pArcsetIn->clone();

	const Node &localNode = nodes[nodeIDMap[localNodeID]];
	const Node &appendNode = pAppendArc->getNodeRef_const(appendNodeID);		// Will do its own index checks

	// Both nodes must have one "open port"
	if(!localNode.isLinkedTo(Linkable::INVALID_ID) || !appendNode.isLinkedTo(Linkable::INVALID_ID))
		throw Exception("BaseArcset::appendSetAtNode: specified nodes are not both open to a new link");

	// Determine if localNode is the origin or terminus of a segment
	// printf("linkToNode has links [%d, %d]\n", localNode.getLink(0), localNode.getLink(1));
	// printf("Choosing segment (ID %d)\n", localNode.getLink(0) == Linkable::INVALID_ID ? localNode.getLink(1) : localNode.getLink(0));
	Segment localSeg = getSeg(localNode.getLink(0) == Linkable::INVALID_ID ? localNode.getLink(1) : localNode.getLink(0));
	bool bLocalNodeIsOrigin = localSeg.getOrigin() == localNode.getID();
	Segment appendSeg = pAppendArc->getSeg(appendNode.getLink(0) == Linkable::INVALID_ID ? appendNode.getLink(1) : appendNode.getLink(0));
	bool bAppendNodeIsOrigin = appendSeg.getOrigin() == appendNode.getID();

	if(!bLocalNodeIsOrigin && !bAppendNodeIsOrigin)
		throw Exception("BaseArcset::appendSetAtNode: neither node is an origin; cannot create segment between them");

	// Store STM between two data sets; if tof != 0, we don't know the path between
	// the two, so linkSTM is initialized to the identity matrix. If tof = 0,
	// then the path is known (and will be deleted) so we store the STM in linkSTM
	int coreSize = pSysData->getDynamicsModel()->getCoreStateSize();
	MatrixXRd linkSTM = MatrixXRd::Identity(coreSize, coreSize);
	std::vector<double> linkSegStates, linkSegTimes;
	double linkTOF = tof;
	ControlLaw *pLinkCtrlLaw = nullptr;
	unsigned int linkStateWidth = 0;

	/* If TOF is zero, we need to connect the two arcsets. If localNode is an origin,
	 * we keep it, delete appendNode and copy appendSeg. However, if localNode is a terminus,
	 * then it makes more sense to keep appendNode (an origin), delete localNode, and copy
	 * localSeg
	 */
	if(tof == 0){
		if(bLocalNodeIsOrigin){
			/* 	Cases:
			 * 	
			 * 				vv Local Node
			 * 	[]<--Local--[] . []<--Append--[]
			 *	[]>--Local--[] . []--Append-->[]
			 *	[]<--Local--[] . []--Append--<[]
			 *					 ^^ Append Node
			 *					 
			 * 	In all these cases, the local node is an origin, so delete 
			 * 	appendNode and copy appendSeg
			 */
			linkTOF = appendSeg.getTOF();
			linkSTM = appendSeg.getSTM();
			linkSegStates = appendSeg.getStateVector();
			linkSegTimes = appendSeg.getTimeVector();
			pLinkCtrlLaw = appendSeg.getCtrlLaw();
			linkStateWidth = appendSeg.getStateWidth();

			// Get the ID of the other node attached to appendSeg before deleting appendSeg
			int otherNodeID = bAppendNodeIsOrigin ? appendSeg.getTerminus() : appendSeg.getOrigin();

			// Delete appendSeg and appendNode
			pAppendArc->deleteSeg(appendSeg.getID());
			pAppendArc->deleteNode(appendNodeID);

			// The node to be appended is now otherNode
			appendNodeID = otherNodeID;
			const Node &otherNode = pAppendArc->getNodeRef_const(otherNodeID);

			// The seg to be appended is now the segment attached to otherNode
			int otherSegID = otherNode.getLink(0) == Linkable::INVALID_ID ? otherNode.getLink(1) : otherNode.getLink(0);

			if(otherSegID == Linkable::INVALID_ID){
				// No segments left, just a node.
				// Leave appendSeg untouched; it is used later to determine the direction of time
				// Make bAppendNodeIsOrigin = true if the TOF is negative, false if TOF is positive
				// to avoid parallel structure problems
				 
				bAppendNodeIsOrigin = linkTOF < 0;
			}else{
				appendSeg = pAppendArc->getSeg(otherSegID);
				bAppendNodeIsOrigin = appendSeg.getOrigin() == otherNode.getID();
			}
		}else{	// appendNode is the origin, local node is terminus
			/*	Cases:
			 *	
			 *				vv Local Node
			 *	[]--Local-->[] . []--Append-->[]
			 *	[]--Local-->[] . []>--Append--[]	(NOT ALLOWED; both localNode and appendNode are terminus nodes)
			 *	[]--Local--<[] . []<--Append--[]	(NOT ALLOWED; both localNode and appendNode are terminus nodes)
			 *					 ^^ Append Node
			 *	
			 *	In this case, appendNode is an origin, so delete localNode
			 *	and copy over localSeg
			 */
			linkTOF = localSeg.getTOF();
			linkSTM = localSeg.getSTM();
			linkSegStates = localSeg.getStateVector();
			linkSegTimes = localSeg.getTimeVector();
			pLinkCtrlLaw = localSeg.getCtrlLaw();
			linkStateWidth = localSeg.getStateWidth();

			// Get the ID of the other node attached to localSeg before deleting localSeg
			int otherNodeID = localSeg.getOrigin();

			// Delete localSeg and localNode
			deleteSeg(localSeg.getID());
			deleteNode(localNodeID);

			// The node to be appended is now otherNode
			localNodeID = otherNodeID;
			const Node &otherNode = nodes[nodeIDMap[localNodeID]];

			// The seg to be appended is now the segment attached to otherNode
			int otherSegID = otherNode.getLink(0) == Linkable::INVALID_ID ? otherNode.getLink(1) : otherNode.getLink(0);

			if(otherSegID == Linkable::INVALID_ID){
				// No segments left, just a node.
				// Leave localSeg untouched; it is used later to determine the direction of time
				// Make bLocalNodeIsOrigin = true if the TOF is negative, false if TOF is positive
				 
				bLocalNodeIsOrigin = linkTOF < 0;
			}else{
				localSeg = segs[segIDMap[otherSegID]];
				bLocalNodeIsOrigin = localSeg.getOrigin() == otherNode.getID();
			}
		}
	}

	/*	Concatenate the two arcsets without connecting anything.
	 *	The mapping vector: index is the old node ID, value is the new node ID
	 *	All new IDs are initialized to the default pAppendArc value
	 */
	std::vector<int> map_oldID_to_newID = concatArcset(pAppendArc.get());

	// ------------------------------------------------------------------------
	// Add a new segment to link the nodes from [pAppendArc] to [this object]
	// ------------------------------------------------------------------------

	int linkOrigin = Linkable::INVALID_ID, linkTerminus = Linkable::INVALID_ID;

	if(!bLocalNodeIsOrigin){
		linkOrigin = localNodeID;
		linkTerminus = map_oldID_to_newID[appendNodeID];
	}else if(!bAppendNodeIsOrigin){
		linkOrigin = map_oldID_to_newID[appendNodeID];
		linkTerminus = localNodeID;
	}else{
		// Both are origins; the only double-linkOrigin node posibility is one where
		// each segment that originates from the node has a different time direction
		
		// localSeg and appendSeg have been updated to be the closest segments to the concatenation interface
		if(localSeg.getTOF() < 0){
			// If the closest segment on the local side is/was reverse time
			linkOrigin = tof > 0 ? localNodeID : map_oldID_to_newID[appendNodeID];
			linkTerminus = tof > 0 ? map_oldID_to_newID[appendNodeID] : localNodeID;
		}else if(appendSeg.getTOF() < 0){
			// If the closest segment on the local side is/was reverse time
			linkOrigin = tof > 0 ? map_oldID_to_newID[appendNodeID] : localNodeID;
			linkTerminus = tof > 0 ? localNodeID : map_oldID_to_newID[appendNodeID];
		}
	}

	if(tof != 0){
		// Put minimum amount of data in the linkSeg state and time vectors
		linkSegTimes.push_back(getEpoch(linkOrigin));
		linkSegTimes.push_back(getEpoch(linkTerminus));

		std::vector<double> q0 = getState(linkOrigin);
		std::vector<double> qf = getState(linkTerminus);
		std::vector<double> extra(pSysData->getDynamicsModel()->getExtraStateSize(), 0);

		// Append state, STM, and extra states for the linkOrigin and terimal nodes that the link segment connects to
		// Do not add any control information; the default control law is set to nullptr
		linkSegStates.insert(linkSegStates.end(), q0.begin(), q0.end());
		linkSegStates.insert(linkSegStates.end(), linkSTM.data(), linkSTM.data() + coreSize*coreSize);
		linkSegStates.insert(linkSegStates.end(), extra.begin(), extra.end());

		linkSegStates.insert(linkSegStates.end(), qf.begin(), qf.end());
		linkSegStates.insert(linkSegStates.end(), linkSTM.data(), linkSTM.data() + coreSize*coreSize);
		linkSegStates.insert(linkSegStates.end(), extra.begin(), extra.end());

		linkStateWidth = q0.size() + coreSize*coreSize + extra.size();
	}
	// print();
	bInChronoOrder = false;
	Segment linkSeg = Segment(linkOrigin, linkTerminus, linkTOF);
	linkSeg.setSTM(linkSTM);
	linkSeg.setCtrlLaw(pLinkCtrlLaw);
	linkSeg.setStateVector(linkSegStates);
	linkSeg.setStateWidth(linkStateWidth);
	linkSeg.setTimeVector(linkSegTimes);
	
	// linkSeg.print();
	int newSegID = addSeg(linkSeg);

	// Update all epochs
	updateEpochs(nodeIDMap[0], getEpoch(0));

	return newSegID;
}//====================================================

/**
 *  @brief Remove all constraints from this arcset object.
 *  @details Note that this does not affect any constraints placed on
 *  individual nodes or segments
 */
void BaseArcset::clearArcConstraints(){ cons.clear(); }


/**
 *  @brief Remove constraints from this arcset object as well as
 *  all its node and segment children
 */
void BaseArcset::clearAllConstraints(){
	for(unsigned int n = 0; n < nodes.size(); n++){ nodes[n].clearConstraints(); }
	for(unsigned int s = 0; s < segs.size(); s++){ segs[s].clearConstraints(); }
	cons.clear();
}//====================================================

/**
 *  @brief Concatenate two arcset objects
 *  @details The nodes, segments and constraints are copied from one arcset to another
 *  without creating or deleting any nodes or segments; i.e., the arcset object will
 *  include two independent "flows" without a segment to connect them
 * 
 *  @param pSet pointer to an arcset object
 *  @return a map relating the nodeIDs in `pSet` to the new IDs of the same nodes
 *  in this object; the index of the vector is the old node ID and the value is the 
 *  new node ID. If a node does not exist for one of the old ID values, a new value 
 *  equivalent to `Linkable::INVALID_ID` is stored in the associated 
 *  vector element.
 *  
 *  @throws Exception if the input arcset does not have the same system data object as this one
 */
std::vector<int> BaseArcset::concatArcset(const BaseArcset *pSet){
	if(pSet->pSysData != pSysData)
		throw Exception("BaseArcset::concatArcset: Cannot concatenate two arcsets with different system data objects");

	// A mapping vector: index is the old node ID, value is the new node ID
	// All new IDs are initialized to the default INVALID_ID value
	std::vector<int> map_oldID_to_newID(pSet->getNextNodeID(), Linkable::INVALID_ID);

	// Add all nodes from set to this object and keep track of new IDs
	for(unsigned int n = 0; n < pSet->getNumNodes(); n++){
		Node node = pSet->getNodeByIx(n);

		// Remove all links to segments; these will be added back when the segments are added to this new arcset object
		node.clearLinks();
		map_oldID_to_newID[node.getID()] = addNode(node);
	}

	// Add all segments from set to this object and update the link IDs
	// The act of adding the segment will update the links in the newly added nodes
	for(unsigned int s = 0; s < pSet->getNumSegs(); s++){
		Segment seg = pSet->getSegByIx(s);
		
		// Remap the origin and terminus to the new IDs
		if(seg.getOrigin() != Linkable::INVALID_ID)
			seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

		if(seg.getTerminus() != Linkable::INVALID_ID)
			seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

		addSeg(seg);
	}

	// Copy all constraints associated with the arc as a whole (node and segment constraints are contained in Node and Segment objects)
	std::vector<Constraint> arcCons = pSet->getArcConstraints();
	cons.insert(cons.end(), arcCons.begin(), arcCons.end());

	// Update tolerance to the larger of the two tolerances
	tol = pSet->getTol() > tol ? pSet->getTol() : tol;

	bInChronoOrder = false;
	return map_oldID_to_newID;
}//====================================================

/**
 *  @brief Delete the node with the specified ID
 *  @details In addition to deleting the desired node, the arcset object is "healed"
 *  so that time continuity is maintained. In the case of a linear-time set, the two segments
 *  on either side of the deleted node are combined. If the deleted node was the origin
 *  of two segments, they are still combined, but subject to a few extra conditions. First,
 *  if one of the segments connects to another segment and, therefore, doesn't link to another
 *  node, then the terminus of the other segment must be used. If both segments that originate
 *  from the deleted node terminate at other segments, an exception is thrown rather than attempting
 *  to step through mulitple segment links to identify the superposition of all the segments.
 * 
 *  @param id the ID of the node to delete; if the ID is out of range, an exception is
 *  thrown. If the ID is in range but doesn't represent an actual node, no action
 *  is taken.
 *  @throws Exception when:
 *  * The ID is out of bounds
 *  * The ID is in bounds but the associated node has been deleted (no longer exists)
 *  * Deleting the node will result in multiple segment interfaces; these must be 
 *  created more explicitely by the user and will not be created automatically by
 *  this function.
 */
void BaseArcset::deleteNode(int id){
	if( nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::deleteNode: id out of bounds");

	int nodeIx = nodeIDMap[id];

	if(nodeIx != Linkable::INVALID_ID){
		
		// printf("Attempting to delete node (ID %d)\n", id);
		// Get the node we're deleting
		const Node *pTheNode = &(nodes[nodeIx]);

		// Get the indices of any segments this node is linked to
		std::vector<int> linkedSegIxs;
		for(int i = 0; i < Linkable::NUM_LINKS; i++){
			if(pTheNode->getLink(i) != Linkable::INVALID_ID){
				int segIx = segIDMap[pTheNode->getLink(i)];
				if(segIx != Linkable::INVALID_ID){
					linkedSegIxs.push_back(segIx);
					// printf("  linked to segment (ID %d) at index %d\n", pTheNode->getLink(i), segIx);
				}
			}
		}

		// Determine what type of node this is
		if(linkedSegIxs.size() == 2){
			/* There are two possibilities: either this node is part of a linear time progression (one 
			 * terminating segment and one originating segment) or this node is a "double source" with 
			 * two nodes originating, each with different time directions
			 */
			if(segs[linkedSegIxs[0]].getTerminus() == id || segs[linkedSegIxs[1]].getTerminus() == id){
				// Get the segment that terminates at this node, and the segment that originates at this node
				int termSegIx = segs[linkedSegIxs[0]].getTerminus() == id ? 0 : 1;
				const Segment *pTermSeg = &(segs[linkedSegIxs[termSegIx]]);
				const Segment *pOrigSeg = &(segs[linkedSegIxs[(termSegIx + 1) % 2]]);

				// printf("  > Segment (ID %d) terminates at node (ID %d)\n", pTermSeg->getID(), id);
				// printf("  > Segment (ID %d) originates at node (ID %d)\n", pOrigSeg->getID(), id);

				// Just to check
				if(pTermSeg->getTOF()*pOrigSeg->getTOF() < 0){
					throw Exception("BaseArcset::deleteNode: I made an incorrect assumption about origin/terminus TOF direction!");
				}

				// Create a new segment
				Segment combo(pTermSeg->getOrigin(), pOrigSeg->getTerminus(), pTermSeg->getTOF() + pOrigSeg->getTOF());

				// Get the IDs of the segments; deleting id1 will change what pOrigSeg points to, which will affect the getID() return value
				int id1 = pTermSeg->getID(), id2 = pOrigSeg->getID();
				
				// Get control laws of the segments
				ControlLaw* law1 = pTermSeg->getCtrlLaw();
				ControlLaw* law2 = pOrigSeg->getCtrlLaw();

				// Use the control law information from pTermSeg, so use its state width info too
				combo.setStateWidth(pTermSeg->getStateWidth());

				// Replace the two segments with the new combined one
				deleteSeg(id1);	// CANNOT USE pTermSeg pointer AFTER THIS LINE
				deleteSeg(id2);	// CANNOT USE pOrigSeg pointer AFTER THIS LINE

				if((law1 && law2) && (*law1 != *law2))
					printWarn("BaseArcset::deleteNode: Node deleted between segments with different control laws. A new segment is created using only one of those control laws!\n");

				combo.setCtrlLaw(law1);
				addSeg(combo);
			}else{
				// Both must originate at node and have opposite time directions
				if(segs[linkedSegIxs[0]].getOrigin() != id || segs[linkedSegIxs[1]].getOrigin() != id){
					throw Exception("BaseArcset::deleteNode: double origin - made incorrect assumption about origin-ness");
				}

				if(segs[linkedSegIxs[0]].getTOF()*segs[linkedSegIxs[1]].getTOF() > 0){
					throw Exception("BaseArcset::deleteNode: double origin - made incorrect assumption about TOF");
				}

				int revSegIx = segs[linkedSegIxs[0]].getTOF() < 0 ? 0 : 1;
				const Segment *pRevSeg = &(segs[linkedSegIxs[revSegIx]]);
				const Segment *pForwardSeg = &(segs[linkedSegIxs[(revSegIx+1) % 2]]);

				/*	It is possible that, in this case, the segment that originates from this node and proceeds
				 * 	in reverse time does not terminate at a node, but links to a forward-propagated segment instead.
				 * 	If this is so, then the combination of the two segments becomes a reverse-time segment. However,
				 * 	if the segments both link to nodes, then the default behavior is to construct a new forward-time
				 * 	segment to replace the reverse and forward time segments that originated from this node.
				 */
				Segment combo;
				if(pRevSeg->getTerminus() != Linkable::INVALID_ID){
					combo = Segment(pRevSeg->getTerminus(), pForwardSeg->getTerminus(), std::abs(pRevSeg->getTOF()) + pForwardSeg->getTOF());
				}else{
					if(pForwardSeg->getTerminus() == Linkable::INVALID_ID)
						throw Exception("BaseArcset::deleteNode: Cannot delete node as both segments terminate at other segments");
					combo = Segment(pForwardSeg->getTerminus(), pRevSeg->getTerminus(), pRevSeg->getTOF() - pForwardSeg->getTOF());
				}

				// Get the IDs of the segments; deleting id1 will change what pForwardSeg points to, which will affect the getID() return value
				int id1 = pRevSeg->getID(), id2 = pForwardSeg->getID();

				// Use the control law information from pRevSeg, so use its state width info too
				combo.setStateWidth(pRevSeg->getStateWidth());

				// Get control laws of the segments
				ControlLaw* law1 = pRevSeg->getCtrlLaw();
				ControlLaw* law2 = pForwardSeg->getCtrlLaw();

				// Replace the two segments with the new one
				deleteSeg(id1);	// CANNOT USE pRevSeg pointer AFTER THIS LINE
				deleteSeg(id2);	// CANNOT USE pForwardSeg pointer AFTER THIS LINE

				if((law1 && law2) && (*law1 != *law2))
					printWarn("BaseArcset::deleteNode: Node deleted between "
						"segments with different control laws. A new segment "
						"is created using only one of those control laws!\n");

				combo.setCtrlLaw(law1);
				addSeg(combo);
			}
		}else if(linkedSegIxs.size() == 1){
			// Must be either the 1st or last point, so update the 1 linked seg
			segs[linkedSegIxs[0]].removeLink(id);
		}else{ /* Not linked to anything, just delete it! */ }


		nodes.erase(nodes.begin() + nodeIx);		// Remove node from vector
		nodeIDMap.erase(id);

		// Update vector that maps ID values to storage indices
		for(auto &index : nodeIDMap){
			if(index.second != Linkable::INVALID_ID && index.second > nodeIx){
				// Decrement all indices after the node we just removed
				index.second--;
			}
		}
		bInChronoOrder = false;
	}else{
		printf("Cannot Delete: Node with ID %02d was not located\n", id);
	}
}//====================================================

/**
 * @brief Delete a node at the specified index
 * @param ix index of the node in the nodes vector. If ix < 0, it will count
 * backward from the end of the vector
 * @throws Exception if `ix` is out of bounds
 */
void BaseArcset::deleteNodeByIx(int ix){
	while(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size())){
		char msg[128];
		sprintf(msg, "BaseArcset::deleteNodeByIx: ix = %d out of bounds", ix);
		throw Exception(msg);
	}else{
		deleteNode(nodes[ix].getID());
	}
}//====================================================

/**
 *  @brief Delete a segment with the specified ID.
 *  @details Additionally, any nodes linked to the specified segment are updated
 *  so that their link arrays no longer include a relationship with the soon-to-be
 *  deceased segment.
 * 
 *  @param id The ID of the segment to delete. If the ID is out of range, an exception
 *  is thrown. If the ID is in range but doesn't represent an existing segment,
 *  no deletion is made
 * 	@throws Exception if `id` is out of bounds
 */
void BaseArcset::deleteSeg(int id){
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::deleteSeg: Invalid ID (out of bounds)");

	int segIx = segIDMap[id];

	if(segIx != Linkable::INVALID_ID){
		// printf("Deleting segment (ID %d)\n", id);

		const Segment *pSeg = &(segs[segIx]);
		// printf("  Retrieved segment (ID %d)\n", pSeg->getID());
		for(int i = 0; i < Linkable::NUM_LINKS; i++){
			if(pSeg->getLink(i) != Linkable::INVALID_ID){
				int nodeIx = nodeIDMap[pSeg->getLink(i)];
				if(nodeIx != Linkable::INVALID_ID){
					// printf("  * Trying to remove a link to segment (ID %d) from node (ID %d)\n", pSeg->getID(), nodes[nodeIx].getID());
					nodes[nodeIx].removeLink(id);
				}else{
					// printf("Unable to remove link to segment (ID %d) from node (ID %d)\n", pSeg->getLink(i), id);
				}
			}
		}

		segs.erase(segs.begin() + segIx);	// CANNOT USE pSeg pointer AFTER THIS LINE
		segIDMap.erase(id);

		for(auto &index : segIDMap){
			if(index.second != Linkable::INVALID_ID && index.second > segIx){
				index.second--;
			}
		}
	}else{
		// printf("Cannot Delete: Segment with ID %02d was not located\n", id);
	}

	bInChronoOrder = false;
}//===========================================

/**
 * @brief Delete a segment at the specified index
 * @param ix index of the segment in the nodes vector. If ix < 0, it will count
 * backward from the end of the vector
 * @throws Exception if `ix` is out of bounds
 */
void BaseArcset::deleteSegByIx(int ix){
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size())){
		char msg[128];
		sprintf(msg, "BaseArcset::deleteSegByIx: ix = %d out of bounds", ix);
		throw Exception(msg);
	}else{
		deleteSeg(segs[ix].getID());
	}
}//====================================================

/**
 *  @brief Retrieve the state derivative vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the acceleration vector
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the node with the specified ID is not located in the nodeIDMap
 */
std::vector<double> BaseArcset::getStateDeriv(int id){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getStateDeriv: Node ID out of range");

	return getStateDerivByIx(nodeIDMap.at(id));
}//====================================================

/**
 *	@brief Retrieve an acceleration on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last acceleration)
 *	@return the acceleration associated with the specified index
 *	@throws Exception if `ix` is out of bounds
 */
std::vector<double> BaseArcset::getStateDerivByIx(int ix){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getStateDerivByIx: node index out of bounds");

	try{
		return nodes[ix].getExtraParamVec(PARAMKEY_STATE_DERIV);
	}catch(Exception &e){
		// Control laws are defined on segments, but state derivatives are stored on nodes...
			// Get the control law associated with the first segment linked to the node (should definitely be one)
			// CAUTION: If two segments link to a node with different control laws, this behavior may be undesirable.
			ControlLaw *pLaw = segs[segIDMap[nodes[ix].getLink(0)]].getCtrlLaw();
			EOM_ParamStruct params(pSysData, pLaw);

			std::vector<double> state = nodes[ix].getState();
			if(pLaw && pLaw->getNumStates() > 0){
				std::vector<double> ctrlState = nodes[ix].getExtraParamVec(PARAMKEY_CTRL);
				state.insert(state.end(), ctrlState.begin(), ctrlState.end());
			}

			std::vector<double> a = pSysData->getDynamicsModel()->getStateDeriv(nodes[ix].getEpoch(), state, &params);
			nodes[ix].setExtraParamVec(PARAMKEY_STATE_DERIV, a);
			return a;
	}
}//====================================================

/**
 *  @brief Retrieve a vector containing all the constraints applied to this arcset object.
 *  @details This vector does not include constraints placed on individual nodes or segments.
 *  @return a vector containing all the constraints applied to this arcset object.
 */
std::vector<Constraint> BaseArcset::getArcConstraints() const { return cons; }

/**
 *  @brief Retrieve a vector containing all the constraints applied to the arcset,
 *  nodes, and segments
 *  @details This vector DOES include constraints placed on individual nodes and segments
 *  @return a vector containing all the constraints applied to the arcset,
 *  nodes, and segments
 */
std::vector<Constraint> BaseArcset::getAllConstraints() const{
	std::vector<Constraint> allCons = cons;

	for(unsigned int n = 0; n < nodes.size(); n++){
		std::vector<Constraint> c = nodes[n].getConstraints();
		allCons.insert(allCons.end(), c.begin(), c.end());
	}
	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<Constraint> c = segs[s].getConstraints();
		allCons.insert(allCons.end(), c.begin(), c.end());
	}

	return allCons;
}//====================================================

/**
 *  @brief Determine what order to place the nodes and segments of this object
 *  into to achieve a chronological progression in forward time.
 * 
 *  @return a vector of ArcPiece objects that represent the chronological 
 *  order of the nodes and segments
 */
std::vector<ArcPiece> BaseArcset::getChronoOrder() const{
	// printf("* Beginning getChronoOrder()\n");
	std::vector<ArcPiece> pieces;
	if(nodes.size() == 0){
		astrohelion::printErr("BaseArcset::getChronoOrder: No nodes... exiting\n");
		return pieces;
	}

	return sortArcset(nodes[0].getID(), pieces);
}//====================================================

/**
 *  @brief Sort the arcset into chronological order beginning at a specified node
 *  @details This function operates recursively to sort all branches of the arcset
 * 
 *  @param ID the ID of a node in the arcset
 *  @param prevPieces already sorted ArcPiece objects from higher levels of the 
 *  recursive process
 *  @return a vector of ArcPiece objects that represent the chronological 
 *  order of the nodes and segments that make up the section of the arcset that
 *  contains the specified node
 *  
 *  @throws Exception if the ID is out of bounds
 */
std::vector<ArcPiece> BaseArcset::sortArcset(int ID, std::vector<ArcPiece> prevPieces) const{
	if(nodeIDMap.count(ID) == 0){
		printErr("Attempted to access ID = %d\n", ID);
		throw Exception("BaseArcset::sortArcset: ID out of bounds");
	}

	std::vector<ArcPiece> pieces;

	Node node0 = nodes[nodeIDMap.at(ID)];
	// printf("Beginning with node ID %d\n", node0.getID());

	pieces.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, node0.getID()));
	for(int dir = 1; dir > -2; dir -= 2){
		// printf("  Direction is %s\n", dir > 0 ? "[FORWARD]" : "[BACKWARD]");

		bool go = true;
		Node node = node0;
		
		while(go){
			// Get segments connected to node
			bool foundNextSeg = false;
			for(int i = 0; i < Linkable::NUM_LINKS; i++){
				if(node.getLink(i) != Linkable::INVALID_ID){

					Segment seg = getSeg(node.getLink(i));
					// printf("   (checking out segment ID %d): ", seg.getID());

					if(dir > 0 && ((seg.getTerminus() == node.getID() && seg.getTOF() < 0) || (seg.getOrigin() == node.getID() && seg.getTOF() > 0))){
						// This segment moves forward in time from node
						pieces.push_back(ArcPiece(ArcPiece::Piece_tp::SEG, seg.getID()));
						foundNextSeg = true;
						// printf("USE THIS\n");
					}else if(dir < 0 && ((seg.getTerminus() == node.getID() && seg.getTOF() > 0) || (seg.getOrigin() == node.getID() && seg.getTOF() < 0))){
						// This segment moves in reverse time from node
						pieces.insert(pieces.begin(), ArcPiece(ArcPiece::Piece_tp::SEG, seg.getID()));
						foundNextSeg = true;
						// printf("USE THIS\n");
					}else{
						foundNextSeg = false;
						// printf("do not use\n");
					}

					if(foundNextSeg){
						// Get next node
						int nextNodeID = Linkable::INVALID_ID;
						if(seg.getTerminus() == node.getID())
							nextNodeID = seg.getOrigin();
						else
							nextNodeID = seg.getTerminus();

						if(nextNodeID != Linkable::INVALID_ID){
							// printf("  Next node has ID %d\n", nextNodeID);
							node = getNode(nextNodeID);

							if(dir > 0)
								pieces.push_back(ArcPiece(ArcPiece::Piece_tp::NODE, node.getID()));
							else
								pieces.insert(pieces.begin(), ArcPiece(ArcPiece::Piece_tp::NODE, node.getID()));
						}else{
							// The segment terminates/originates without a node; Look for a link to
							// another segment OR end this section

							int linkedSegID = Linkable::INVALID_ID;
							for(unsigned int c = 0; c < cons.size(); c++){
								if(cons[c].getType() == Constraint_tp::SEG_CONT_PV){
									// Check to see if this constraint is linked to this node at all
									int ID0 = cons[c].getID(), ID1 = Linkable::INVALID_ID;
									std::vector<double> data = cons[c].getData();
									for(unsigned int d = 0; d < data.size(); d++){
										if(!std::isnan(data[d])){
											ID1 = static_cast<int>(data[d]);
											break;
										}
									}
									if(ID0 == seg.getID()){
										linkedSegID = ID1;
										break;
									}else if(ID1 == seg.getID()){
										linkedSegID = ID0;
										break;
									}
									// Didn't find the correct link, keep looking
								}
							}

							if(linkedSegID != Linkable::INVALID_ID){
								// printf(">Found seg-2-seg link with segment %d\n", linkedSegID);
								
								// Check to see if the identified segment has already been counted
								bool alreadyUsed = false;
								for(unsigned int p = 0; p < prevPieces.size(); p++){
									if(prevPieces[p].type == ArcPiece::Piece_tp::SEG && 
										prevPieces[p].id == linkedSegID){
										
										alreadyUsed = true;
										break;
									}
								}

								// Also check the current segment
								if(!alreadyUsed){
									for(unsigned int p = 0; p < pieces.size(); p++){
										if(pieces[p].type == ArcPiece::Piece_tp::SEG && pieces[p].id == linkedSegID){
											alreadyUsed = true;
											break;
										}
									}
								}

								if(!alreadyUsed){
									// Located a segment that links to this one, go recursive!
									Segment linkedSeg = segs[segIDMap.at(linkedSegID)];
									int linkedNodeID = linkedSeg.getOrigin();
									if(linkedNodeID != Linkable::INVALID_ID){
										// Concatenate the pieces I've arranged so far and the previously arranged ones; order is unimportant here
										std::vector<ArcPiece> allPieces = pieces;
										allPieces.insert(allPieces.begin(), prevPieces.begin(), prevPieces.end());
										std::vector<ArcPiece> section = sortArcset(linkedNodeID, allPieces);

										// Append newly arranged section of pieces to the total vector (order IS important here)
										pieces.insert(dir > 0 ? pieces.end() : pieces.begin(), section.begin(), section.end());
									}
								}
							}
							go = false;
						}

						break;	// Exit the loop
					}
				}
			}

			if(!foundNextSeg)
				go = false;		// Exit the iterations and switch directions
		}
	}

	return pieces;
}//====================================================

/**
 *	@brief Get a vector of one coordinate for all nodes
 *	@param ix the index of the state coordinate
 *	@return a vector containing the specified coordinate for all
 *	nodes (not necessarily in chronological order)
 *	@throws Exception if `ix` is out of bounds
 */
std::vector<double> BaseArcset::getCoord(unsigned int ix) const{
	if(nodes.size() > 0 && ix >= nodes[0].getState().size())
		throw Exception("BaseArcset::getCoord: Index out of range");

	std::vector<double> coord;
	for(unsigned int n = 0; n < nodes.size(); n++){
		coord.push_back(nodes[n].getState()[ix]);
	}

	return coord;
}//====================================================

/**
 *  @brief Retrieve the control law ID for a segment at the specified index
 *  @details [long description]
 * 
 *  @param ix Index of the segment within the storage vector. If ix < 0, it 
 *  will count backwards from the end of the storage vector.
 *  s
 *  @return control law ID for the specified segment
 *  @throws Exception if ix is out of bounds.
 */
const ControlLaw* BaseArcset::getCtrlLawByIx(int ix) const{
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("BaseArcset::getCtrlLawIDByIx: Index out of range");

	return segs[ix].getCtrlLaw();
}//====================================================

/**
 *  @brief Retrieve the epoch associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the epoch
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the node with the specified ID is not located in the nodeIDMap
 */
double BaseArcset::getEpoch(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getEpoch: Node ID out of range");

	int ix = nodeIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(nodes.size()) && ix >= 0){
		return nodes[ix].getEpoch();
	}else{
		throw Exception("BaseArcset::getEpoch: Could not locate the node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve the epoch of a specific node
 * 
 *  @param ix the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `n` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @return The epoch associated with the specified node
 *  @throws Exception if `ix` is out of bounds
 */
double BaseArcset::getEpochByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getEpochByIx: node index out of bounds");

	return nodes[ix].getEpoch();
}//=====================================================

/**
 *  @brief Retrieve a vector of all epoch values for the nodes
 *  @details Epochs are returned in the order corresponding
 *  to the nodes vector; to ensure chronological order, it is best
 *  to sort the arcset first.
 *  @return A vector with the epochs for all the nodes
 */
std::vector<double> BaseArcset::getEpochs() const{
	std::vector<double> time;
	for(unsigned int n = 0; n < nodes.size(); n++){
		time.push_back(nodes[n].getEpoch());
	}

	return time;
}//====================================================

/**
 *	@brief Retrieve a set of extra parameters for the specified node
 *	@param n the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `n` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@param key string that identifies the extra parameter
 *	@return a vector containing the extra parameter at the specified step and index
 *	@throws Exception if `n` is out of bounds
 */
double BaseArcset::getExtraParamByIx(int n, std::string key) const{
	if(n < 0)
		n += nodes.size();

	if(n < 0 || n >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getExtraParam: node index out of bounds");

	return nodes[n].getExtraParam(key);
}//====================================================

/**
 *  @brief Retrieve a set of extra parameters for the specified node
 * 
 *  @param n node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `n` is negative, this index will
 *	cound backwards from the end of the array.
 *  @param key string that identifies the extra parameter
 * 
 *  @return the set of extra parameters
 */
std::vector<double> BaseArcset::getExtraParamVecByIx(int n, std::string key) const{
	if(n < 0)
		n += nodes.size();

	if(n < 0 || n >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getExtraParam: node index out of bounds");

	return nodes[n].getExtraParamVec(key);
}//====================================================

/**
 *  @brief Retrieve the value of the ID that will be assigned to the next 
 *  node added to this arcset object
 *  @return the next node ID
 */
int BaseArcset::getNextNodeID() const { return nextNodeID; }

/**
 *  @brief Retrieve the value of the ID that will be assigned to the next 
 *  segment added to this arcset object
 *  @return the next segment ID
 */
int BaseArcset::getNextSegID() const { return nextSegID; }

/**
 *	@brief Retrieve the number of nodes
 *	@return the number of nodes
 */
unsigned int BaseArcset::getNumNodes() const { return nodes.size(); }

/**
 *	@brief Retrieve the number of segments
 *	@return the number of segments
 */
unsigned int BaseArcset::getNumSegs() const { return segs.size(); }

/**
 *  @brief Retrieve the total number of constraints contained by all nodes, segments,
 *  and the arcset object itself
 *  @return the total number of constraints applied to this object and its children
 */
unsigned int BaseArcset::getNumCons() const {
	unsigned int conCount = cons.size();
	for(unsigned int n = 0; n < nodes.size(); n++){ conCount += nodes[n].getNumCons(); }
	for(unsigned int s = 0; s < segs.size(); s++){ conCount += segs[s].getNumCons(); }

	return conCount;
}//===================================================

/**
 *  @brief Retrieve a specific node
 * 
 *  @param id the ID of the desired node.
 *	
 *  @return the node located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no node exists with the specified ID
 */
Node BaseArcset::getNode(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getNode: Invalid ID (out of bounds)");

	int ix = nodeIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(nodes.size()) && ix >= 0){
		return nodes[ix];
	}else{
		throw Exception("BaseArcset::getNode: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a node based on its index in the storage array
 * 
 *  @param ix The index of the node; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a node at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
Node BaseArcset::getNodeByIx(int ix) const{
	while(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size())){
		printErr("Attempted to access index %d\n", ix);
		throw Exception("BaseArcset::getNodeByIx: Index out of bounds");
	}

	return nodes[ix];
}//=====================================================

/**
 *  @brief Retrieve a reference to a specific node
 * 
 *  @param id the ID of the desired node.
 *	
 *  @return the node located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no node exists with the specified ID
 */
Node& BaseArcset::getNodeRef(int id){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getNodeRef: Invalid ID (out of bounds)");

	int ix = nodeIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(nodes.size()) && ix >= 0){
		bInChronoOrder = false;
		return nodes[ix];
	}else{
		throw Exception("BaseArcset::getNodeRef: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a reference to a node based on its index in the storage array
 * 
 *  @param ix The index of the node; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a reference to a node at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
const Node& BaseArcset::getNodeRefByIx_const(int ix) const{
	while(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getNodeRefByIx: Index out of bounds");

	return nodes[ix];
}//=====================================================

/**
 *  @brief Retrieve a reference to a specific node
 * 
 *  @param id the ID of the desired node.
 *	
 *  @return the node located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no node exists with the specified ID
 */
const Node& BaseArcset::getNodeRef_const(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getNodeRef: Invalid ID (out of bounds)");

	int ix = nodeIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(nodes.size()) && ix >= 0){
		return nodes[ix];
	}else{
		throw Exception("BaseArcset::getNodeRef: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a reference to a node based on its index in the storage array
 * 
 *  @param ix The index of the node; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a reference to a node at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
Node& BaseArcset::getNodeRefByIx(int ix){
	while(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getNodeRefByIx: Index out of bounds");

	bInChronoOrder = false;
	return nodes[ix];
}//=====================================================

/**
 *  @brief Retrieve the index of a specific node within the node storage vector
 * 
 *  @param id node ID
 *  @return the index of the node with the specified ID within the storage vector
 *  @throws Exception if `id` is out of bounds
 */
int BaseArcset::getNodeIx(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getNodeIx: Inavlid ID; out of bounds");

	return nodeIDMap.at(id);
}//====================================================

/**
 *  @brief Retrieve a specific node
 * 
 *  @param id the ID of the desired node
 *	
 *  @return the node located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no segment exists with the specified ID
 */
Segment BaseArcset::getSeg(int id) const{
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getSeg: Invalid ID (out of bounds)");

	int ix = segIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(segs.size()) && ix >= 0){
		return segs[ix];
	}else{
		throw Exception("BaseArcset::getSeg: Could not locate a segment with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a segment based on its index in the storage array
 * 
 *  @param ix The index of the segment; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a segment at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
Segment BaseArcset::getSegByIx(int ix) const{
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("BaseArcset::getSegByIx: Index out of bounds");

	return segs[ix];
}//=====================================================

/**
 *  @brief Retrieve a reference to a specific segment
 * 
 *  @param id the ID of the desired segment.
 *	
 *  @return the segment located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no segment exists with the specified ID
 */
Segment& BaseArcset::getSegRef(int id){
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getSegRef: Invalid ID (out of bounds)");

	int ix = segIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(segs.size()) && ix >= 0){
		bInChronoOrder = false;
		return segs[ix];
	}else{
		throw Exception("BaseArcset::getSegRef: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a reference to a segment based on its index in the storage array
 * 
 *  @param ix The index of the segment; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a reference to a segment at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
Segment& BaseArcset::getSegRefByIx(int ix){
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("BaseArcset::getSegRefByIx: Index out of bounds");

	bInChronoOrder = false;
	return segs[ix];
}//=====================================================

/**
 *  @brief Retrieve a reference to a specific segment
 * 
 *  @param id the ID of the desired segment.
 *	
 *  @return the segment located with the specified ID
 *  @throws Exception if `id` is out of bounds or if no segment exists with the specified ID
 */
const Segment& BaseArcset::getSegRef_const(int id) const{
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getSegRef: Invalid ID (out of bounds)");

	int ix = segIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(segs.size()) && ix >= 0){
		return segs[ix];
	}else{
		throw Exception("BaseArcset::getSegRef: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a reference to a segment based on its index in the storage array
 * 
 *  @param ix The index of the segment; if `ix` is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a reference to a segment at the specified index
 *  @throws Exception if `ix` is out of bounds
 */
const Segment& BaseArcset::getSegRefByIx_const(int ix) const{
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("BaseArcset::getSegRefByIx: Index out of bounds");

	return segs[ix];
}//=====================================================

/**
 *  @brief Retrieve the index of a specific node within the node storage vector
 * 
 *  @param id node ID
 *  @return the index of the node with the specified ID within the storage vector
 *  @throws Exception if `id` is out of bounds
 */
int BaseArcset::getSegIx(int id) const{
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getSegIx: Inavlid ID; out of bounds");

	return segIDMap.at(id);
}//====================================================

/**
 *  @brief Retrieve the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the state vector
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the node with the specified ID is not located in the nodeIDMap
 */
std::vector<double> BaseArcset::getState(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::getState: Node ID out of range");

	int ix = nodeIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(nodes.size()) && ix >= 0){
		return nodes[ix].getState();
	}else{
		throw Exception("BaseArcset::getState: Could not locate the node with the specified ID");
	}
}//====================================================

/**
 *	@brief Retrieve a position-velocity state on the arc
 *	@param ix the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `n` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@return the state associated with the specified index
 *	@throws Exception if `ix` is out of bounds
 */
std::vector<double> BaseArcset::getStateByIx(int ix) const{
	while(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::getStateByIx: node index out of bounds");

	return nodes[ix].getState();
}//====================================================

/**
 *  @brief Retrieve the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the STM
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the segment with the specified ID is not located in the segIDMap
 */
MatrixXRd BaseArcset::getSTM(int id) const{
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getSTM: Node ID out of range");

	int ix = segIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(segs.size()) && ix >= 0){
		return segs[ix].getSTM();
	}else{
		throw Exception("BaseArcset::getSTM: Could not locate the segment with the specified ID");
	}
}//====================================================

/**
 *	@brief Retrieve an STM on the arc
 *	@param ix the segment index. If it is negative, the index will count backwards
 *	from the end of the `segs` storage array
 *	
 *	@return the STM associated with the specified index
 *	@throws Exception if `ix` is out of bounds
 */
MatrixXRd BaseArcset::getSTMByIx(int ix) const{
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size())){
		printErr("Attempted to reach STM with ix = %d (max = %zu)\n", ix, segs.size());
		print();
		throw Exception("BaseArcset::getSTMByIx: segment index out of bounds");
	}

	return segs[ix].getSTM();
}//====================================================

/**
 *	@brief Retrieve the a pointer to the system data object associated with this arc
 *	@return a pointer to the system data object associated with this arc
 */
const SysData* BaseArcset::getSysData() const { return pSysData; }

/**
 *  @brief Retrieve the time-of-flight associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the time-of-flight
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the segment with the specified ID is not located in the segIDMap
 */
double BaseArcset::getTOF(int id) const{
	if(segIDMap.count(id) == 0)
		throw Exception("BaseArcset::getTOF: Segment ID out of range");

	int ix = segIDMap.at(id);
	if(ix != Linkable::INVALID_ID && ix < static_cast<int>(segs.size()) && ix >= 0){
		return segs[ix].getTOF();
	}else{
		throw Exception("BaseArcset::getTOF: Could not locate the segment with the specified ID");
	}
}//====================================================
/**
 *	@brief Get the time-of-flight for a specific segment
 *	@param ix node index (NOT the ID); if less than 0, the index counts
 *	backwards from the end of the nodeset
 *	@return non-dimensional time-of-flight along the specified segment
 *	@throws Exception if `ix` is out of bounds
 */
double BaseArcset::getTOFByIx(int ix) const {
	while(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("Arcset::getTOFByIx: invalid segment index");

	return segs[ix].getTOF();
}//====================================================

/**
 *	@brief Retrieve the tolerance with which data in this object was computed
 *	@return the tolerance with which data in this object was computed
 */
double BaseArcset::getTol() const { return tol; }

/**
 *  @brief Determine the total time-of-flight along this arc.
 *  @details This function sums the TOF along each segment; derived
 *  classes may override this function to use different methods.
 *  @return the total time-of-flight along this arc, units consistent
 *  with the SysData object
 */
double BaseArcset::getTotalTOF() const{
	double total = 0;
	for(unsigned int s = 0; s < segs.size(); s++){
		// total += std::abs(segs[s].getTOF());
		total += segs[s].getTOF();
	}
	return total;
}//=================================================

/**
 *  @brief Determine if the arcset is arranged in chronological order
 *  @details This is sufficient to prove that the arcset has
 *  been sorted, but not necessary; i.e., even if the flag is
 *  false, the arcset may indeed be in chronological order, but
 *  the oposite can not be true.
 *  @return a flag indicating if the arcset has been sorted into chronological order
 */
bool BaseArcset::isInChronoOrder() const{ return bInChronoOrder; }

/**
 *  @brief Rearrange the nodes and segments so that they are listed
 *  in chronological order in their storage arrays.
 *  
 *  @details This does not change the ID of any of the nodes or segments,
 *  only their index within the storage array. After calling this function,
 *  accessing the -1 node or segment will return the latest (in time) object.
 *  
 *  @param force if true, the arcset will be sorted regardless of the value
 *  of the isInChronoOrder() flag.
 *  
 *  @throws Exception if the getChronoOrder() sorting algorithm returns a 
 *  set of ArcPiece objects that has a different size than the combined 
 *  node and segment vectors, the function is aborted as it is likely a node or 
 *  segment was skipped and we don't want to lose information.
 */
void BaseArcset::putInChronoOrder(bool force){
	if(bInChronoOrder && !force)
		return;	// already sorted

	// First, determine the chronological order
	// print();
	std::vector<ArcPiece> pieces = getChronoOrder();

	if(pieces.size() != nodes.size() + segs.size()){
		astrohelion::printErr("Pieces has %zu elements, but there are %zu nodes and %zu segs\n", pieces.size(), nodes.size(), segs.size());
		saveToMat("ChronoOrderErr.mat");
		throw Exception("BaseArcset::putInChronoOrder: The sorted vector does not include all nodes and segments; aborting to avoid losing data\n");
	}

	std::vector<Node> newNodes;
	std::vector<Segment> newSegs;
	std::map<int, int> newNodeIDMap;
	std::map<int, int> newSegIDMap;

	for(unsigned int i = 0; i < pieces.size(); i++){
		if(pieces[i].type == ArcPiece::Piece_tp::NODE){
			Node node = getNode(pieces[i].id);

			// Update the links to include the segment before and after this node
			node.clearLinks();
			if(i > 0 && pieces[i-1].type == ArcPiece::Piece_tp::SEG)
				node.addLink(pieces[i-1].id);

			if(i+1 < pieces.size() && pieces[i+1].type == ArcPiece::Piece_tp::SEG)
				node.addLink(pieces[i+1].id);

			newNodeIDMap[node.getID()] = newNodes.size();
			newNodes.push_back(node);
		}else if(pieces[i].type == ArcPiece::Piece_tp::SEG){
			Segment seg = getSeg(pieces[i].id);
			seg.setTOF(std::abs(seg.getTOF()));	// Time should always progress forwards!
			
			// Update the links to include the node before and after this segment
			seg.clearLinks();
			if(i > 0 && pieces[i-1].type == ArcPiece::Piece_tp::NODE)
				seg.setOrigin(pieces[i-1].id);

			if(i+1 < pieces.size() && pieces[i+1].type == ArcPiece::Piece_tp::NODE)
				seg.setTerminus(pieces[i+1].id);

			newSegIDMap[seg.getID()] = newSegs.size();
			newSegs.push_back(seg);
		}
	}

	nodes = newNodes;
	segs = newSegs;
	nodeIDMap = newNodeIDMap;
	segIDMap = newSegIDMap;

	bInChronoOrder = true;
}//=============================================

/**
 *  @brief Set the epoch time associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param epoch the epoch time
 *  @throws Exception if `id` is out of bounds
 */
void BaseArcset::setEpoch(int id, double epoch){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::setEpoch: Node ID out of range");

	nodes[nodeIDMap[id]].setEpoch(epoch);
}//====================================================

/**
 *  @brief Set the epoch time for a specific node
 * 
 *  @param ix the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `ix` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param epoch the epoch time
 *  @throws Exception if `ix` is out of bounds
 */
void BaseArcset::setEpochByIx(int ix, double epoch){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::setEpochByIx: Node index out of bounds");

	nodes[ix].setEpoch(epoch);
}//====================================================

/**
 *  @brief Set the state derivative vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param qdot the acceleration vector
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if `qdot` does not have the same size as the 
 *  state vector
 */
void BaseArcset::setStateDeriv(int id, std::vector<double> qdot){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::setStateDeriv: Node ID out of range");

	if(qdot.size() != pSysData->getDynamicsModel()->getCoreStateSize())
		throw Exception("BaseArcset::setStateDeriv: state derivative vector must have the same dimension as the state vector");

	nodes[nodeIDMap[id]].setExtraParamVec(PARAMKEY_STATE_DERIV, qdot);
}//====================================================

/**
 *  @brief Set the state derivative vector for a specific step/node
 * 
 *  @param ix the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `n` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param derivVec state derivative vector
 *  @throws Exception if `ix` is out of bounds
 *  @throws Exception if `qdot` does not have the same size as the 
 *  state vector
 */
void BaseArcset::setStateDerivByIx(int ix, std::vector<double> derivVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::setStateDerivByIx: node index out of bounds");

	if(derivVec.size() != pSysData->getDynamicsModel()->getCoreStateSize())
		throw Exception("BaseArcset::setStateDeriv: state derivative vector must have the same dimension as the state vector");

	nodes[ix].setExtraParamVec(PARAMKEY_STATE_DERIV, derivVec);
}//=================================================

/**
 *  @brief Set the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param state the state vector
 *  @throws Exception if `id` is out of bounds
 */
void BaseArcset::setState(int id, std::vector<double> state){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::setState: Node ID out of range");

	nodes[nodeIDMap[id]].setState(state);
}//====================================================

/**
 *  @brief Set the state vector for a specific node
 * 
 *  @param ix the node index within the `nodes` storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If `ix` is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param stateVec vector of non-dimensional state values
 *  @throws Exception if `ix` is out of bounds
 */
void BaseArcset::setStateByIx(int ix, std::vector<double> stateVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= static_cast<int>(nodes.size()))
		throw Exception("BaseArcset::setStateByIx: node index out of bounds");

	nodes[ix].setState(stateVec);
}//=================================================

/**
 *  @brief Set the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @param stm the STM
 *  @throws Exception if `id` is out of bounds
 *  @throws Exception if the STM is not the size specified by the DynamicalModel
 */
void BaseArcset::setSTM(int id, MatrixXRd stm){
	if(nodeIDMap.count(id) == 0)
		throw Exception("BaseArcset::setSTM: Node ID out of range");

	int stateSize = pSysData->getDynamicsModel()->getCoreStateSize();
	if(stm.rows() != stateSize || stm.cols() != stateSize)
		throw Exception("BaseArcset::setSTMByIx: STM size does not match the state size for this dynamical system");

	segs[segIDMap[id]].setSTM(stm);
}//====================================================

/**
 *  @brief Set the STM for a specific step/node
 * 
 *  @param ix index of the segment with the `segs` storage array; if it is negative,
 *  it will count backwards from the end of the array.
 *  
 *  @param stm a matrix containing the STM
 *  @throws Exception if `ix` is out of bounds
 *  @throws Exception if the STM is not the size specified by the DynamicalModel
 */
void BaseArcset::setSTMByIx(int ix, MatrixXRd stm){
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= static_cast<int>(segs.size()))
		throw Exception("BaseArcset::setSTMByIx: node index out of bounds");

	int stateSize = pSysData->getDynamicsModel()->getCoreStateSize();
	if(stm.rows() != stateSize || stm.cols() != stateSize)
		throw Exception("BaseArcset::setSTMByIx: STM size does not match the state size for this dynamical system");

	segs[ix].setSTM(stm);
}//=================================================

/**
 *	@brief Set the computational tolerance for this data object
 *	@param d the tolerance
 */
void BaseArcset::setTol(double d){ tol = d; }

/**
 *  @brief Update the epochs of all nodes such that time is continuous.
 *  @details By specifying the epoch of one node in the set, all other
 *  nodes are updated using the segment times-of-flight between them.
 * 
 *  @param nodeID the ID of a node
 *  @param epoch the epoch of the node with the specified ID.
 */
void BaseArcset::updateEpochs(int nodeID, double epoch){
	if(nodeIDMap.count(nodeID) == 0)
		throw Exception("BaseArcset::updateEpochs: Invalide node ID");

	std::vector<ArcPiece> pieces = getChronoOrder();

	// Update epoch of the specified node
	nodes[nodeIDMap[nodeID]].setEpoch(epoch);

	std::vector<ArcPiece>::iterator pieceIt = std::find(pieces.begin(), pieces.end(), ArcPiece(ArcPiece::Piece_tp::NODE, nodeID));

	if(pieceIt == pieces.end()){
		astrohelion::printErr("BaseArcset::updateEpochs: Could not find the node with the specified ID... INVESTIGATE THIS, ANDREW!\n");
		return;
	}

	int ixInPieces = pieceIt - pieces.begin();

	// Move outward from that node in the chronological map and
	// update all other nodes
	double ellapsed;
	for(int stepDir = 1; stepDir > -2; stepDir -= 2){
		ellapsed = 0;
		for(int i = ixInPieces + stepDir; i < static_cast<int>(pieces.size()) && i >= 0; i += stepDir){
			if(pieces[i].type == ArcPiece::Piece_tp::SEG){
				ellapsed += segs[segIDMap[pieces[i].id]].getTOF();
			}else if(pieces[i].type == ArcPiece::Piece_tp::NODE){
				nodes[nodeIDMap[pieces[i].id]].setEpoch(epoch+ellapsed);
			}
		}
	}
}//====================================================

//-------------------------------------------------------------------------------------
//      General Utility Functions
//-------------------------------------------------------------------------------------

/**
 *	@brief Copy all data from the input arc data to this one
 *	@param d an arc data object reference
 */
void BaseArcset::copyMe(const BaseArcset &d){
	nodes = d.nodes;
	nodeIDMap = d.nodeIDMap;
	segs = d.segs;
	segIDMap = d.segIDMap;
	cons = d.cons;
	tol = d.tol;
	nextNodeID = d.nextNodeID;
	nextSegID = d.nextSegID;
	bInChronoOrder = d.bInChronoOrder;
}//====================================================

/**
 *  @brief Delete all data from this object and reset all parameters
 */
void BaseArcset::reset(){
	nodes.clear();
	nodeIDMap.clear();
	segs.clear();
	segIDMap.clear();
	cons.clear();
	tol = 0;
	nextNodeID = 0;
	nextSegID = 0;
	bInChronoOrder = false;
}//====================================================

/**
 *  @brief Print a ASCII graphic of the arcset in chronological order
 */
void BaseArcset::printInChrono() const{
	std::vector<ArcPiece> pieces = getChronoOrder();

	for(unsigned int i = 0; i < pieces.size(); i++){
		ArcPiece p = pieces[i];

		if(p.type == ArcPiece::Piece_tp::NODE){
			printf("[%02d]", p.id);
		}else if(p.type == ArcPiece::Piece_tp::SEG){
			if(getSeg(p.id).getTOF() > 0)
				printf("--(%02d)->", p.id);
			else
				printf(">-(%02d)--", p.id);
		}
	}
	printf("\n");
}//====================================================

/**
 *  @brief Print nodeIDMap to standard output
 */
void BaseArcset::printNodeIDMap() const{
	int count = 0;
	for(const auto &index : nodeIDMap){
		if(count % 20 == 0)
			printf("----------------\n%4s -> %4s\n----------------\n", "ID", "Ix");

		printf("%4d -> %4d\n", index.first, index.second);

		count++;
	}
}//====================================================

/**
 *  @brief Print segIDMap to standard output
 */
void BaseArcset::printSegIDMap() const{
	int count = 0;
	for(const auto &index : segIDMap){
		if(count % 20 == 0)
			printf("----------------\n%4s -> %4s\n----------------\n", "ID", "Ix");

		printf("%4d -> %4d\n", index.first, index.second);

		count++;
	}
}//====================================================

//-------------------------------------------------------------------------------------
//      File I/O Utility Functions
//-------------------------------------------------------------------------------------

/**
 *  @brief Create a matio variable for the link table
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_LinkTable(const char *pVarName) const{
	unsigned int numSegs = segs.size();
	std::vector<int> segTable(numSegs*4, Linkable::INVALID_ID);

	// Store data in column-major order
	for(unsigned int s = 0; s < numSegs; s++){
		segTable[0*numSegs + s] = segs[s].getID();

		if(segs[s].getOrigin() != Linkable::INVALID_ID)
			segTable[1*numSegs + s] = nodeIDMap.at(segs[s].getOrigin());

		if(segs[s].getTerminus() != Linkable::INVALID_ID)
			segTable[2*numSegs + s] = nodeIDMap.at(segs[s].getTerminus());

		segTable[3*numSegs + s] = astrohelion::sign(segs[s].getTOF());
	}

	// DO copy the data (it is freed at the end of the function!)
	size_t dims[2] = {numSegs, 4};
	return Mat_VarCreate(pVarName, MAT_C_INT32, MAT_T_INT32, 2, dims, 
		&(segTable.front()), 0);
}//====================================================

/**
 *  @brief Create a matio variable for the constraints
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_Constraints(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	// Step 1: Gather all the constraints
	std::vector<Constraint> allCons = cons;
	matvar_t* pMatVar = nullptr;

	for(unsigned int n = 0; n < nodes.size(); n++){
		std::vector<Constraint> temp = nodes[n].getConstraints();
		allCons.insert(allCons.end(), temp.begin(), temp.end());
	}

	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<Constraint> temp = segs[s].getConstraints();
		allCons.insert(allCons.end(), temp.begin(), temp.end());
	}

	if(allCons.size() > 0){
		matvar_t *cell_element = nullptr;

		// Step 2: Create the cell array
		size_t dims[2] = {allCons.size(), 1};
		pMatVar = Mat_VarCreate(pVarName, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
		if(pMatVar == nullptr){
			return pMatVar; 	// Can't save any data... exit
		}

		dims[0] = 1;
		for(unsigned int c = 0; c < allCons.size(); c++){
			// Get the constraint data
			std::vector<double> conData = allCons[c].getData();
			// Append the type and parent index to the beginning of the data vector: [type, parent_index, data0, data1, ...]
			int parentIx = -1;
			if(allCons[c].getAppType() != ConstraintApp_tp::APP_TO_ARC){
				parentIx = allCons[c].getAppType() == ConstraintApp_tp::APP_TO_NODE ? nodeIDMap.at(allCons[c].getID()) : segIDMap.at(allCons[c].getID());
			}else{
				parentIx = allCons[c].getID();
			}
			conData.insert(conData.begin(), static_cast<double>(parentIx));
			conData.insert(conData.begin(), static_cast<double>(to_underlying(allCons[c].getType())));

			dims[1] = conData.size();
			cell_element = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(conData.front()), 0);	// Do copy the data (it is freed at the end of the function!)
			if(cell_element != nullptr){
				Mat_VarSetCell(pMatVar, c, cell_element);
			}else{
				Mat_VarFree(pMatVar);
				throw Exception("BaseArcset::saveConstraints: Could not create cell array variable\n");
			}
		}
	}

	return pMatVar;
}//====================================================

/**
 *  @brief Create a matio variable for the node states
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeState(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	unsigned int stateSize = pSysData->getDynamicsModel()->getCoreStateSize();
	unsigned int numNodes = nodes.size();
	std::vector<double> posVel(stateSize*numNodes, NAN);

	for(unsigned int r = 0; r < numNodes; r++){
		std::vector<double> state = nodes[r].getState();

		if(state.size() < stateSize){
			char msg[256];
			sprintf(msg, "BaseArcset::createVar_NodeStates: "
				"State vector, length %zu, is less than core state size: %u",
				state.size(), stateSize);
			throw Exception(msg);
		}

		for(unsigned int c = 0; c < stateSize; c++){
			posVel[c*numNodes + r] = state[c];
		}
	}

	// Next, create a matlab variable for the state and save it to the file
	/*	Create a matlab variable. Arguments are:
	 *	const char *name 	- pVarName, the name of the variable
	 *	enum matio_classes 	- MAT_C_DOUBLE, Matlab double-precision variable class
	 *	enum matio_types 	- MAT_T_DOUBLE, Matlab IEEE 754 double precision data type
	 * 	int rank 			- 2 - the variable rank. Must be 2 or more; not really sure what this does
	 *	size_t dims 		- dims - the dimensions of the variable (e.g. matrix size) {rows, cols}
	 *	void *data 			- data - the variable we're saving. The algorithm assumes data is in column-major 
	 *							format
	 *	int opt 			- 0, or bit-wise OR of the following options:
	 *							MAT_F_DONT_COPY_DATA: just use the pointer to the data, don't copy it. 
	 *								Note that the pointer should not be freed until you are done with 
	 *								the matvar. The Mat_VarFree function will NOT free data that was
	 *								created with MAT_F_DONT_COPY_DATA, so free it yourself.
	 *							MAT_F_COMPLEX: specify that the data is complex
	 *							MAT_F_GLOBAL: make the matlab variable global
	 *							MAT_F_LOGICAL: this variable is a logical variable
	 */
	size_t dims[2] = {numNodes, stateSize};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(posVel.front()), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the node epoch times
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeEpoch(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	std::vector<double> allEpochs(nodes.size());

	for(unsigned int n = 0; n < nodes.size(); n++){
		allEpochs[n] = nodes[n].getEpoch();
	}
	
	size_t dims[2] = {allEpochs.size(), 1};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allEpochs[0]), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the node state derivatives
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeStateDeriv(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	unsigned int stateSize = pSysData->getDynamicsModel()->getCoreStateSize();

	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> deriv_colMaj(stateSize*nodes.size());

	for(unsigned int r = 0; r < nodes.size(); r++){
		std::vector<double> deriv(stateSize, NAN);
		try{
			deriv = nodes[r].getExtraParamVec(PARAMKEY_STATE_DERIV);
		}catch(Exception &e){
			// printErr("Unable to get acceleration vector for node %u\n", r);
		}

		for(unsigned int c = 0; c < deriv.size(); c++){
			deriv_colMaj[c*nodes.size() + r] = deriv[c];
		}
	}
	
	size_t dims[2] = {nodes.size(), stateSize};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(deriv_colMaj[0]), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the vector node extra parameters
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param varKey the key (i.e., the name) of the vector parameter
 * 	@param len number of elements in the extra parameter vector
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeExtraParamVec(std::string varKey, size_t len,
	Save_tp saveTp, const char *pVarName) const{

	(void) saveTp;
	// Get the specified coordinate
	std::vector<double> param(nodes.size()*len);
	for(unsigned int r = 0; r < nodes.size(); r++){
		// Save NAN (rather than un-allocated memory) if the node does not have the specified parameter
		std::vector<double> vec(len, NAN);

		try{
			vec = nodes[r].getExtraParamVec(varKey);
			for(unsigned int c = 0; c < vec.size(); c++){
				if(c >= len)
					break;

				param[c*nodes.size() + r] = vec[c];
			}
		}catch(Exception &e){
			// Save NAN (rather than un-allocated memory) if the node does not have the specified parameter
			for(unsigned int c = 0; c < len; c++)
				param[c*nodes.size() + r] = vec[c];
		}
	}

	size_t dims[2] = {nodes.size(), len};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(param[0]), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the scalar node extra parameters
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 	
 * 	@param varKey the unique key (i.e., name) of the scalar paremeter
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeExtraParam(std::string varKey, 
	Save_tp saveTp, const char *pVarName) const{

	(void) saveTp;
	// Get the specified coordinate
	std::vector<double> param(nodes.size());
	for(unsigned int r = 0; r < nodes.size(); r++){
		try{
			param[r] = nodes[r].getExtraParam(varKey);
		}catch(Exception &e){
			// Save NAN (rather than un-allocated memory) if the node does not have the specified parameter
			param[r] = NAN;
		}
	}

	size_t dims[2] = {nodes.size(), 1};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(param[0]), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the node control data
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_NodeCtrl(Save_tp saveTp, const char *pVarName) const{
	matvar_t *pMatVar = nullptr, *cell_element = nullptr;

	(void) saveTp;
	// Create the cell array
	size_t dims[2] = {nodes.size(), 1};
	pMatVar = Mat_VarCreate(pVarName, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
	if(pMatVar == nullptr){
		return pMatVar;	// Can't save any data... exit
	}

	dims[1] = 1;	 // only one column
	for(unsigned int n = 0; n < nodes.size(); n++){
		std::vector<double> ctrl;
		try{
			ctrl = nodes[n].getExtraParamVec(PARAMKEY_CTRL);
		}catch(Exception &e){}

		dims[0] = ctrl.size();

		// Save the data to an element of the cell array 		
		cell_element = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(ctrl.front()), 0);	// Do copy the data (it is freed at the end of the function!)
		if(cell_element != nullptr)
			Mat_VarSetCell(pMatVar, n, cell_element);
		else{
			Mat_VarFree(pMatVar);
			throw Exception("BaseArcset::saveNodeCtrl: "
				"Could not create cell array variable\n");
		}
	}

	return pMatVar;
}//====================================================

/**
 *  @brief Create a matio variable for the segment states
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_SegState(Save_tp saveTp, const char *pVarName) const{
	matvar_t *pMatVar = nullptr, *cell_element = nullptr;

	// Create the cell array
	size_t dims[2] = {segs.size(), 1};
	pMatVar = Mat_VarCreate(pVarName, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
	if(pMatVar == nullptr){
		return pMatVar;	// Can't save any data... exit
	}
	
	const unsigned int core_size = pSysData->getDynamicsModel()->getCoreStateSize();
	const unsigned int extra_size = pSysData->getDynamicsModel()->getExtraStateSize();

	for(unsigned int s = 0; s < segs.size(); s++){
		unsigned int ctrl_size = segs[s].getCtrlLaw() ? segs[s].getCtrlLaw()->getNumStates() : 0;
		unsigned int full_size = core_size + (core_size + ctrl_size)*(core_size + ctrl_size) + extra_size + ctrl_size;

		dims[1] = saveTp == Save_tp::SAVE_ALL ? full_size : core_size + ctrl_size;

		std::vector<double> segStates = segs[s].getStateVector();
		if(segStates.size() % full_size != 0){
			Mat_VarFree(pMatVar);
			char msg[256];
			sprintf(msg, "BaseArcset::saveSegState: "
				"segStates has length %zu; expecting full_size = %u\n"
				"\tRemainder = %zu THUS segStates is not multiple of full_size\n"
				"\tSegState width = %u",
				segStates.size(), full_size, segStates.size() % full_size,
				segs[s].getStateWidth());
			throw Exception(msg);
		}

		std::vector<double> segStates_trans;
		if(saveTp == Save_tp::SAVE_FRAME){
			dims[0] = 2;	// rows
			segStates_trans.assign(dims[0]*dims[1], 0);
			unsigned int lastRow = segStates.size()/full_size - 1;

			// Transpose data into column-major order
			for(unsigned int c = 0; c < dims[1]; c++){
				segStates_trans[c*dims[0] + 0] = segStates[0 + c];
				segStates_trans[c*dims[0] + 1] = segStates[lastRow*full_size + c];

			}
		}else{
			dims[0] = segStates.size()/full_size;	// rows
			segStates_trans.assign(dims[0]*dims[1], 0);

			// Transpose data into column-major order
			for(unsigned int r = 0; r < dims[0]; r++){
				for(unsigned int c = 0; c < dims[1]; c++){
					segStates_trans[c*dims[0] + r] = segStates[r*full_size + c];
				}
			}
		}

		// Save the data to an element of the cell array 		
		cell_element = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(segStates_trans.front()), 0);	// Do copy the data (it is freed at the end of the function!)
		if(cell_element != nullptr)
			Mat_VarSetCell(pMatVar, s, cell_element);
		else{
			Mat_VarFree(pMatVar);
			throw Exception("BaseArcset::saveSegStates: Could not create cell array variable\n");
		}
	}

	return pMatVar;
}//====================================================

/**
 *  @brief Create a matio variable for the segment times
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_SegTime(Save_tp saveTp, const char *pVarName) const{
	matvar_t *pMatVar = nullptr, *cell_element = nullptr;

	// Create the cell array
	size_t dims[2] = {segs.size(), 1};
	pMatVar = Mat_VarCreate(pVarName, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
	if(pMatVar == nullptr){
		return pMatVar;	// Can't save any data... exit
	}

	dims[1] = 1;	 // only one column
	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<double> segTimes = segs[s].getTimeVector();

		if(saveTp == Save_tp::SAVE_FRAME){
			// Only save the first and last time
			dims[0] = 2;
			segTimes.erase(segTimes.begin()+1, segTimes.end()-1);
		}else{
			dims[0] = segTimes.size();
		}

		// Save the data to an element of the cell array 		
		cell_element = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(segTimes.front()), 0);	// Do copy the data (it is freed at the end of the function!)
		if(cell_element != nullptr)
			Mat_VarSetCell(pMatVar, s, cell_element);
		else{
			Mat_VarFree(pMatVar);
			throw Exception("BaseArcset::saveSegStates: Could not create cell array variable\n");
		}
	}

	return pMatVar;
}//====================================================

/**
 *  @brief Create a matio variable for the segment times-of-flight
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_SegTOF(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	std::vector<double> allTOFs(segs.size());

	for(unsigned int s = 0; s < segs.size(); s++){
		allTOFs[s] = segs[s].getTOF();
	}
	
	size_t dims[2] = {allTOFs.size(), 1};
	return Mat_VarCreate(pVarName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allTOFs[0]), 0);	// Do copy the data (it is freed at the end of the function!)
}//====================================================

/**
 *  @brief Create a matio variable for the segment STMs
 *  @details The STM is copied from each segment state vector, thus,
 *  the STM represents the evolution of each individual segment 
 *  regardless of whether or not the arcset has been set to store
 *  cumulative STMs via setSTM_cumulative().
 *  
 *  The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_SegSTM(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	matvar_t *pMatVar = nullptr, *cell_element = nullptr;

	size_t dims[2] = {segs.size(), 1};
	pMatVar = Mat_VarCreate(pVarName, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
	if(pMatVar == nullptr){
		return pMatVar;	// can't save any data... exit
	}

	unsigned int count = 0;
	unsigned int core_dim = pSysData->getDynamicsModel()->getCoreStateSize(), ctrl_dim = 0;
	for(const Segment &seg : segs){
		ctrl_dim = seg.getCtrlLaw() ? seg.getCtrlLaw()->getNumStates() : 0;

		MatrixXRd P = seg.getSTM_fromStates(core_dim, ctrl_dim).transpose();
		dims[0] = P.cols();
		dims[1] = P.rows();

		cell_element = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, P.data(), 0);	// Do copy the data (it is freed at the end of the function!)
		if(cell_element != nullptr){
			Mat_VarSetCell(pMatVar, count++, cell_element);
		}else{
			Mat_VarFree(pMatVar);
			throw Exception("BaseArcset::saveSegSTMs: Could not create cell array variable\n");
		}
	}

	return pMatVar;
}//====================================================

/**
 *  @brief Create a matio variable for the segment control law data
 *  @details The data is copied into a matvar_t pointer, which is
 *  allocated on the stack
 * 
 * 	@param saveTp describes how much data to save
 *  @param pVarName variable name
 *  @return pointer to the matio variable (must be freed by MatVar_Free())
 */
matvar_t* BaseArcset::createVar_SegCtrlLaw(Save_tp saveTp, const char *pVarName) const{
	(void) saveTp;
	size_t struct_dims[] = {segs.size(),1};
	
	unsigned int nfields = 3;
	const char *fieldnames[3] = {"Type", "NumStates", "Params"};

	matvar_t *pMatVar = Mat_VarCreateStruct(pVarName, 2, struct_dims, fieldnames, nfields);
	if(pMatVar == nullptr){
		printErr("ControlLaw::createMatStruct: Error creating ControlLaw structure variable\n");
		return pMatVar;
	}

	matvar_t *field = nullptr;

	for(unsigned int s = 0; s < segs.size(); s++){
		ControlLaw *pLaw = segs[s].getCtrlLaw();

		// Initialize variables to represent the no control case
		unsigned int id = 0, numStates = 0;
		std::vector<double> params;
		
		// Update values if a control law has been implemented
		if(pLaw){
			id = pLaw->getType();
			numStates = pLaw->getNumStates();
			params = pLaw->getParams();
		}

		// Save structure fields
		size_t field_dims[] = {1,1};
		field = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, field_dims, &id, 0);	// Do copy the data (it is freed at the end of the function!)
		Mat_VarSetStructFieldByName(pMatVar, fieldnames[0], s, field);
		// No need to free `field`; Mat_VarSetStructFieldByName frees the variable
		
		field = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, field_dims, &numStates, 0);	// Do copy the data (it is freed at the end of the function!)
		Mat_VarSetStructFieldByName(pMatVar, fieldnames[1], s, field);

		field_dims[0] = params.size();
		field = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, field_dims, &(params.front()), 0);	// Do copy the data (it is freed at the end of the function!)
		Mat_VarSetStructFieldByName(pMatVar, fieldnames[2], s, field);
	}

	return pMatVar;
}//====================================================


bool BaseArcset::readVar_LinkTable(matvar_t *pVar){
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_LinkTable: Could not read link table from file");
		return false;
	}else{
		unsigned int numSegs = pVar->dims[0];

		if(pVar->dims[1] != 4){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_LinkTable: Link table width is not the expected dimension (4)");
		}

		if(pVar->class_type == MAT_C_INT32 && pVar->data_type == MAT_T_INT32){
			int *data = static_cast<int *>(pVar->data);

			if(data != nullptr){
				// Step 0: Clear all variables
				nodes.clear();
				segs.clear();
				nodeIDMap.clear();
				segIDMap.clear();
				nextNodeID = 0;
				nextSegID = 0;

				// Step 1: Count the unique node indices to get a count of the number of nodes
				std::vector<int> nodeIxs;
				for(unsigned int s = 0; s < numSegs; s++){
					int oldOriginIx = data[1*numSegs + s];
					int oldTerminusIx = data[2*numSegs + s];

					if(oldOriginIx != Linkable::INVALID_ID && std::find(nodeIxs.begin(), nodeIxs.end(), oldOriginIx) == nodeIxs.end())
						nodeIxs.push_back(oldOriginIx);

					if(oldTerminusIx != Linkable::INVALID_ID && std::find(nodeIxs.begin(), nodeIxs.end(), oldTerminusIx) == nodeIxs.end())
						nodeIxs.push_back(oldTerminusIx);
				}

				if(nodeIxs.size() == 0)
					throw Exception("BaseArcset::readVar_LinkTable: No nodes linked from segments... not built to handle this");

				unsigned int numNodes = *std::max_element(nodeIxs.begin(), nodeIxs.end()) + 1;

				// Step 2: Create the nodes; index = id in this case
				for(unsigned int n = 0; n < numNodes; n++){
					addNode(Node());	// updates the nextNodeID counter and the nodeIDMap
				}

				// Step 3: Create the segments: Use the link indices from the link table to get
				// the correct link network.
				for(unsigned int s = 0; s < numSegs; s++){
					// The fourth column lists +1 for forward-time and -1 for reverse-time
					// Updates the nextSegID counter and the segIDMap
					addSeg(Segment(data[1*numSegs + s], data[2*numSegs + s], data[3*numSegs + s]));
				}

				// std::cout << "tempNodeIDMap:\n";
				// for(auto it = tempNodeIDMap.cbegin(); it != tempNodeIDMap.cend(); ++ it){
				// 	std::cout << it->first << " -> " << it->second << std::endl;
				// }
			}
		}
	}
	return true;
}//====================================================

/**
 *  @brief Read constraints from a matio variable
 * 
 *  @param pVar a pointer to a cell array variable that stores the constraints
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_Constraints(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;

	if(pVar != nullptr){
		unsigned int numCons = pVar->dims[0];

		if(pVar->class_type != MAT_C_CELL || pVar->data_type != MAT_T_CELL){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_Constraints: Constraint cell array variable is not stored as a cell array");
		}

		matvar_t **cell_elements = static_cast<matvar_t **>(pVar->data);

		for(unsigned int c = 0; c < numCons; c++){
			if(cell_elements[c] != nullptr && 
				cell_elements[c]->class_type == MAT_C_DOUBLE && cell_elements[c]->data_type == MAT_T_DOUBLE){

				unsigned int width = cell_elements[c]->dims[1];
				
				if(width < 2 || cell_elements[c]->dims[0] == 0){
					Mat_VarFree(pVar);
					throw Exception("BaseArcset::readVar_Constraints: Constraint data has less than 2 (the minimum) data values");
				}

				double *data = static_cast<double *>(cell_elements[c]->data);
				if(data != nullptr){
					Constraint_tp tp = static_cast<Constraint_tp>(static_cast<int>(data[0]));
					int parentIx = static_cast<int>(data[1]);
					addConstraint(Constraint(tp, parentIx, data+2, width-2));
				}
			}
		}
		return true;
	}else{
		return false;
	}
}//====================================================

/**
 *  @brief Read node states from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node states
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeState(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;

	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeState: Could not read state data vector");
		return false;
	}else{
		unsigned int core_dim = pSysData->getDynamicsModel()->getCoreStateSize();
		unsigned int numSteps = pVar->dims[0];
		
		if(nodes.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeState: Node vector has not been initialized!");
		}

		if(numSteps != nodes.size()){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeState: State vector has a different size than the initialized node vector");
		}

		if(pVar->dims[1] != core_dim){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeState: Incompatible data file: State width is not consistent with DynamicalModel definition.");
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);

			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					std::vector<double> state;
					for(unsigned int s = 0; s < core_dim; s++){
						state.push_back(data[s*numSteps + i]);
					}

					nodes[i].setState(state);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeState: Incompatible data file: unsupported data type/class");
		}
		return true;
	}
}//====================================================

/**
 *  @brief Read node state derivatives from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node state derivatives
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeStateDeriv(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;

	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeStateDeriv: Could not read data vector");
		return false;
	}else{
		unsigned int core_dim = pSysData->getDynamicsModel()->getCoreStateSize();
		unsigned int numSteps = pVar->dims[0];
		
		if(nodes.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeStateDeriv: Node vector has not been initialized!");
		}

		if(numSteps != nodes.size()){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeStateDeriv: Derivative vector has a different size than the initialized node vector");
		}

		if(pVar->dims[1] != core_dim){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeStateDeriv: Incompatible data file: Vector width is not equal to the state size.");
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);

			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					std::vector<double> deriv(core_dim, 0);

					for(unsigned int s = 0; s < core_dim; s++){
						deriv[s] = data[s*numSteps + i];
					}
					nodes[i].setExtraParamVec(PARAMKEY_STATE_DERIV, deriv);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeStateDeriv: Incompatible data file: unsupported data type/class");
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read node epochs from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node epochs
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeEpoch(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;

	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeEpoch: Could not read data vector");
		return false;
	}else{
		unsigned int numSteps = pVar->dims[0];

		if(nodes.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeEpoch: Node vector has not been initialized");
		}

		if(numSteps != nodes.size()){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeEpoch: Epoch vector has different size than the initialized node evctor");
		}

		if(pVar->dims[1] != 1){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeEpoch: Incompatible data file: Epoch vector has more than one column");
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);

			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					nodes[i].setEpoch(data[i]);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeEpoch: Incompatible data file: unsupported data type or class");
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read node control data from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node control data
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeCtrl(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeCtrl: Could not read state data vector");
		return false;
	}else{
		if(pVar->class_type == MAT_C_CELL && pVar->data_type == MAT_T_DOUBLE){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeCtrl: Node control variable is not a cell array.");
		}

		unsigned int numNodes = pVar->dims[0];
		
		if(nodes.size() != numNodes){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeCtrl: Node vector has been initialized to a different size than the file has data for");
		}

		matvar_t **cell_elements = static_cast<matvar_t **>(pVar->data);

		for(unsigned int n = 0; n < numNodes; n++){
			if(cell_elements[n]){
				if(cell_elements[n]->class_type == MAT_C_DOUBLE && cell_elements[n]->data_type == MAT_T_DOUBLE){
					unsigned int numSteps = cell_elements[n]->dims[0];	
					if(numSteps > 0){
						double *data = static_cast<double *>(cell_elements[n]->data);
						std::vector<double> ctrl(data, data+numSteps);
						nodes[n].setExtraParamVec(PARAMKEY_CTRL, ctrl);
					}
				}else{
					Mat_VarFree(pVar);
					throw Exception("BaseArcset::readVar_NodeCtrl: Cell element is not a double array.");
				}
			}
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read node extra parameter data from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node extra parameters
 *  @param varKey the key (i.e., name) of the extra parameter vector
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeExtraParam(matvar_t *pVar, std::string varKey, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeExtraParam: Could not read data vector");
		return false;
	}else{
		unsigned int numSteps = pVar->dims[0];
		
		if(nodes.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeExtraParam: Step vector has not been initialized!");
		}

		if(pVar->dims[1] != 1){
			char message[128];
			const char * pVarName = pVar->name ? pVar->name : "NULL";
			sprintf(message, "BaseArcset::readVar_NodeExtraParam: Incompatible data file: %s width is not %d", pVarName, 1);
			Mat_VarFree(pVar);
			throw Exception(message);
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);
			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					nodes[i].setExtraParam(varKey, data[i]);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeExtraParam: Incompatible data file: unsupported data type/class");
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read node extra parameter vectors from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the node extra parameter vectors
 *  @param varKey the key (i.e., name) of the extra parameter vector
 *  @param len the length of the extra parameter vector
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_NodeExtraParamVec(matvar_t *pVar, std::string varKey, size_t len, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_NodeExtraParamVec: Could not read data vector");
		return false;
	}else{
		unsigned int numSteps = pVar->dims[0];
		
		if(nodes.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeExtraParamVec: Step vector has not been initialized!");
		}

		if(pVar->dims[1] != len){
			char msg[128];
			const char *varName = pVar->name ? pVar->name : "NULL";
			sprintf(msg, "BaseArcset::readVar_NodeExtraParamVec: Incompatible data file: %s width is not %zu", varName, len);
			Mat_VarFree(pVar);
			throw Exception(msg);
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);
			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					std::vector<double> vec(len,0);
					for(unsigned int c = 0; c < len; c++){
						vec[c] = data[c*numSteps + i];
					}
					nodes[i].setExtraParamVec(varKey, vec);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_NodeExtraParamVec: Incompatible data file: unsupported data type/class");
		}
		return true;
	}
}//====================================================

/**
 *  @brief Read segment state data from a matio variable
 * 
 *  @param pVar a pointer to a cell array variable that stores the segment states
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_SegState(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_SegState: "
			"Could not read state data vector");
		return false;
	}else{
		unsigned int coreDim = pSysData->getDynamicsModel()->getCoreStateSize();
		unsigned int extraDim = pSysData->getDynamicsModel()->getExtraStateSize();
		unsigned int numSegs = pVar->dims[0];
		
		if(segs.size() != numSegs){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegState:\n"
				"\tSegment vector has been initialized to a different size than "
				"the file has data for");
		}

		if(pVar->class_type != MAT_C_CELL || pVar->data_type != MAT_T_CELL){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegState: "
				"Segment state variable is not a cell array.");
		}

		matvar_t **cell_elements = static_cast<matvar_t **>(pVar->data);

		for(unsigned int s = 0; s < numSegs; s++){
			if(cell_elements[s] != nullptr && 
				cell_elements[s]->class_type == MAT_C_DOUBLE && 
				cell_elements[s]->data_type == MAT_T_DOUBLE){

				unsigned int ctrlDim = segs[s].getCtrlLaw() ? \
					segs[s].getCtrlLaw()->getNumStates() : 0;

				// If all data is read, this is the expected width
				unsigned int expectedWidth = coreDim + 
					(coreDim+ctrlDim)*(coreDim+ctrlDim) + extraDim + ctrlDim;
				std::vector<double> fillerStates;


				unsigned int width = cell_elements[s]->dims[1];
				unsigned int numSteps = cell_elements[s]->dims[0];
				double *data = static_cast<double *>(cell_elements[s]->data);

				if(expectedWidth < width){
					char msg[128];
					sprintf(msg, "BaseArcset::readVar_SegState: "
						"expected width, %u, is less than data width, %u. "
						"Cannot proceed",
						expectedWidth, width);
					throw Exception(msg);
				}

				// If some states are missing from the data file (e.g., if a 
				// conservative save method was used) append a "filler state" 
				// full of zeros
				if(expectedWidth - width > 0){
					fillerStates.assign(expectedWidth - width, 0);
				}

				segs[s].setStateWidth(expectedWidth);
				if(data != nullptr){
					for(unsigned int i = 0; i < numSteps; i++){
						std::vector<double> state(width);
						for(unsigned int c = 0; c < width; c++)
							state[c] = data[c*numSteps + i];

						segs[s].appendState(state);

						if(expectedWidth - width > 0)
							segs[s].appendState(fillerStates);
					}
				}
			}else{
				Mat_VarFree(pVar);
				throw Exception("BaseArcset::readVar_SegState: "
					"Cell element is not a double array.");
			}
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read segment time data from a matio variable
 * 
 *  @param pVar a pointer to a cell array variable that stores the segment times
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_SegTime(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;

	if(pVar == nullptr){
		printErr("BaseArcset::readSegStatesFromMat: Could not read state data vector");
		return false;
	}else{
		unsigned int numSegs = pVar->dims[0];
		
		if(segs.size() != numSegs){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readSegStatesFromMat: Segment vector has been initialized to a different size than the file has data for");
		}

		if(pVar->class_type != MAT_C_CELL || pVar->data_type != MAT_T_CELL){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readSegStatesFromMat: Segment state variable is not a cell array.");
		}

		matvar_t **cell_elements = static_cast<matvar_t **>(pVar->data);

		for(unsigned int s = 0; s < numSegs; s++){
			if(cell_elements[s]->class_type == MAT_C_DOUBLE && cell_elements[s]->data_type == MAT_T_DOUBLE){
				unsigned int numSteps = cell_elements[s]->dims[0];	
				double *data = static_cast<double *>(cell_elements[s]->data);
				std::vector<double> times(data, data+numSteps);
				segs[s].setTimeVector(times);
			}else{
				Mat_VarFree(pVar);
				throw Exception("BaseArcset::readSegStatesFromMat: Cell element is not a double array.");
			}
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read segment control law data from a matio variable
 * 
 *  @param pVar a pointer to a structure array variable that stores the segment 
 *  control laws
 *  @param refLaws Reference to a vector of ControlLaw pointers. As control laws 
 *  are read from the matio variable, unique control laws are constructed and 
 *  allocated on the stack. The user must manually delete the ControlLaw objects 
 *  to avoid memory leaks.
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_SegCtrlLaw(matvar_t *pVar, 
	std::vector<ControlLaw*> &refLaws, Save_tp saveTp){

	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_SegCtrlLaw: Could not read data vector");
		return false;
	}else{
		const unsigned int fieldsPerController = 3;
		unsigned int numStructs = (pVar->dims[0])*(pVar->dims[1]);

		if(segs.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegCtrlLaw: "
				"Segment vector has not been initialized");
		}

		if(numStructs != segs.size()){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegCtrlLaw: "
				"Control Law vector has different size"
				" than the initialized segment evctor");
		}

		if(pVar->class_type == MAT_C_STRUCT && pVar->data_type == MAT_T_STRUCT){
			matvar_t **fields = static_cast<matvar_t **>(pVar->data);
			
			if(fields){

				// Loop through the controllers; should be same as number of 
				// segs (checked above)
				for(unsigned int s = 0; s < numStructs; s++){
					unsigned int id = 0;
					std::vector<double> params;

					// Loop through fields within one controller structure
					for(unsigned int f = 0; f < fieldsPerController; f++){

						matvar_t *oneField = fields[s*fieldsPerController + f];
						if(oneField){
							switch(f){
								case 0: // lawType
								{
									if(oneField->class_type == MAT_C_INT32 && 
										oneField->data_type == MAT_T_INT32){
										
										unsigned int *intData = \
											static_cast<unsigned int *>(
												oneField->data);

										if(intData)
											id = intData[0];
										else{
											Mat_VarFree(pVar);
											throw Exception("BaseArcset::"
												"readVar_SegCtrlLaw: controller "
												"ID data is nullptr");
										}
									}else{
										Mat_VarFree(pVar);
										throw Exception("BaseArcset::"
											"readVar_SegCtrlLaw: controller ID "
											"field has wrong class type or "
											"data type");
									}
									break;
								}//-----------------------------------
								case 1: // numStates
									// No need to read this in; will be 
									// populated when ControlLaw is constructed
									break;
								//-----------------------------------
								case 2: // params
								{
									if(oneField->class_type == MAT_C_DOUBLE && 
										oneField->data_type == MAT_T_DOUBLE){
										
										double *doubleData = 
											static_cast<double *>(oneField->data);

										if(doubleData){
											unsigned int len = 
												oneField->dims[0] * 
												oneField->dims[1];
											params.insert(params.begin(), 
												doubleData, doubleData + len);	
										}else{
											Mat_VarFree(pVar);
											throw Exception("BaseArcset::readVar_SegCtrlLaw: "
												"controller params data is nullptr");
										}
									}else{
										Mat_VarFree(pVar);
										throw Exception("BaseArcset::readVar_SegCtrlLaw: "
											"controller params field has wrong "
											"class type or data type");
									}
									break;
								}//-----------------------------------
							}// End of Switch/Case
						}
					}// End of loop through controller fields
					
					// Only proceed if the control is nontrivial
					if(id != ControlLaw::NO_CTRL){

						if(id < (1<<11)){
							// likely saved with old IDs
							ControlLaw::convertID(id);
						}
						// Allocate a new control law on the stack; by using a 
						// function in the DynamicsModel, we ensure that the 
						// system-specific derived class is constructed rather 
						// than the base class ControlLaw.
						ControlLaw *newLaw = 
							pSysData->getDynamicsModel()->createControlLaw(id, 
								params);

						// Check to see if the controller has been loaded already
						bool foundDuplicate = false;
						for(auto &law : refLaws){
							if(law && *law == *newLaw){
								delete(newLaw);
								newLaw = law;
								foundDuplicate = true;
								break;
							}
						}

						// Save the pointer to return back to the top level for
						// further use
						if(!foundDuplicate)
							refLaws.push_back(newLaw);

						// Assign the control law to the segment
						segs[s].setCtrlLaw(newLaw);
					}
				}// End loop through structs
			}else{
				Mat_VarFree(pVar);
				throw Exception("BaseArcset::readVar_SegCtrlLaw: Structure array is nullptr");
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegCtrlLaw: Incompatible data file: unsupported data type or class");
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read segment TOF data from a matio variable
 * 
 *  @param pVar a pointer to a variable that stores the segment TOFs
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_SegTOF(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_SegTOF: Could not read data vector");
		return false;
	}else{
		unsigned int numSteps = pVar->dims[0];

		if(segs.size() == 0){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegTOF: Segment vector has not been initialized");
		}

		if(numSteps != segs.size()){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegTOF: TOF vector has different size than the initialized segment evctor");
		}

		if(pVar->dims[1] != 1){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegTOF: Incompatible data file: TOF vector has more than one column");
		}

		if(pVar->class_type == MAT_C_DOUBLE && pVar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pVar->data);

			if(data != nullptr){
				for(unsigned int i = 0; i < numSteps; i++){
					segs[i].setTOF(data[i]);
				}
			}
		}else{
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegTOF: Incompatible data file: unsupported data type or class");
		}

		return true;
	}
}//====================================================

/**
 *  @brief Read segment STM data from a matio variable
 * 
 *  @param pVar a pointer to a cell array variable that stores the segment STMs
 *  @param saveTp Describes the amount of detail that was saved to the file
 * 
 *  @return whether or not the matio variable needs to be freed
 *  @throws Exception if there are any errors while importing the data
 */
bool BaseArcset::readVar_SegSTM(matvar_t *pVar, Save_tp saveTp){
	(void) saveTp;
	if(pVar == nullptr){
		printErr("BaseArcset::readVar_SegSTM: Could not read STM data vector");
		return false;
	}else{
		unsigned int numSegs = pVar->dims[0];
		
		if(segs.size() != numSegs){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegSTM: Segment vector has been initialized to a different size than the file has data for");
		}

		if(pVar->class_type != MAT_C_CELL || pVar->data_type != MAT_T_CELL){
			Mat_VarFree(pVar);
			throw Exception("BaseArcset::readVar_SegSTM: Segment STM variable is not a cell array.");
		}

		matvar_t **cell_elements = static_cast<matvar_t **>(pVar->data);

		for(unsigned int s = 0; s < numSegs; s++){
			if(cell_elements[s]->class_type == MAT_C_DOUBLE && cell_elements[s]->data_type == MAT_T_DOUBLE){				
				double *data = static_cast<double *>(cell_elements[s]->data);

				MatrixXRd P = Eigen::Map<MatrixXRd>(data, cell_elements[s]->dims[1], cell_elements[s]->dims[0]);
				segs[s].setSTM(P.transpose());
			}else{
				Mat_VarFree(pVar);
				throw Exception("BaseArcset::readVar_SegSTM: Cell element is not a double array.");
			}
		}

		return true;
	}
}//====================================================

/**
 *  @brief Initialize the vectors of node and segment objects from a *.mat file
 *  @details DEPRECATED; replaced by readLinkTable().
 *  
 *  THIS FUNCTION MUST BE THE FIRST READ_DATA-TYPE FUNCTION CALLED because
 *  it clears the vectors and then initializes them by calculating the number
 *  of steps in the arcset object from the state vector. Individual nodes and segments are
 *  able to be called by index after this, though they will not contain data
 *  until another function is called to populate the data fields with values from 
 *  the *.mat file
 * 
 *  @param pMatFile pointer to an open matlab data file
 *  @param pVarName the name of a variable that has as many rows as there are
 *  steps along the data object. Valid variables typically include the time vector,
 *  state matrix, or acceleration matrix
 *  @throws Exception if the state vector variable cannot be read from the data file
 */
void BaseArcset::initNodesSegsFromMat(mat_t *pMatFile, const char* pVarName){
	matvar_t *pStateMat = Mat_VarRead(pMatFile, pVarName);
	if(pStateMat == nullptr){
		throw Exception("BaseArcset::initNodesSegsFromMat: Could not read state data vector");
	}else{
		unsigned int numSteps = pStateMat->dims[0];
		nodes.clear();
		segs.clear();
		nodeIDMap.clear();
		segIDMap.clear();
		
		// Create a set of nodes and segments all linked together in linear time
		for(unsigned int i = 0; i < numSteps; i++){
			Node node;
			node.setID(i);

			if(i > 0)
				node.addLink(i-1);
			if(i < numSteps-1)
				node.addLink(i);

			nodes.push_back(node);
			nodeIDMap[i] = i;
			
			if(i > 0){
				Segment seg(i-1, i, NAN);
				seg.setID(i-1);
				segs.push_back(seg);
				segIDMap[i-1] = i-1;
			}
		}
	}
	Mat_VarFree(pStateMat);
}//======================================================

}// END of Astrohelion namespace