/**
 *  @file tpat_base_arcset.cpp
 *	@brief Data object that stores information about an integrated arc
 *	
 *	@author Andrew Cox
 *	@version April 28, 2016
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

#include "tpat_base_arcset.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

#include <algorithm>
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Constructor (requires system data object)
 *	@param sys a pointer to a system data object that describes
 *	the system this trajectory is integrated in
 */
tpat_base_arcset::tpat_base_arcset(const tpat_sys_data *sys) : sysData(sys){}

/**
 *	@brief Copy constructor
 *	@param d an arcset reference
 */
tpat_base_arcset::tpat_base_arcset(const tpat_base_arcset &d) : sysData(d.sysData){
	copyMe(d);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_base_arcset::~tpat_base_arcset(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Set this object equal to another
 *	@param d an arcset reference
 *	@return a reference to this arcset object
 */
tpat_base_arcset& tpat_base_arcset::operator =(const tpat_base_arcset &d){
	copyMe(d);
	sysData = d.sysData;
	return *this;
}//====================================================

/**
 *  @brief Sum two arcset objects.
 *  @details This function returns <tt>result</tt> = <tt>lhs</tt> + <tt>rhs</tt>; 
 *  Both <tt>lhs</tt> and <tt>rhs</tt> are copied and sorted into chronological 
 *  order. The <tt>rhs</tt> is then appended to the end of <tt>lhs</tt> with no
 *  time-of-flight between the end of <tt>lhs</tt> and the beginning of <tt>rhs</tt>.
 *  This will result in the deletion of the first node of <tt>rhs</tt> in the 
 *  concatenated set. To avoid losing data, make sure the first node of <tt>rhs</tt>
 *  is equivalent to the final node (chronologically) of <tt>lhs</tt>.
 * 
 *  @param lhs pointer to an arcset object
 *  @param rhs pointer to an arcset object
 *  @param result a pointer to the an arcset object that will store the result of
 *  the summing operation.
 */
void tpat_base_arcset::sum(const tpat_base_arcset *lhs, const tpat_base_arcset *rhs, tpat_base_arcset *result){
	baseArcsetPtr lhs_cpy = lhs->clone();
	baseArcsetPtr rhs_cpy = rhs->clone();

	lhs_cpy->putInChronoOrder();
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
 *  @see tpat_constraint::tpat_conApp_tp
 *  @throws tpat_exception if node or segment ID is out of bounds, or if the
 *  constraint application type is unknown
 */
void tpat_base_arcset::addConstraint(tpat_constraint con){
	int id = con.getID();
	switch(con.getAppType()){
		case tpat_constraint::APP_TO_NODE:
			if(id < 0 || id >= (int)(nodeIDMap.size()))
				throw tpat_exception("tpat_base_arcset::addConstraint: Node ID out of bounds");

			nodes[nodeIDMap[id]].addConstraint(con);
			break;
		case tpat_constraint::APP_TO_SEG:
			if(id < 0 || id >= (int)(segIDMap.size()))
				throw tpat_exception("tpat_base_arcset::addConstraint: Segment ID out of bounds");

			segs[segIDMap[id]].addConstraint(con);
			break;
		case tpat_constraint::APP_TO_ARC:
			cons.push_back(con);
			break;
		default:
			throw tpat_exception("tpat_base_arcset::addConstraint: Constraint application type is unknown");
	}
}//=====================================================

/**
 *  @brief Add a node to this data object
 *  @details A unique key is assigned to the node when it is added
 * 
 *  @param n the node to add
 *  @return the ID assigned to the node
 */
int tpat_base_arcset::addNode(tpat_node n){
	n.clearLinks();			// Cannot have any links coming in
	n.setID(nextNodeID);
	nodes.push_back(n);
	nodeIDMap.push_back(nodes.size()-1);

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
 *  @throws tpat_exception if adding the segment will result in a time direction
 *  collision or a parallel structure. tpat_exception is also thrown if the segment
 *  is linked to a non-existant node other than the placeholder 
 *  tpat_linkable::INVALID_ID
 */
int tpat_base_arcset::addSeg(tpat_segment s){
	s.setID(nextSegID);
	
	if(s.getOrigin() == tpat_linkable::INVALID_ID)
		throw tpat_exception("tpat_base_arcset::addSeg: Segment must have a valid origin node\n");
	
	bool foundValidLink = false;
	for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
		// Get the ID of one of the nodes this segment is linked to
		int linkedNodeID = s.getLink(i);
		if(linkedNodeID != tpat_linkable::INVALID_ID){

			if(linkedNodeID < 0 || linkedNodeID >= (int)(nodeIDMap.size()))
				throw tpat_exception("tpat_base_arcset::addSeg: Segment has link to an invalid node ID");

			// Get the index of that node within the storage array
			int linkedNodeIx = nodeIDMap[linkedNodeID];
			if(linkedNodeIx != tpat_linkable::INVALID_ID){

				// See if the node is linked to any other segments
				tpat_node *theNode = &(nodes[linkedNodeIx]);
				int secondaryLinks = 0;
				for(int j = 0; j < tpat_linkable::NUM_LINKS; j++){
					// Check to make sure the link is to a real segment
					if(theNode->getLink(j) != tpat_linkable::INVALID_ID){
						int nearSegIx = segIDMap[theNode->getLink(j)];
						// Make sure the segment is real
						if(nearSegIx != tpat_linkable::INVALID_ID){
							secondaryLinks++;
							// If the node is linked to another segment, get that segment and compare it to this one
							tpat_segment *nearSeg = &(segs[nearSegIx]);
							bool sameLinkType = nearSeg->getLink(i) == linkedNodeID;
							bool sameTimeDir = nearSeg->getTOF()*s.getTOF() > 0;

							if(sameLinkType && i == tpat_segment::TERM_IX){
								throw tpat_exception("tpat_base_arcset::addSeg: would create node with two terminating segments");
							}else if(sameLinkType && sameTimeDir){
								// either a time collision (both terminate) or parallel structure (both originate)
								printf("Nearby segment w/ ID %d originates at node %d and terminates at node %d\n", nearSeg->getID(), nearSeg->getOrigin(), nearSeg->getTerminus());
								printf("The new seg (ID %d) originates at node %d and terminates at node %d\n", s.getID(), s.getOrigin(), s.getTerminus());
								throw tpat_exception("tpat_base_arcset::addSeg: either time collision or parallel structure");
							}else if(!sameLinkType && !sameTimeDir){
								// parallel structure
								throw tpat_exception("tpat_base_arcset::addSeg: parallel structure");
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
				throw tpat_exception("tpat_base_arcset::addSeg: Linked node is not part of this arcset object");
			}
		}
	}

	if(foundValidLink){
		segs.push_back(s);
		segIDMap.push_back(segs.size()-1);
		return nextSegID++;	
	}else{
		return tpat_linkable::INVALID_ID;
	}
}//===========================================

/**
 *  @brief Append an arcset object (i.e., a set of nodes and segments) to this one
 * 
 *  @param arcset a pointer to the arcset derivative object to append
 *  @param linkTo_ID the ID of the node in *this* arcset object to link to
 *  @param linkFrom_ID the ID of the node in <tt>set</tt> to link from
 *  @param tof time-of-flight between linkFrom_ID to linkTo_ID; if set to zero, 
 *  it is assumed that the two nodes are identical, and the original (the one in *this* object)
 *  will be retained, the other deleted, and segments rerouted accordingly.
 *  
 *  @return the ID of a new segment that links the old and new arcset objects
 *  @throws tpat_exception if either ID is out of bounds
 *  @throws tpat_exception if one or both of the identifies nodes does not have
 *  a free link slot
 */
int tpat_base_arcset::appendSetAtNode(const tpat_base_arcset *arcset, int linkTo_ID, int linkFrom_ID, double tof){
	// First, check to make sure the specified nodes are valid
	if(linkTo_ID < 0 || linkTo_ID >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::appendSetAtNode: linkTo_ID is out of bounds");

	// Create a copy so we don't affect the original
	baseArcsetPtr set = arcset->clone();

	tpat_node linkTo_node = nodes[nodeIDMap[linkTo_ID]];
	tpat_node linkFrom_node = set->getNode(linkFrom_ID);		// Will do its own index checks

	// Both nodes must have one "open port"
	if(!linkTo_node.isLinkedTo(tpat_linkable::INVALID_ID) || !linkFrom_node.isLinkedTo(tpat_linkable::INVALID_ID))
		throw tpat_exception("tpat_base_arcset::appendSetAtNode: specified nodes are not both open to a new link");

	// Determine if linkTo_node is the origin or terminus of a segment
	// printf("linkToNode has links [%d, %d]\n", linkTo_node.getLink(0), linkTo_node.getLink(1));
	// printf("Choosing segment (ID %d)\n", linkTo_node.getLink(0) == tpat_linkable::INVALID_ID ? linkTo_node.getLink(1) : linkTo_node.getLink(0));
	tpat_segment linkTo_seg = getSeg(linkTo_node.getLink(0) == tpat_linkable::INVALID_ID ? linkTo_node.getLink(1) : linkTo_node.getLink(0));
	bool linkTo_isOrigin = linkTo_seg.getOrigin() == linkTo_node.getID();
	tpat_segment linkFrom_seg = set->getSeg(linkFrom_node.getLink(0) == tpat_linkable::INVALID_ID ? linkFrom_node.getLink(1) : linkFrom_node.getLink(0));
	bool linkFrom_isOrigin = linkFrom_seg.getOrigin() == linkFrom_node.getID();

	if(!linkTo_isOrigin && !linkFrom_isOrigin)
		throw tpat_exception("tpat_base_arcset::appendSetAtNode: neither node is an origin; cannot create segment between them\n");

	// if TOF is zero, then linkFrom_node is assumed to be the same as linkTo_node
	// To avoid having a segment with a TOF of zero, we delete one and update the
	// TOF and linkFrom_ID
	if(tof == 0){
		tof = linkFrom_seg.getTOF();						// Update tof

		// Get the next node down the line
		int new_linkFrom_ID = linkFrom_isOrigin ? linkFrom_seg.getTerminus() : linkFrom_seg.getOrigin();

		// Delete the end node and the segment that connects to it
		set->deleteSeg(linkFrom_seg.getID());
		set->deleteNode(linkFrom_ID);
		
		// Update objects and variables that depend on linkFrom_ID
		linkFrom_ID = new_linkFrom_ID;
		linkFrom_node = set->getNode(linkFrom_ID);
		int new_linkFrom_segIx = linkFrom_node.getLink(0) == tpat_linkable::INVALID_ID ? linkFrom_node.getLink(1) : linkFrom_node.getLink(0);
		
		if(new_linkFrom_segIx != tpat_linkable::INVALID_ID){
			linkFrom_seg = set->getSeg(new_linkFrom_segIx);
			linkFrom_isOrigin = linkFrom_seg.getOrigin() == linkFrom_node.getID();
		}else{
			// No segments left, just a node
			// Leave linkFrom_seg the same; this is used later to determine the direction of time
			// make linkFrom_isOrigin = false; no more segments originating from the final node
			linkFrom_isOrigin = false;
		}
	}

	// A mapping vector: index is the old node ID, value is the new node ID
	// All new IDs are initialized to the default INVALID_ID value
	std::vector<int> map_oldID_to_newID(set->getNextNodeID(), tpat_linkable::INVALID_ID);

	// Add all nodes from set to this object and keep track of new IDs
	for(int n = 0; n < set->getNumNodes(); n++){
		tpat_node node = set->getNodeByIx(n);

		// Remove all links to segments; these will be added back when the segments are added to this new arcset object
		node.clearLinks();
		map_oldID_to_newID[node.getID()] = addNode(node);
	}

	// Add all segments from set to this object and update the link IDs
	// The act of adding the segment will update the links in the newly added nodes
	for(int s = 0; s < set->getNumSegs(); s++){
		tpat_segment seg = set->getSegByIx(s);
		
		// Remap the origin and terminus to the new IDs
		if(seg.getOrigin() != tpat_linkable::INVALID_ID)
			seg.setOrigin(map_oldID_to_newID[seg.getOrigin()]);

		if(seg.getTerminus() != tpat_linkable::INVALID_ID)
			seg.setTerminus(map_oldID_to_newID[seg.getTerminus()]);

		addSeg(seg);
	}

	// Add a new segment to link the nodes from [set] to [this object]
	int origin = tpat_linkable::INVALID_ID, terminus = tpat_linkable::INVALID_ID;
	if(!linkTo_isOrigin){
		origin = linkTo_ID;
		terminus = map_oldID_to_newID[linkFrom_ID];
	}else if(!linkFrom_isOrigin){
		origin = map_oldID_to_newID[linkFrom_ID];
		terminus = linkTo_ID;
	}else{
		// Both are origins; the only double-origin node posibility is one where
		// each segment that originates from the node has a different time direction
		if(linkTo_seg.getTOF() < 0){
			origin = tof > 0 ? linkTo_ID : map_oldID_to_newID[linkFrom_ID];
			terminus = tof > 0 ? map_oldID_to_newID[linkFrom_ID] : linkTo_ID;
		}else if(linkFrom_seg.getTOF() < 0){
			origin = tof > 0 ? map_oldID_to_newID[linkFrom_ID] : linkTo_ID;
			terminus = tof > 0 ? linkTo_ID : map_oldID_to_newID[linkFrom_ID];
		}
	}

	// print();
	return addSeg(tpat_segment(origin, terminus, tof));
}//====================================================

/**
 *  @brief Remove all constraints from this arcset object.
 *  @details Note that this does not affect any constraints placed on
 *  individual nodes or segments
 */
void tpat_base_arcset::clearArcConstraints(){ cons.clear(); }


/**
 *  @brief Remove constraints from this arcset object as well as
 *  all its node and segment children
 */
void tpat_base_arcset::clearAllConstraints(){
	for(size_t n = 0; n < nodes.size(); n++){ nodes[n].clearConstraints(); }
	for(size_t s = 0; s < segs.size(); s++){ segs[s].clearConstraints(); }
	cons.clear();
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
 *  @throws tpat_exception when:
 *  * The ID is out of bounds
 *  * The ID is in bounds but the associated node has been deleted (no longer exists)
 *  * Deleting the node will result in multiple segment interfaces; these must be 
 *  created more explicitely by the user and will not be created automatically by
 *  this function.
 */
void tpat_base_arcset::deleteNode(int id){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::deleteNode: id out of bounds");

	int nodeIx = nodeIDMap[id];

	if(nodeIx != tpat_linkable::INVALID_ID){
		
		// printf("Attempting to delete node (ID %d)\n", id);
		// Get the node we're deleting
		tpat_node theNode = nodes[nodeIx];

		// Get the indices of any segments this node is linked to
		std::vector<int> linkedSegIxs;
		for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
			if(theNode.getLink(i) != tpat_linkable::INVALID_ID){
				int segIx = segIDMap[theNode.getLink(i)];
				if(segIx != tpat_linkable::INVALID_ID){
					linkedSegIxs.push_back(segIx);
					// printf("  linked to segment (ID %d) at index %d\n", theNode.getLink(i), segIx);
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
				tpat_segment *termSeg = &(segs[linkedSegIxs[termSegIx]]);
				tpat_segment *origSeg = &(segs[linkedSegIxs[(termSegIx + 1) % 2]]);

				// printf("  > Segment (ID %d) terminates at node (ID %d)\n", termSeg.getID(), id);
				// printf("  > Segment (ID %d) originates at node (ID %d)\n", origSeg.getID(), id);

				// Just to check
				if(termSeg->getTOF()*origSeg->getTOF() < 0){
					throw tpat_exception("tpat_base_arcset::deleteNode: I made an incorrect assumption about origin/terminus TOF direction!");
				}

				// Create a new segment
				tpat_segment combo(termSeg->getOrigin(), origSeg->getTerminus(), termSeg->getTOF() + origSeg->getTOF());

				// Replace the two segments with the new combined one
				deleteSeg(termSeg->getID());
				// print();
				deleteSeg(origSeg->getID());
				// print();
				addSeg(combo);
			}else{
				// Both must originate at node and have opposite time directions
				if(segs[linkedSegIxs[0]].getOrigin() != id || segs[linkedSegIxs[1]].getOrigin() != id){
					throw tpat_exception("tpat_base_arcset::deleteNode: double origin - made incorrect assumption about origin-ness");
				}

				if(segs[linkedSegIxs[0]].getTOF()*segs[linkedSegIxs[1]].getTOF() > 0){
					throw tpat_exception("tpat_base_arcset::deleteNode: double origin - made incorrect assumption about TOF");
				}

				int revSegIx = segs[linkedSegIxs[0]].getTOF() < 0 ? 0 : 1;
				tpat_segment *revSeg = &(segs[linkedSegIxs[revSegIx]]);
				tpat_segment *forwardSeg = &(segs[linkedSegIxs[(revSegIx+1) % 2]]);

				/*	It is possible that, in this case, the segment that originates from this node and proceeds
				 * 	in reverse time does not terminate at a node, but links to a forward-propagated segment instead.
				 * 	If this is so, then the combination of the two segments becomes a reverse-time segment. However,
				 * 	if the segments both link to nodes, then the default behavior is to construct a new forward-time
				 * 	segment to replace the reverse and forward time segments that originated from this node.
				 */
				tpat_segment combo;
				if(revSeg->getTerminus() != tpat_linkable::INVALID_ID){
					combo = tpat_segment(revSeg->getTerminus(), forwardSeg->getTerminus(), std::abs(revSeg->getTOF()) + forwardSeg->getTOF());
				}else{
					if(forwardSeg->getTerminus() == tpat_linkable::INVALID_ID)
						throw tpat_exception("tpat_base_arcset::deleteNode: Cannot delete node as both segments terminate at other segments");
					combo = tpat_segment(forwardSeg->getTerminus(), revSeg->getTerminus(), revSeg->getTOF() - forwardSeg->getTOF());
				}

				// Replace the two segments with the new one
				deleteSeg(revSeg->getID());
				deleteSeg(forwardSeg->getID());
				addSeg(combo);
			}
		}else if(linkedSegIxs.size() == 1){
			// Must be either first or last point, cannot delete!
			throw tpat_exception("tpat_base_arcset::deleteNode: Node is only linked to one segment and must, therefore, be the first or last node and cannot be deleted");
		}else{ /* Not linked to anything, just delete it! */ }


		nodes.erase(nodes.begin() + nodeIx);					// Remove node from vector
		nodeIDMap[id] = tpat_linkable::INVALID_ID;				// This ID will never be used again, set to INVALID_ID

		// Update vector that maps ID values to storage indices
		for(size_t n = id+1; n < nodeIDMap.size(); n++){
			if(nodeIDMap[n] != tpat_linkable::INVALID_ID)
				nodeIDMap[n]--;										// Decrement all indices after the node we just removed
		}
	}else{
		printf("Cannot Delete: Node with ID %02d was not located\n", id);
	}
}//===========================================

/**
 *  @brief Delete a segment with the specified ID.
 *  @details Additionally, any nodes linked to the specified segment are updated
 *  so that their link arrays no longer include a relationship with the soon-to-be
 *  deceased segment.
 * 
 *  @param id The ID of the segment to delete. If the ID is out of range, an exception
 *  is thrown. If the ID is in range but doesn't represent an existing segment,
 *  no deletion is made
 * 	@throws tpat_exception if <tt>id</tt> is out of bounds
 */
void tpat_base_arcset::deleteSeg(int id){
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_base_arcset::deleteSeg: Invalid ID (out of bounds)");

	int segIx = segIDMap[id];

	if(segIx != tpat_linkable::INVALID_ID){
		// printf("Deleting segment (ID %d)\n", id);

		tpat_segment *seg = &(segs[segIx]);
		// printf("  Retrieved segment (ID %d)\n", seg.getID());
		for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
			if(seg->getLink(i) != tpat_linkable::INVALID_ID){
				int nodeIx = nodeIDMap[seg->getLink(i)];
				if(nodeIx != tpat_linkable::INVALID_ID){
					// printf("  * Trying to remove a link to segment (ID %d) from node (ID %d)\n", seg->getID(), nodes[nodeIx].getID());
					nodes[nodeIx].removeLink(id);
				}else{
					// printf("Unable to remove link to segment (ID %d) from node (ID %d)\n", seg->getLink(i), id);
				}
			}
		}

		segs.erase(segs.begin() + segIx);

		segIDMap[id] = tpat_linkable::INVALID_ID;
		for(size_t s = id+1; s < segIDMap.size(); s++){
			if(segIDMap[s] != tpat_linkable::INVALID_ID)
				segIDMap[s]--;
		}
	}else{
		// printf("Cannot Delete: Segment with ID %02d was not located\n", id);
	}
}//===========================================

/**
 *  @brief Retrieve the acceleration vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the acceleration vector
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getAccel(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getAccel: Node ID out of range");

	return nodes[nodeIDMap[id]].getAccel();
}//====================================================

/**
 *	@brief Retrieve an acceleration on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last acceleration)
 *	@return the acceleration associated with the specified index
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getAccelByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::getAccelByIx: node index out of bounds");

	return nodes[ix].getAccel();
}//====================================================

/**
 *  @brief Retrieve a vector containing all the constraints applied to this arcset object.
 *  @details This vector does not include constraints placed on individual nodes or segments.
 *  @return a vector containing all the constraints applied to this arcset object.
 */
std::vector<tpat_constraint> tpat_base_arcset::getArcConstraints() const { return cons; }

/**
 *  @brief Determine what order to place the nodes and segments of this object
 *  into to achieve a chronological progression in forward time.
 * 
 *  @param set the arcset object to order
 *  @return a vector of tpat_arc_piece objects that represent the chronological 
 *  order of the nodes and segments
 */
std::vector<tpat_arc_piece> tpat_base_arcset::getChronoOrder() const{
	
	std::vector<tpat_arc_piece> pieces;
	if(nodes.size() == 0){
		printErr("tpat_base_arcset::getChronoOrder: First node is invalid... exiting\n");
		return pieces;
	}

	tpat_node node0 = nodes[0];
	// printf("Beginning with node ID %d\n", node0.getID());

	pieces.push_back(tpat_arc_piece(tpat_arc_piece::NODE, node0.getID()));
	for(int dir = 1; dir > -2; dir -= 2){
		// printf("Direction is %s\n", dir > 0 ? "[FORWARD]" : "[BACKWARD]");

		bool go = true;
		tpat_node node = node0;
		
		while(go){
			// Get segments connected to node
			bool foundNextSeg = false;
			for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
				if(node.getLink(i) != tpat_linkable::INVALID_ID){

					tpat_segment seg = getSeg(node.getLink(i));
					// printf(" (checking out segment ID %d): ", seg.getID());

					if(dir > 0 && ((seg.getTerminus() == node.getID() && seg.getTOF() < 0) || (seg.getOrigin() == node.getID() && seg.getTOF() > 0))){
						// This segment moves forward in time from node
						pieces.push_back(tpat_arc_piece(tpat_arc_piece::SEG, seg.getID()));
						foundNextSeg = true;
						// printf("USE THIS\n");
					}else if(dir < 0 && ((seg.getTerminus() == node.getID() && seg.getTOF() > 0) || (seg.getOrigin() == node.getID() && seg.getTOF() < 0))){
						// This segment moves in reverse time from node
						pieces.insert(pieces.begin(), tpat_arc_piece(tpat_arc_piece::SEG, seg.getID()));
						foundNextSeg = true;
						// printf("USE THIS\n");
					}else{
						foundNextSeg = false;
						// printf("do not use\n");
					}

					if(foundNextSeg){
						// Get next node
						int nextNodeID = tpat_linkable::INVALID_ID;
						if(seg.getTerminus() == node.getID())
							nextNodeID = seg.getOrigin();
						else
							nextNodeID = seg.getTerminus();

						if(nextNodeID != tpat_linkable::INVALID_ID){
							// printf("Next node has ID %d\n", nextNodeID);
							node = getNode(nextNodeID);

							if(dir > 0)
								pieces.push_back(tpat_arc_piece(tpat_arc_piece::NODE, node.getID()));
							else
								pieces.insert(pieces.begin(), tpat_arc_piece(tpat_arc_piece::NODE, node.getID()));
						}else{
							// printf("No next node, switching direction...\n");
							// The segment terminates/originates without a node; end this string
							// TODO - This is where we would implement a segment-to-segment link
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
 *	@param ix the index of the coordinate: 0 = x, 1 = y, etc.
 *	@return a vector containing the specified coordinate for all
 *	nodes (not necessarily in chronological order)
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getCoord(int ix) const{
	if(ix >= 6)
		throw tpat_exception("tpat_base_arcset::getCoord: Index Out of Range");

	std::vector<double> coord;
	for(size_t n = 0; n < nodes.size(); n++){
		coord.push_back(nodes[n].getState()[ix]);
	}

	return coord;
}//====================================================

/**
 *  @brief Retrieve the epoch associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the epoch
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
double tpat_base_arcset::getEpoch(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getEpoch: Node ID out of range");

	return nodes[nodeIDMap[id]].getEpoch();
}//====================================================

/**
 *  @brief Retrieve the epoch of a specific node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @return The epoch associated with the specified node
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
double tpat_base_arcset::getEpochByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::getEpochByIx: node index out of bounds");

	return nodes[ix].getEpoch();
}//=====================================================

/**
 *	@brief Retrieve a set of extra parameters for the specified node
 *	@param n the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@param ix the index of the extra parameter
 *	@return a vector containing the extra parameter at the specified step and index
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getExtraParam(int n, int ix) const{
	if(n < 0)
		n += nodes.size();

	if(ix < 0 || ix >= (int)(extraParamRowSize.size()))
		throw tpat_exception("tpat_base_arcset::getExtraParam: parameter index out of bounds");

	int startIx = 0;
	for(int i = 0; i < ix; i++)
		startIx += extraParamRowSize[i];

	int size = extraParamRowSize[ix];
	std::vector<double> extraParam = nodes[n].getExtraParams();
	return std::vector<double>(extraParam.begin()+startIx, extraParam.begin()+startIx + size);
}//====================================================

/**
 *  @brief Retrieve the value of the ID that will be assigned to the next 
 *  node added to this arcset object
 *  @return the next node ID
 */
int tpat_base_arcset::getNextNodeID() const { return nextNodeID; }

/**
 *  @brief Retrieve the value of the ID that will be assigned to the next 
 *  segment added to this arcset object
 *  @return the next segment ID
 */
int tpat_base_arcset::getNextSegID() const { return nextSegID; }

/**
 *	@brief Retrieve the number of nodes
 *	@return the number of nodes
 */
int tpat_base_arcset::getNumNodes() const { return (int)(nodes.size()); }

/**
 *	@brief Retrieve the number of segments
 *	@return the number of segments
 */
int tpat_base_arcset::getNumSegs() const { return (int)(segs.size()); }

/**
 *  @brief Retrieve the total number of constraints contained by all nodes, segments,
 *  and the arcset object itself
 *  @return the total number of constraints applied to this object and its children
 */
int tpat_base_arcset::getNumCons() const {
	int conCount = (int)cons.size();
	for(size_t n = 0; n < nodes.size(); n++){ conCount += nodes[n].getNumCons(); }
	for(size_t s = 0; s < segs.size(); s++){ conCount += segs[s].getNumCons(); }

	return conCount;
}//===================================================

/**
 *  @brief Retrieve a specific node
 * 
 *  @param id the ID of the desired node.
 *	
 *  @return the node located with the specified ID
 *  @throws tpat_exception if <tt>id</tt> is out of bounds or if no node exists with the specified ID
 */
tpat_node tpat_base_arcset::getNode(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getNode: Invalid ID (out of bounds)");

	int ix = nodeIDMap[id];
	if(ix != tpat_linkable::INVALID_ID && ix < (int)(nodes.size()) && ix >= 0){
		return nodes[ix];
	}else{
		throw tpat_exception("tpat_base_arcset::getNode: Could not locate a node with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a node based on its index in the storage array
 * 
 *  @param ix The index of the node; if <tt>ix</tt> is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a segment at the specified index
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
tpat_node tpat_base_arcset::getNodeByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::getNodeByIx: Index out of bounds");

	return nodes[ix];
}//=====================================================

/**
 *  @brief Retrieve a specific node
 * 
 *  @param id the ID of the desired node
 *	
 *  @return the node located with the specified ID
 *  @throws tpat_exception if <tt>id</tt> is out of bounds or if no segment exists with the specified ID
 */
tpat_segment tpat_base_arcset::getSeg(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getSeg: Invalid ID (out of bounds)");

	int ix = segIDMap[id];
	if(ix != tpat_linkable::INVALID_ID && ix < (int)(segs.size()) && ix >= 0){
		return segs[ix];
	}else{
		throw tpat_exception("tpat_base_arcset::getSeg: Could not locate a segment with the specified ID");
	}
}//====================================================

/**
 *  @brief Retrieve a segment based on its index in the storage array
 * 
 *  @param ix The index of the segment; if <tt>ix</tt> is negative, the index will
 *  count backwards from the end of the storage array.
 *  @return a segment at the specified index
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
tpat_segment tpat_base_arcset::getSegByIx(int ix) const{
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_base_arcset::getSegByIx: Index out of bounds");

	return segs[ix];
}//=====================================================

/**
 *  @brief Retrieve the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the state vector
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getState(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getState: Node ID out of range");

	return nodes[nodeIDMap[id]].getState();
}//====================================================

/**
 *	@brief Retrieve a position-velocity state on the arc
 *	@param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@return the state associated with the specified index
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> tpat_base_arcset::getStateByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::getStateByIx: node index out of bounds");

	return nodes[ix].getState();
}//====================================================

/**
 *  @brief Retrieve the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the STM
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
MatrixXRd tpat_base_arcset::getSTM(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getSTM: Node ID out of range");

	return segs[segIDMap[id]].getSTM();
}//====================================================

/**
 *	@brief Retrieve an STM on the arc
 *	@param ix the segment index. If it is negative, the index will count backwards
 *	from the end of the <tt>segs</tt> storage array
 *	
 *	@return the STM associated with the specified index
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
MatrixXRd tpat_base_arcset::getSTMByIx(int ix) const{
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_base_arcset::getSTMByIx: segment index out of bounds");

	return segs[ix].getSTM();
}//====================================================

/**
 *	@brief Retrieve the a pointer to the system data object associated with this arc
 *	@return a pointer to the system data object associated with this arc
 */
const tpat_sys_data* tpat_base_arcset::getSysData() const { return sysData; }

/**
 *  @brief Retrieve the time-of-flight associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the time-of-flight
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
double tpat_base_arcset::getTOF(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_base_arcset::getTOF: Segment ID out of range");

	return segs[segIDMap[id]].getTOF();
}//====================================================
/**
 *	@brief Get the time-of-flight for a specific segment
 *	@param ix node index (NOT the ID); if less than 0, the index counts
 *	backwards from the end of the nodeset
 *	@return non-dimensional time-of-flight along the specified segment
 *	@throws tpat_exception if <tt>ix</tt> is out of bounds
 */
double tpat_base_arcset::getTOFByIx(int ix) const {
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= ((int)segs.size()))
		throw tpat_exception("tpat_nodeset::getTOFByIx: invalid segment index");

	return segs[ix].getTOF();
}//====================================================

/**
 *	@brief Retrieve the tolerance with which data in this object was computed
 *	@return the tolerance with which data in this object was computed
 */
double tpat_base_arcset::getTol() const { return tol; }

/**
 *  @brief Determine the total time-of-flight along this arc.
 *  @details This function sums the TOF along each segment; derived
 *  classes may override this function to use different methods.
 *  @return the total time-of-flight along this arc, units consistent
 *  with the tpat_sys_data object
 */
double tpat_base_arcset::getTotalTOF() const{
	double total = 0;
	for(size_t s = 0; s < segs.size(); s++){
		total += std::abs(segs[s].getTOF());
	}
	return total;
}//=================================================

/**
 *  @brief Rearrange the nodes and segments so that they are listed
 *  in chronological order in their storage arrays.
 *  
 *  @details This does not change the ID of any of the nodes or segments,
 *  only their index within the storage array. After calling this function,
 *  accessing the -1 node or segment will return the latest (in time) object.
 *  
 *  @throws tpat_exception if the getChronoOrder() sorting algorithm returns a 
 *  set of tpat_arc_piece objects that has a different size than the combined 
 *  node and segment vectors, the function is aborted as it is likely a node or 
 *  segment was skipped and we don't want to lose information.
 */
void tpat_base_arcset::putInChronoOrder(){
	// First, determine the chronological order
	// print();
	std::vector<tpat_arc_piece> pieces = getChronoOrder();

	if(pieces.size() != nodes.size() + segs.size()){
		printErr("Pieces has %zu elements, but there are %zu nodes and %zu segs\n", pieces.size(), nodes.size(), segs.size());
		throw tpat_exception("tpat_base_arcset::putInChronoOrder: The sorted vector does not include all nodes and segments; aborting to avoid losing data\n");
	}

	std::vector<tpat_node> newNodes;
	std::vector<tpat_segment> newSegs;
	std::vector<int> newNodeIDMap(nodeIDMap.size(), tpat_linkable::INVALID_ID);
	std::vector<int> newSegIDMap(segIDMap.size(), tpat_linkable::INVALID_ID);

	for(size_t i = 0; i < pieces.size(); i++){
		if(pieces[i].type == tpat_arc_piece::NODE){
			tpat_node node = getNode(pieces[i].ID);
			newNodeIDMap[node.getID()] = newNodes.size();
			newNodes.push_back(node);
		}else if(pieces[i].type == tpat_arc_piece::SEG){
			tpat_segment seg = getSeg(pieces[i].ID);
			newSegIDMap[seg.getID()] = newSegs.size();
			newSegs.push_back(seg);
		}
	}

	nodes = newNodes;
	segs = newSegs;
	nodeIDMap = newNodeIDMap;
	segIDMap = newSegIDMap;
}//=============================================

/**
 *  @brief Set the acceleration vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param accel the acceleration vector
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
void tpat_base_arcset::setAccel(int id, std::vector<double> accel){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::setAccel: Node ID out of range");

	nodes[nodeIDMap[id]].setAccel(accel);
}//====================================================

/**
 *  @brief Set the acceleration vector for a specific step/node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param accelVec 3-element (at least) vector of non-dimensional acceleration 
 *  values (ax, ay, az, ...); only the first three are used
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
void tpat_base_arcset::setAccelByIx(int ix, std::vector<double> accelVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::setAccelByIx: node index out of bounds");

	nodes[ix].setAccel(accelVec);
}//=================================================

/**
 *  @brief Set the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param state the state vector
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
void tpat_base_arcset::setState(int id, std::vector<double> state){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_base_arcset::setState: Node ID out of range");

	nodes[nodeIDMap[id]].setState(state);
}//====================================================

/**
 *  @brief Set the state vector for a specific step/node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arcset object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param stateVec 6-element (at least) vector of non-dimensional state 
 *  values (x, y, z, vx, vy, vz, ...); only the first six are used
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
void tpat_base_arcset::setStateByIx(int ix, std::vector<double> stateVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_base_arcset::setStateByIx: node index out of bounds");

	nodes[ix].setState(stateVec);
}//=================================================

/**
 *  @brief Set the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @param stm the STM
 *  @throws tpat_exception if <tt>id</tt> is out of bounds
 */
void tpat_base_arcset::setSTM(int id, MatrixXRd stm){
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_base_arcset::setSTM: Node ID out of range");

	segs[segIDMap[id]].setSTM(stm);
}//====================================================

/**
 *  @brief Set the STM for a specific step/node
 * 
 *  @param ix index of the segment with the <tt>segs</tt> storage array; if it is negative,
 *  it will count backwards from the end of the array.
 *  
 *  @param stm a 6x6 matrix containing the STM
 *  @throws tpat_exception if <tt>ix</tt> is out of bounds
 */
void tpat_base_arcset::setSTMByIx(int ix, MatrixXRd stm){
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_base_arcset::setSTMByIx: node index out of bounds");

	segs[ix].setSTM(stm);
}//=================================================

/**
 *	@brief Set the computational tolerance for this data object
 *	@param d the tolerance
 */
void tpat_base_arcset::setTol(double d){ tol = d; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy all data from the input arc data to this one
 *	@param d an arc data object reference
 */
void tpat_base_arcset::copyMe(const tpat_base_arcset &d){
	nodes = d.nodes;
	nodeIDMap = d.nodeIDMap;
	segs = d.segs;
	segIDMap = d.segIDMap;
	cons = d.cons;
	// sysData = d.sysData; // Copying ADDRESS of sys_data object
	numExtraParam = d.numExtraParam;
	extraParamRowSize = d.extraParamRowSize;
	tol = d.tol;
	nextNodeID = d.nextNodeID;
	nextSegID = d.nextSegID;
}//====================================================

/**
 *  @brief Initialize the vectors of node and segment objects from a *.mat file
 *  @details THIS FUNCTION MUST BE THE FIRST READ_DATA-TYPE FUNCTION CALLED because
 *  it clears the vectors and then initializes them by calculating the number
 *  of steps in the arcset object from the state vector. Individual nodes and segments are
 *  able to be called by index after this, though they will not contain data
 *  until another function is called to populate the data fields with values from 
 *  the *.mat file
 * 
 *  @param matFile pointer to an open matlab data file
 *  @param varName the name of a variable that has as many rows as there are
 *  steps along the data object. Valid variables typically include the time vector,
 *  state matrix, or acceleration matrix
 *  @throws tpat_exception if the state vector variable cannot be read from the data file
 */
void tpat_base_arcset::initNodesSegsFromMat(mat_t *matFile, const char* varName){
	matvar_t *stateMat = Mat_VarRead(matFile, varName);
	if(stateMat == NULL){
		throw tpat_exception("tpat_base_arcset::initNodeSegsFromMat: Could not read state data vector");
	}else{
		int numSteps = stateMat->dims[0];
		nodes.clear();
		segs.clear();
		nodeIDMap.clear();
		segIDMap.clear();
		tpat_node blank_node;
		tpat_segment blank_seg;
		
		// Create a set of nodes and segments all linked together in linear time
		for(int i = 0; i < numSteps; i++){
			nodes.push_back(blank_node);
			nodeIDMap.push_back(i);
			
			if(i > 0){
				segs.push_back(tpat_segment(i-1, i, NAN));
				segIDMap.push_back(i-1);
			}
		}
		// nodes.assign(numSteps, blank_node);	// Initialize array with a bunch of default objects
		// segs.assign(numSteps-1, blank_seg);
	}
	Mat_VarFree(stateMat);
}//======================================================

/**
 *  @brief Print a ASCII graphic of the arcset in chronological order
 */
void tpat_base_arcset::printInChrono() const{
	std::vector<tpat_arc_piece> pieces = getChronoOrder();

	for(size_t i = 0; i < pieces.size(); i++){
		tpat_arc_piece p = pieces[i];

		if(p.type == tpat_arc_piece::NODE){
			printf("[%02d]", p.ID);
		}else if(p.type == tpat_arc_piece::SEG){
			if(getSeg(p.ID).getTOF() > 0)
				printf("--(%02d)->", p.ID);
			else
				printf(">-(%02d)--", p.ID);
		}
	}
	printf("\n");
}//====================================================

/**
 *  @brief Read the state vector for this arcset object from a matlab data file
 *  @details This function must be called after initNodeSegsFromMat() as it
 *  populates the step vector objects with state data
 * 
 *  @param matFile pointer to an open matlab data file
 *  @param varName the name of the state variable (e.g., "State" or "Nodes")
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readStateFromMat(mat_t *matFile, const char* varName){
	matvar_t *stateMat = Mat_VarRead(matFile, varName);
	if(stateMat == NULL){
		throw tpat_exception("tpat_base_arcset::readStateFromMat: Could not read state data vector");
	}else{
		int numSteps = stateMat->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_base_arcset::readStateFromMat: Step vector has not been initialized!");
		}

		if(numSteps != (int)nodes.size()){
			throw tpat_exception("tpat_base_arcset::readStateFromMat: State vector has a different size than the initialized step vector");
		}

		if(stateMat->dims[1] != 6){
			throw tpat_exception("tpat_base_arcset::readStateFromMat: Incompatible data file: State width is not 6.");
		}

		if(stateMat->class_type == MAT_C_DOUBLE && stateMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(stateMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					double state[] = {0,0,0,0,0,0};
					state[0] = data[0*numSteps + i];
					state[1] = data[1*numSteps + i];
					state[2] = data[2*numSteps + i];
					state[3] = data[3*numSteps + i];
					state[4] = data[4*numSteps + i];
					state[5] = data[5*numSteps + i];

					nodes[i].setState(state);
				}
			}
		}else{
			throw tpat_exception("tpat_base_arcset::readStateFromMat: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(stateMat);
}//===============================================

/**
 *  @brief Read acceleration values from a matlab file
 * 
 *  @param matFile pointer to an open Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readAccelFromMat(mat_t *matFile){
	matvar_t *accelMat = Mat_VarRead(matFile, "Accel");
	if(accelMat == NULL){
		throw tpat_exception("tpat_base_arcset::readAccelFromMat: Could not read data vector");
	}else{
		int numSteps = accelMat->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_base_arcset::readAccelFromMat: Node vector has not been initialized!");
		}

		if(numSteps != (int)nodes.size()){
			throw tpat_exception("tpat_base_arcset::readAccelFromMat: Accel vector has a different size than the initialized node vector");
		}

		if(accelMat->dims[1] != 3){
			throw tpat_exception("tpat_base_arcset::readAccelFromMat: Incompatible data file: Accel width is not 3.");
		}

		if(accelMat->class_type == MAT_C_DOUBLE && accelMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(accelMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					double accel[] = {0,0,0};
					accel[0] = data[0*numSteps + i];
					accel[1] = data[1*numSteps + i];
					accel[2] = data[2*numSteps + i];

					nodes[i].setAccel(accel);
				}
			}
		}else{
			throw tpat_exception("tpat_base_arcset::readAccelFromMat: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(accelMat);
}//===============================================

/**
 *  @brief Read epoch times from a matlab file in a variable with the specified name
 * 
 *  @param matFile pointer to an open Matlab file
 *  @param varName The name of the variable within the Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readEpochFromMat(mat_t *matFile, const char* varName){
	matvar_t *epochMat = Mat_VarRead(matFile, varName);
	if(epochMat == NULL){
		throw tpat_exception("tpat_base_arcset::readEpochFromMat: Could not read data vector");
	}else{
		int numSteps = epochMat->dims[0];

		if(nodes.size() == 0)
			throw tpat_exception("tpat_base_arcset::readEpochFromMat: Node vector has not been initialized");

		if(numSteps != (int)nodes.size())
			throw tpat_exception("tpat_base_arcset::readEpochFromMat: Epoch vector has different size than the initialized node evctor");

		if(epochMat->dims[1] != 1)
			throw tpat_exception("tpat_base_arcset::readEpochFromMat: Incompatible data file: Epoch vector has more than one column");

		if(epochMat->class_type == MAT_C_DOUBLE && epochMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(epochMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					nodes[i].setEpoch(data[i]);
				}
			}
		}else{
			throw tpat_exception("tpat_base_arcset::readEpochFromMat: Incompatible data file: unsupported data type or class");
		}
	}
}//================================================

/**
 *  @brief Read State Transition Matrices from a matlab file
 * 
 *  @param matFile pointer to an open Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readSTMFromMat(mat_t *matFile){
	matvar_t *allSTM = Mat_VarRead(matFile, "STM");
	if(allSTM == NULL){
		throw tpat_exception("tpat_base_arcset::readSTMFromMat: Could not read data vector");
	}else{
		int numSteps = allSTM->dims[2];

		if(segs.size() == 0){
			throw tpat_exception("tpat_base_arcset::readSTMFromMat: Step vector has not been initialized!");
		}

		if(numSteps != (int)segs.size()){
			throw tpat_exception("tpat_base_arcset::readSTMFromMat: STM vector has a different size than the initialized step vector");
		}

		if(allSTM->dims[0] != 6 || allSTM->dims[1] != 6){
			throw tpat_exception("tpat_base_arcset::readSTMFromMat: Incompatible data file: STM is not 6x6.");
		}

		if(allSTM->class_type == MAT_C_DOUBLE && allSTM->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(allSTM->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					double stmEl[36];
					for(int j = 0; j < 36; j++){
						stmEl[j] = data[36*i + j];
					}

					MatrixXRd P = Eigen::Map<MatrixXRd>(stmEl, 6, 6);
					segs[i].setSTM(P.transpose());
				}
			}
		}else{
			throw tpat_exception("ttpat_arc_data::readSTMFromMat: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(allSTM);
}//===============================================

/**
 *  @brief Read times-of-flight from a matlab file in a variable with the specified name
 * 
 *  @param matFile pointer to an open Matlab file
 *  @param varName The name of the variable within the Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readTOFFromMat(mat_t *matFile, const char* varName){
	matvar_t *tofMat = Mat_VarRead(matFile, varName);
	if(tofMat == NULL){
		throw tpat_exception("tpat_base_arcset::readTOFFromMat: Could not read data vector");
	}else{
		int numSteps = tofMat->dims[0];

		if(segs.size() == 0)
			throw tpat_exception("tpat_base_arcset::readTOFFromMat: Node vector has not been initialized");

		if(numSteps != (int)segs.size())
			throw tpat_exception("tpat_base_arcset::readTOFFromMat: Epoch vector has different size than the initialized segment evctor");

		if(tofMat->dims[1] != 1)
			throw tpat_exception("tpat_base_arcset::readTOFFromMat: Incompatible data file: Epoch vector has more than one column");

		if(tofMat->class_type == MAT_C_DOUBLE && tofMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(tofMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					segs[i].setTOF(data[i]);
				}
			}
		}else{
			throw tpat_exception("tpat_base_arcset::readTOFFromMat: Incompatible data file: unsupported data type or class");
		}
	}
}//================================================

/**
 *  @brief Read values of the specified extra paramter from a matlab file
 * 
 *  @param matFile pointer to an open Matlab file
 *  @param varIx the index of the extra parameter variable within the <tt>extraParamRowSize</tt> array
 *  @param varName the name of the storage variable within the Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_base_arcset::readExtraParamFromMat(mat_t *matFile, int varIx, const char *varName){
	if(varIx >= numExtraParam || varIx < 0)
		throw tpat_exception("tpat_base_arcset::readExtraParamFromMat: Could not read extra parameter; index out of bounds");

	// Get starting index of this extra param within a arc step's extra parameter vector
	int ix0 = 0;
	for(int i = 0; i < varIx; i++){ ix0 += extraParamRowSize[i]; }

	matvar_t *matvar = Mat_VarRead(matFile, varName);
	if(matvar == NULL){
		throw tpat_exception("tpat_base_arcset::readExtraParamFromMat: Could not read data vector");
	}else{
		int numSteps = matvar->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_base_arcset::readExtraParamFromMat: Step vector has not been initialized!");
		}

		if(matvar->dims[1] != ((size_t)extraParamRowSize[varIx])){
			char message[64];
			sprintf(message, "tpat_base_arcset::readExtraParamFromMat: Incompatible data file: %s width is not %d", varName, extraParamRowSize[varIx]);
			throw tpat_exception(message);
		}

		if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(matvar->data);
			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					for(int c = 0; c < extraParamRowSize[varIx]; c++){
						nodes[i].setExtraParam(ix0+c, data[c*numSteps + i]);
					}
				}
			}
		}else{
			throw tpat_exception("tpat_base_arcset::readExtraParamFromMat: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//===============================================

/**
 *	@brief Save the acceleration vector to file
 *	@param matFile a pointer to the destination mat-file
 */
void tpat_base_arcset::saveAccel(mat_t *matFile) const{
	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> accel_colMaj(3*nodes.size());

	for(size_t r = 0; r < nodes.size(); r++){
		std::vector<double> accel = nodes[r].getAccel();
		for(int c = 0; c < 3; c++){
			accel_colMaj[c*nodes.size() + r] = accel[c];
		}
	}
	
	size_t dims[2] = {nodes.size(), 3};
	matvar_t *matvar = Mat_VarCreate("Accel", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(accel_colMaj[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Accel", MAT_COMPRESSION_NONE);
}//=====================================================

/**
 *	@brief Save all node epochs to file
 *	@param matFile a pointer to the destination mat-file
 */
void tpat_base_arcset::saveEpoch(mat_t *matFile) const{
	saveEpoch(matFile, "Epoch");
}//=====================================================

/**
 *	@brief Save all node epochs to file with a specified variable name
 *	@param matFile a pointer to the destination mat-file
 *	@param varName the name of the variable
 */
void tpat_base_arcset::saveEpoch(mat_t *matFile, const char* varName) const{
	std::vector<double> allEpochs(nodes.size());

	for(size_t n = 0; n < nodes.size(); n++){
		allEpochs[n] = nodes[n].getEpoch();
	}
	
	size_t dims[2] = {allEpochs.size(), 1};
	matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allEpochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, varName, MAT_COMPRESSION_NONE);
}//=====================================================

/**
 *	@brief Save one of the extra parameters to file
 *	@param matFile a pointer to the destination mat-file
 *	@param varIx the index of the parameter
 *	@param name the name of the variable being saved
 *	@throws tpat_exception if <tt>varIx</tt> is out of bounds
 */
void tpat_base_arcset::saveExtraParam(mat_t *matFile, int varIx, const char *name) const{
	if(varIx > numExtraParam || varIx < 0)
		throw tpat_exception("Could not save extra parameter; index out of bounds");

	// Get starting index of this extra param within a arc step's extra parameter vector
	int ix0 = 0;
	for(int i = 0; i < varIx; i++){ ix0 += extraParamRowSize[i]; }

	// Get the specified coordinate
	std::vector<double> param(extraParamRowSize[varIx]*nodes.size());
	for(size_t r = 0; r < nodes.size(); r++){
		std::vector<double> ep  = nodes[r].getExtraParams();
		for(int c = 0; c < extraParamRowSize[varIx];c++){
			// Save NAN (rather than un-allocated memeory) if the index is out of bounds
			if(ix0 + c < (int)(ep.size()))
				param[c*nodes.size() + r] = ep[ix0+c];
			else
				param[c*nodes.size() + r] = NAN;
		}
	}

	size_t dims[2] = {nodes.size(), static_cast<size_t>(extraParamRowSize[varIx])};
	matvar_t *matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(param[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, name, MAT_COMPRESSION_NONE);
}//======================================================

/**
 *	@brief Save the state vector [pos, vel] to a file with variable name "State"
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_base_arcset::saveState(mat_t *matFile) const{
	saveState(matFile, "State");
}//==================================================

/**
 *	@brief Save the state vector [pos, vel] to a file
 *	@param matFile a pointer to the destination matlab file 
 *	@param varName the name of the variable (e.g. "State" or "Nodes")
 */
void tpat_base_arcset::saveState(mat_t *matFile, const char* varName) const{
	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> posVel(6*nodes.size());

	for(size_t r = 0; r < nodes.size(); r++){
		std::vector<double> state = nodes[r].getState();
		for(int c = 0; c < 6; c++){
			posVel[c*nodes.size() + r] = state[c];
		}
	}

	// Next, create a matlab variable for the state and save it to the file
	/*	Create a matlab variable. Arguments are:
	 *	const char *name 	- varName, the name of the variable
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
	size_t dims[2] = {nodes.size(), 6};
	matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(posVel[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, varName, MAT_COMPRESSION_NONE);
}//======================================================

/**
 *	@brief Save the STMs to a file; STMs are stored in a 6x6xn array for 
 *	compatibility with existing MATLAB scripts
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_base_arcset::saveSTMs(mat_t *matFile) const{
	// Create one large vector to put all the STM elements in
	std::vector<double> allSTMEl(segs.size()*36);

	for (size_t n = 0; n < segs.size(); n++){
		// get the transpose of the STM matrix; we need to store it in column-major order
		// and it's currently in row-major order
		MatrixXRd P = segs[n].getSTM().transpose();
		// Retrieve the data from the matrix
		double *matData = P.data();
		// Store that data in our huge vector
		std::copy(matData, matData+36, &(allSTMEl[0]) + n*36);
	}

	size_t dims[3] = {6, 6, segs.size()};
	matvar_t *matvar = Mat_VarCreate("STM", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, &(allSTMEl[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "STM", MAT_COMPRESSION_NONE);
}//======================================================

/**
 *	@brief Save all segment times-of-flight to file with a specified variable name
 *	@param matFile a pointer to the destination mat-file
 *	@param varName the name of the variable
 */
void tpat_base_arcset::saveTOF(mat_t *matFile, const char* varName) const{
	std::vector<double> allTOFs(segs.size());

	for(size_t n = 0; n < segs.size(); n++){
		allTOFs[n] = segs[n].getTOF();
	}
	
	size_t dims[2] = {allTOFs.size(), 1};
	matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allTOFs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, varName, MAT_COMPRESSION_NONE);
}//=====================================================




