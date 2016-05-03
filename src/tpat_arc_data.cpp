/**
 *  @file tpat_arc_data.cpp
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

#include "tpat_arc_data.hpp"

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
tpat_arc_data::tpat_arc_data(const tpat_sys_data *sys) : sysData(sys){}

/**
 *	@brief Copy constructor
 *	@param d an arc_data reference
 */
tpat_arc_data::tpat_arc_data(const tpat_arc_data &d){
	copyMe(d);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_arc_data::~tpat_arc_data(){
	nodes.clear();
	segs.clear();
	nodeIDMap.clear();
	segIDMap.clear();
	cons.clear();
	extraParamRowSize.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Set this object equal to another
 *	@param d an arc_data reference
 *	@return a reference to this arc_data object
 */
tpat_arc_data& tpat_arc_data::operator =(const tpat_arc_data &d){
	copyMe(d);
	return *this;
}//====================================================

/**
 *	@brief Concatenate two arc_data objects
 *
 * 	When adding A + B, if the final state of A and initial state
 *	of B are the same, this algorithm will skip the initial state
 *	of B in the concatenation to avoid duplicating a node or state.
 *	
 *	The process begins by creating a new arc_data object and setting
 *	it equal to A. B is then appended to the new object and the STM's
 *	from B are updated to be continuous with those found in A. Note 
 *	that if A and B are not conitnuous, the new STMs will be incorrect
 *	and will have no practical meaning
 *
 *	@param rhs the right-hand-side of the addition operation
 *	@return a reference to the concatenated arc_data object
 */
tpat_arc_data& tpat_arc_data::operator +=(const tpat_arc_data &rhs){
	if( *sysData != *(rhs.sysData))
		throw tpat_exception("tpat_arc_data::+=: Cannot concatenate data sets from different systems");

	throw tpat_exception("tpat_arc_data::+=: Not implemented");
	// if(steps.size() == 0){
	// 	copyMe(rhs);
	// 	return *this;
	// }

	// if(rhs.steps.size() == 0)
	// 	return *this;

	// int skipShift = 0;
	// double newTol = tol > rhs.tol ? tol : rhs.tol;
	// if(newTol == 0)
	// 	tol = 1e-9;

	// if(steps[0] == rhs.steps[0])
	// 	skipShift = 1;

	// // Delete data from the lhs if anything is duplicated; chose LHS because it will serve nodesets better
	// if(skipShift > 0)
	// 	steps.erase(steps.end()-skipShift, steps.end());

	// size_t lhs_numSteps = steps.size();

	// // Copy data from rhs
	// steps.insert(steps.end(), rhs.steps.begin(), rhs.steps.end());

	// // Adjust STMs (Assuming arcs are continuous)
	// MatrixXRd lhs_lastSTM = steps[lhs_numSteps-1].getSTM(); 	// PHI(t1, t0)
	// for(size_t n = lhs_numSteps; n < steps.size(); n++){
	// 	MatrixXRd oldSTM = steps[n].getSTM();	// PHI(t2, t1)

	// 	// Multiply PHI(t2, t1)*PHI(t1, t0) to get PHI(t2, t0)
	// 	steps[n].setSTM(oldSTM*lhs_lastSTM);

	// 	// Update constraint step numbers
	// 	steps[n].setConstraintNodeNum(n);
	// }

	// return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Add a constraint to the arc_data object.
 *  @details The constraint application type determines whether the constraint
 *  is applied to a node, segment, or the arc as a whole.
 * 
 *  @param con the constraint
 *  @see tpat_constraint::tpat_conApp_tp
 */
void tpat_arc_data::addConstraint(tpat_constraint con){
	int id = con.getNode();
	switch(con.getAppType()){
		case tpat_constraint::APP_TO_NODE:
			if(id < 0 || id >= (int)(nodeIDMap.size()))
				throw tpat_exception("tpat_arc_data::addConstraint: Node ID out of bounds");

			nodes[nodeIDMap[id]].addConstraint(con);
			break;
		case tpat_constraint::APP_TO_SEG:
			if(id < 0 || id >= (int)(segIDMap.size()))
				throw tpat_exception("tpat_arc_data::addConstraint: Segment ID out of bounds");

			segs[segIDMap[id]].addConstraint(con);
			break;
		case tpat_constraint::APP_TO_ARC:
			cons.push_back(con);
			break;
		default:
			throw tpat_exception("tpat_arc_data::addConstraint: Constraint application type is unknown");
	}
}//=====================================================

/**
 *  @brief Add a node to this data object
 *  @details A unique key is assigned to the node when it is added
 * 
 *  @param n the node to add
 *  @return the ID assigned to the node
 */
int tpat_arc_data::addNode(tpat_node n){
	n.clearLinks();			// Cannot have any links coming in
	n.setID(nextNodeID);
	nodes.push_back(n);
	nodeIDMap.push_back(nodes.size()-1);

	return nextNodeID++;
}//=====================================================


int tpat_arc_data::addSeg(tpat_segment s){
	s.setID(nextSegID);
	
	bool foundValidLink = false;
	for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
		// Get the ID of one of the nodes this segment is linked to
		int linkedNodeID = s.getLink(i);
		if(linkedNodeID != tpat_linkable::INVALID_ID){

			// Get the index of that node within the storage array
			int linkedNodeIx = nodeIDMap[linkedNodeID];
			if(linkedNodeIx != tpat_linkable::INVALID_ID){

				// See if the node is linked to any other segments
				tpat_node theNode = nodes[linkedNodeIx];
				int secondaryLinks = 0;
				for(int j = 0; j < tpat_linkable::NUM_LINKS; j++){
					if(theNode.getLink(j) != tpat_linkable::INVALID_ID){
						int nearSegIx = segIDMap[theNode.getLink(j)];
						if(nearSegIx != tpat_linkable::INVALID_ID){
							secondaryLinks++;
							// If the node is linked to another segment, get that segment and compare it to this one
							tpat_segment nearSeg = segs[nearSegIx];
							bool sameLinkType = nearSeg.getLink(i) == linkedNodeID;
							bool sameTimeDir = nearSeg.getTOF()*s.getTOF() > 0;

							if(sameLinkType && sameTimeDir){
								// either a time collision (both terminate) or parallel structure (both originate)
								throw tpat_exception("tpat_arc_data::addSeg: either time collision or parallel structure");
							}else if(!sameLinkType && !sameTimeDir){
								// parallel structure
								throw tpat_exception("tpat_arc_data::addSeg: parallel structure");
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
					nodes[linkedNodeIx].addLink(nextSegID);		// OK, looks good from here; node will check for duplicate linkage
					// printf("Node (ID %d) is now linked to segment (ID %d)\n", theNode.getID(), nextSegID);
				}
			}else{
				// Linked node has valid ID but isn't part of this arcset
				throw tpat_exception("tpat_arc_data::addSeg: Linked node is not part of this arc_data object");
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
 *  @brief Remove all constraints from this arc_data object.
 *  @details Note that this does not affect any constraints placed on
 *  individual nodes or segments
 */
void tpat_arc_data::clearArcConstraints(){ cons.clear(); }


/**
 *  @brief Remove constraints from this arc_data object as well as
 *  all its node and segment children
 */
void tpat_arc_data::clearAllConstraints(){
	for(size_t n = 0; n < nodes.size(); n++){ nodes[n].clearConstraints(); }
	for(size_t s = 0; s < segs.size(); s++){ segs[s].clearConstraints(); }
	cons.clear();
}//====================================================

/**
 *  @brief Delete the node with the specified ID
 *  @details In addition to deleting the desired node, the arc_data object is "healed"
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
 *  @throws tpat_exception when the ID is out of bounds, or if the conditions discussed above
 *  are not met
 */
void tpat_arc_data::deleteNode(int id){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::deleteNode: Invalid ID (out of bounds)");

	int nodeIx = nodeIDMap[id];

	if(nodeIx != tpat_linkable::INVALID_ID){
		
		// printf("Attempting to delete node (ID %d)\n", id);
		// Get the node we're deleting
		tpat_node theNode = nodes[nodeIx];

		// Get the indices of any segments this node is linked to
		std::vector<int> linkedSegIxs;
		for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
			int segIx = segIDMap[theNode.getLink(i)];
			if(segIx != tpat_linkable::INVALID_ID){
				linkedSegIxs.push_back(segIx);
				// printf("  linked to segment (ID %d) at index %d\n", theNode.getLink(i), segIx);
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
				tpat_segment termSeg = segs[linkedSegIxs[termSegIx]];
				tpat_segment origSeg = segs[linkedSegIxs[(termSegIx + 1) % 2]];

				// printf("  > Segment (ID %d) terminates at node (ID %d)\n", termSeg.getID(), id);
				// printf("  > Segment (ID %d) originates at node (ID %d)\n", origSeg.getID(), id);

				// Just to check
				if(termSeg.getTOF()*origSeg.getTOF() < 0){
					throw tpat_exception("tpat_arc_data::deleteNode: I made an incorrect assumption about origin/terminus TOF direction!");
				}

				// Create a new segment
				tpat_segment combo(termSeg.getOrigin(), origSeg.getTerminus(), termSeg.getTOF() + origSeg.getTOF());

				// Replace the two segments with the new combined one
				deleteSeg(termSeg.getID());
				// print();
				deleteSeg(origSeg.getID());
				// print();
				addSeg(combo);
			}else{
				// Both must originate at node and have opposite time directions
				if(segs[linkedSegIxs[0]].getOrigin() != id || segs[linkedSegIxs[1]].getOrigin() != id){
					throw tpat_exception("tpat_arc_data::deleteNode: double origin - made incorrect assumption about origin-ness");
				}

				if(segs[linkedSegIxs[0]].getTOF()*segs[linkedSegIxs[1]].getTOF() > 0){
					throw tpat_exception("tpat_arc_data::deleteNode: double origin - made incorrect assumption about TOF");
				}

				int revSegIx = segs[linkedSegIxs[0]].getTOF() < 0 ? 0 : 1;
				tpat_segment revSeg = segs[linkedSegIxs[revSegIx]];
				tpat_segment forwardSeg = segs[linkedSegIxs[(revSegIx+1) % 2]];

				/*	It is possible that, in this case, the segment that originates from this node and proceeds
				 * 	in reverse time does not terminate at a node, but links to a forward-propagated segment instead.
				 * 	If this is so, then the combination of the two segments becomes a reverse-time segment. However,
				 * 	if the segments both link to nodes, then the default behavior is to construct a new forward-time
				 * 	segment to replace the reverse and forward time segments that originated from this node.
				 */
				tpat_segment combo;
				if(revSeg.getTerminus() != tpat_linkable::INVALID_ID){
					combo = tpat_segment(revSeg.getTerminus(), forwardSeg.getTerminus(), std::abs(revSeg.getTOF()) + forwardSeg.getTOF());
				}else{
					if(forwardSeg.getTerminus() == tpat_linkable::INVALID_ID)
						throw tpat_exception("tpat_arc_data::deleteNode: Cannot delete node as both segments terminate at other segments");
					combo = tpat_segment(forwardSeg.getTerminus(), revSeg.getTerminus(), revSeg.getTOF() - forwardSeg.getTOF());
				}

				// Replace the two segments with the new one
				deleteSeg(revSeg.getID());
				deleteSeg(forwardSeg.getID());
				addSeg(combo);
			}
		}else if(linkedSegIxs.size() == 1){
			// Must be either first or last point, cannot delete!
			// printf("basic_arcset::deleteNode ERROR: Node (ID %d) is only linked to one segment, i.e., it is the first or last node and cannot be deleted\n", id);
			throw std::exception();
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
 *  @throw tpat_exception
 */
void tpat_arc_data::deleteSeg(int id){
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_arc_data::deleteSeg: Invalid ID (out of bounds)");

	int segIx = segIDMap[id];

	if(segIx != tpat_linkable::INVALID_ID){
		// printf("Deleting segment (ID %d)\n", id);

		tpat_segment seg = segs[segIx];
		// printf("  Retrieved segment (ID %d)\n", seg.getID());
		for(int i = 0; i < tpat_linkable::NUM_LINKS; i++){
			int nodeIx = nodeIDMap[seg.getLink(i)];
			if(nodeIx != tpat_linkable::INVALID_ID){
				// printf("  * Trying to remove a link to segment (ID %d) from node (ID %d)\n", seg.getID(), nodes[nodeIx].getID());
				nodes[nodeIx].removeLink(id);
			}else{
				// printf("Unable to remove link to segment (ID %d) from node (ID %d)\n", seg.getLink(i), id);
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
 */
std::vector<double> tpat_arc_data::getAccel(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::getAccel: Node ID out of range");

	return nodes[nodeIDMap[id]].getAccel();
}//====================================================

/**
 *	@brief Retrieve an acceleration on the arc
 *	@param ix the step index. If it is negative, the index will count backwards
 *	from the end of the arc (e.g. ix = -1 will return the last acceleration)
 *	@return the acceleration associated with the specified index
 */
std::vector<double> tpat_arc_data::getAccelByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::getAccelByIx: node index out of bounds");

	return nodes[ix].getAccel();
}//====================================================

/**
 *  @brief Retrieve a vector containing all the constraints applied to this arc_data object.
 *  @details This vector does not include constraints placed on individual nodes or segments.
 *  @return a vector containing all the constraints applied to this arc_data object.
 */
std::vector<tpat_constraint> tpat_arc_data::getArcConstraints() const { return cons; }

/**
 *	@brief Get a vector of one coordinate for all nodes
 *	@param ix the index of the coordinate: 0 = x, 1 = y, etc.
 *	@return a vector containing the specified coordinate for all
 *	nodes (not necessarily in chronological order)
 */
std::vector<double> tpat_arc_data::getCoord(int ix) const{
	if(ix >= 6)
		throw tpat_exception("tpat_arc_data::getCoord: Index Out of Range");

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
 */
double tpat_arc_data::getEpoch(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::getEpoch: Node ID out of range");

	return nodes[nodeIDMap[id]].getEpoch();
}//====================================================

/**
 *  @brief Retrieve the epoch of a specific node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arc_data object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @return The epoch associated with the specified node
 */
double tpat_arc_data::getEpochByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::getEpochByIx: node index out of bounds");

	return nodes[ix].getEpoch();
}//=====================================================

/**
 *	@brief Retrieve a set of extra parameters for the specified node
 *	@param n the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arc_data object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@param ix the index of the extra parameter
 *	@return a vector containing the extra parameter at the specified step and index
 */
std::vector<double> tpat_arc_data::getExtraParam(int n, int ix) const{
	if(n < 0)
		n += nodes.size();

	if(ix < 0 || ix >= (int)(extraParamRowSize.size()))
		throw tpat_exception("tpat_arc_data::getExtraParam: parameter index out of bounds");

	int startIx = 0;
	for(int i = 0; i < ix; i++)
		startIx += extraParamRowSize[i];

	int size = extraParamRowSize[ix];
	std::vector<double> extraParam = nodes[n].getExtraParams();
	return std::vector<double>(extraParam.begin()+startIx, extraParam.begin()+startIx + size);
}//====================================================

/**
 *	@brief Retrieve the number of nodes
 *	@return the number of nodes
 */
int tpat_arc_data::getNumNodes() const { return (int)(nodes.size()); }

/**
 *	@brief Retrieve the number of segments
 *	@return the number of segments
 */
int tpat_arc_data::getNumSegs() const { return (int)(segs.size()); }

/**
 *  @brief Retrieve the total number of constraints contained by all nodes, segments,
 *  and the arc_data object itself
 *  @return the total number of constraints applied to this object and its children
 */
int tpat_arc_data::getNumCons() const {
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
tpat_node tpat_arc_data::getNode(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::getNode: Invalid ID (out of bounds)");

	int ix = nodeIDMap[id];
	if(ix != tpat_linkable::INVALID_ID){
		return nodes[ix];
	}else{
		throw tpat_exception("tpat_arc_data::getNode: Could not locate a node with the specified ID");
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
tpat_node tpat_arc_data::getNodeByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::getNodeByIx: Index out of bounds");

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
tpat_segment tpat_arc_data::getSeg(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_arc_data::getSeg: Invalid ID (out of bounds)");

	int ix = segIDMap[id];
	if(ix != tpat_linkable::INVALID_ID){
		return segs[ix];
	}else{
		throw tpat_exception("tpat_arc_data::getSeg: Could not locate a segment with the specified ID");
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
tpat_segment tpat_arc_data::getSegByIx(int ix) const{
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_arc_data::getSegByIx: Index out of bounds");

	return segs[ix];
}//=====================================================

/**
 *  @brief Retrieve the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the state vector
 */
std::vector<double> tpat_arc_data::getState(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::getState: Node ID out of range");

	return nodes[nodeIDMap[id]].getState();
}//====================================================

/**
 *	@brief Retrieve a position-velocity state on the arc
 *	@param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arc_data object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *	@return the state associated with the specified index
 */
std::vector<double> tpat_arc_data::getStateByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::getStateByIx: node index out of bounds");

	return nodes[ix].getState();
}//====================================================

/**
 *  @brief Retrieve the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the STM
 */
MatrixXRd tpat_arc_data::getSTM(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_arc_data::getSTM: Node ID out of range");

	return segs[segIDMap[id]].getSTM();
}//====================================================

/**
 *	@brief Retrieve an STM on the arc
 *	@param ix the segment index. If it is negative, the index will count backwards
 *	from the end of the <tt>segs</tt> storage array
 *	
 *	@return the STM associated with the specified index
 */
MatrixXRd tpat_arc_data::getSTMByIx(int ix) const{
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_arc_data::getSTMByIx: segment index out of bounds");

	return segs[ix].getSTM();
}//====================================================

/**
 *	@brief Retrieve the a pointer to the system data object associated with this arc
 *	@return a pointer to the system data object associated with this arc
 */
const tpat_sys_data* tpat_arc_data::getSysData() const { return sysData; }

/**
 *  @brief Retrieve the time-of-flight associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @return the time-of-flight
 */
double tpat_arc_data::getTOF(int id) const{
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_arc_data::getTOF: Segment ID out of range");

	return segs[segIDMap[id]].getTOF();
}//====================================================
/**
 *	@brief Get the time-of-flight for a specific node
 *	@param ix node index (NOT the ID); if less than 0, the index counts
 *	backwards from the end of the nodeset
 *	@return non-dimensional time-of-flight between nodes ix and ix+1
 */
double tpat_arc_data::getTOFByIx(int ix) const {
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
double tpat_arc_data::getTol() const { return tol; }

/**
 *  @brief Determine the total time-of-flight along this arc.
 *  @details This function sums the TOF along each segment; derived
 *  classes may override this function to use different methods.
 *  @return the total time-of-flight along this arc, units consistent
 *  with the tpat_sys_data object
 */
double tpat_arc_data::getTotalTOF() const{
	double total = 0;
	for(size_t s = 0; s < segs.size(); s++){
		total += std::abs(segs[s].getTOF());
	}
	return total;
}//=================================================

/**
 *  @brief Set the acceleration vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param accel the acceleration vector
 */
void tpat_arc_data::setAccel(int id, std::vector<double> accel){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::setAccel: Node ID out of range");

	nodes[nodeIDMap[id]].setAccel(accel);
}//====================================================

/**
 *  @brief Set the acceleration vector for a specific step/node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arc_data object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param accelVec 3-element (at least) vector of non-dimensional acceleration 
 *  values (ax, ay, az, ...); only the first three are used
 */
void tpat_arc_data::setAccelByIx(int ix, std::vector<double> accelVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::setAccelByIx: node index out of bounds");

	nodes[ix].setAccel(accelVec);
}//=================================================

/**
 *  @brief Set the state vector associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param state the state vector
 */
void tpat_arc_data::setState(int id, std::vector<double> state){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_arc_data::setState: Node ID out of range");

	nodes[nodeIDMap[id]].setState(state);
}//====================================================

/**
 *  @brief Set the state vector for a specific step/node
 * 
 *  @param ix the node index within the <tt>nodes</tt> storage array; This value
 *	is not necessarily the same as the unique ID assigned to the node when it 
 *	was added to the arc_data object. If <tt>n</tt> is negative, this index will
 *	cound backwards from the end of the array.
 *	
 *  @param stateVec 6-element (at least) vector of non-dimensional state 
 *  values (x, y, z, vx, vy, vz, ...); only the first six are used
 */
void tpat_arc_data::setStateByIx(int ix, std::vector<double> stateVec){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix >= (int)(nodes.size()))
		throw tpat_exception("tpat_arc_data::setStateByIx: node index out of bounds");

	nodes[ix].setState(stateVec);
}//=================================================

/**
 *  @brief Set the STM associated with a segment
 *  with the specified ID
 * 
 *  @param id the ID of a segment
 *  @param stm the STM
 */
void tpat_arc_data::setSTM(int id, MatrixXRd stm){
	if(id < 0 || id >= (int)(segIDMap.size()))
		throw tpat_exception("tpat_arc_data::setSTM: Node ID out of range");

	segs[segIDMap[id]].setSTM(stm);
}//====================================================

/**
 *  @brief Set the STM for a specific step/node
 * 
 *  @param ix index of the segment with the <tt>segs</tt> storage array; if it is negative,
 *  it will count backwards from the end of the array.
 *  
 *  @param stm a 6x6 matrix containing the STM
 */
void tpat_arc_data::setSTMByIx(int ix, MatrixXRd stm){
	if(ix < 0)
		ix += segs.size();

	if(ix < 0 || ix >= (int)(segs.size()))
		throw tpat_exception("tpat_arc_data::setSTMByIx: node index out of bounds");

	segs[ix].setSTM(stm);
}//=================================================

/**
 *	@brief Set the computational tolerance for this data object
 *	@param d the tolerance
 */
void tpat_arc_data::setTol(double d){ tol = d; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy all data from the input arc data to this one
 *	@param d an arc data object reference
 */
void tpat_arc_data::copyMe(const tpat_arc_data &d){
	nodes = d.nodes;
	nodeIDMap = d.nodeIDMap;
	segs = d.segs;
	cons = d.cons;
	sysData = d.sysData; // Copying ADDRESS of sys_data object
	numExtraParam = d.numExtraParam;
	extraParamRowSize = d.extraParamRowSize;
	tol = d.tol;
	nextNodeID = d.nextNodeID;
}//====================================================

/**
 *	@brief Save the acceleration vector to file
 *	@param matFile a pointer to the destination mat-file
 */
void tpat_arc_data::saveAccel(mat_t *matFile) const{
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
void tpat_arc_data::saveEpoch(mat_t *matFile) const{
	saveEpoch(matFile, "Epoch");
}//=====================================================

/**
 *	@brief Save all node epochs to file with a specified variable name
 *	@param matFile a pointer to the destination mat-file
 *	@param varName the name of the variable
 */
void tpat_arc_data::saveEpoch(mat_t *matFile, const char* varName) const{
	std::vector<double> allEpochs(nodes.size());

	for(size_t n = 0; n < nodes.size(); n++){
		allEpochs.push_back(nodes[n].getEpoch());
	}
	
	size_t dims[2] = {nodes.size(), 1};
	matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allEpochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, varName, MAT_COMPRESSION_NONE);
}//=====================================================

/**
 *	@brief Save one of the extra parameters to file
 *	@param matFile a pointer to the destination mat-file
 *	@param varIx the index of the parameter
 *	@param name the name of the variable being saved
 */
void tpat_arc_data::saveExtraParam(mat_t *matFile, int varIx, const char *name) const{
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
void tpat_arc_data::saveState(mat_t *matFile) const{
	saveState(matFile, "State");
}//==================================================

/**
 *	@brief Save the state vector [pos, vel] to a file
 *	@param matFile a pointer to the destination matlab file 
 *	@param varName the name of the variable (e.g. "State" or "Nodes")
 */
void tpat_arc_data::saveState(mat_t *matFile, const char* varName) const{
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
void tpat_arc_data::saveSTMs(mat_t *matFile) const{
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
void tpat_arc_data::saveTOF(mat_t *matFile, const char* varName) const{
	std::vector<double> allTOFs(segs.size());

	for(size_t n = 0; n < segs.size(); n++){
		allTOFs.push_back(segs[n].getTOF());
	}
	
	size_t dims[2] = {segs.size(), 1};
	matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allTOFs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, varName, MAT_COMPRESSION_NONE);
}//=====================================================

/**
 *  @brief Initialize the vectors of node and segment objects from a *.mat file
 *  @details THIS FUNCTION MUST BE THE FIRST READ_DATA-TYPE FUNCTION CALLED because
 *  it clears the vectors and then initializes them by calculating the number
 *  of steps in the arc_data object from the state vector. Individual nodes and segments are
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
void tpat_arc_data::initNodesSegsFromMat(mat_t *matFile, const char* varName){
	matvar_t *stateMat = Mat_VarRead(matFile, varName);
	if(stateMat == NULL){
		throw tpat_exception("tpat_arc_data::initNodeSegsFromMat: Could not read state data vector");
	}else{
		int numSteps = stateMat->dims[0];
		nodes.clear();
		segs.clear();
		tpat_node blank_node;
		tpat_segment blank_seg;
		
		nodes.assign(numSteps, blank_node);	// Initialize array with a bunch of default objects
		segs.assign(numSteps-1, blank_seg);
	}
	Mat_VarFree(stateMat);
}//======================================================

/**
 *  @brief Read the state vector for this arc_data object from a matlab data file
 *  @details This function must be called after initNodeSegsFromMat() as it
 *  populates the step vector objects with state data
 * 
 *  @param matFile pointer to an open matlab data file
 *  @param varName the name of the state variable (e.g., "State" or "Nodes")
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_arc_data::readStateFromMat(mat_t *matFile, const char* varName){
	matvar_t *stateMat = Mat_VarRead(matFile, varName);
	if(stateMat == NULL){
		throw tpat_exception("tpat_arc_data::readStateFromMat: Could not read state data vector");
	}else{
		int numSteps = stateMat->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_arc_data::readStateFromMat: Step vector has not been initialized!");
		}

		if(numSteps != (int)nodes.size()){
			throw tpat_exception("tpat_arc_data::readStateFromMat: State vector has a different size than the initialized step vector");
		}

		if(stateMat->dims[1] != 6){
			throw tpat_exception("tpat_arc_data::readStateFromMat: Incompatible data file: State width is not 6.");
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
			throw tpat_exception("tpat_arc_data::readStateFromMat: Incompatible data file: unsupported data type/class");
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
void tpat_arc_data::readAccelFromMat(mat_t *matFile){
	matvar_t *accelMat = Mat_VarRead(matFile, "Accel");
	if(accelMat == NULL){
		throw tpat_exception("tpat_arc_data::readAccelFromMat: Could not read data vector");
	}else{
		int numSteps = accelMat->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_arc_data::readAccelFromMat: Node vector has not been initialized!");
		}

		if(numSteps != (int)nodes.size()){
			throw tpat_exception("tpat_arc_data::readAccelFromMat: Accel vector has a different size than the initialized node vector");
		}

		if(accelMat->dims[1] != 3){
			throw tpat_exception("tpat_arc_data::readAccelFromMat: Incompatible data file: Accel width is not 3.");
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
			throw tpat_exception("tpat_arc_data::readAccelFromMat: Incompatible data file: unsupported data type/class");
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
void tpat_arc_data::readEpochFromMat(mat_t *matFile, const char* varName){
	matvar_t *epochMat = Mat_VarRead(matFile, varName);
	if(epochMat == NULL){
		throw tpat_exception("tpat_arc_data::readEpochFromMat: Could not read data vector");
	}else{
		int numSteps = epochMat->dims[0];

		if(nodes.size() == 0)
			throw tpat_exception("tpat_arc_data::readEpochFromMat: Node vector has not been initialized");

		if(numSteps != (int)nodes.size())
			throw tpat_exception("tpat_arc_data::readEpochFromMat: Epoch vector has different size than the initialized node evctor");

		if(epochMat->dims[1] != 1)
			throw tpat_exception("tpat_arc_data::readEpochFromMat: Incompatible data file: Epoch vector has more than one column");

		if(epochMat->class_type == MAT_C_DOUBLE && epochMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(epochMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					nodes[i].setEpoch(data[i]);
				}
			}
		}else{
			throw tpat_exception("tpat_arc_data::readEpochFromMat: Incompatible data file: unsupported data type or class");
		}
	}
}//================================================

/**
 *  @brief Read State Transition Matrices from a matlab file
 * 
 *  @param matFile pointer to an open Matlab file
 *  @throws tpat_exception if there are any issues importing the data
 */
void tpat_arc_data::readSTMFromMat(mat_t *matFile){
	matvar_t *allSTM = Mat_VarRead(matFile, "STM");
	if(allSTM == NULL){
		throw tpat_exception("tpat_arc_data::readSTMFromMat: Could not read data vector");
	}else{
		int numSteps = allSTM->dims[2];

		if(segs.size() == 0){
			throw tpat_exception("tpat_arc_data::readSTMFromMat: Step vector has not been initialized!");
		}

		if(numSteps != (int)segs.size()){
			throw tpat_exception("tpat_arc_data::readSTMFromMat: STM vector has a different size than the initialized step vector");
		}

		if(allSTM->dims[0] != 6 || allSTM->dims[1] != 6){
			throw tpat_exception("tpat_arc_data::readSTMFromMat: Incompatible data file: STM is not 6x6.");
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
void tpat_arc_data::readTOFFromMat(mat_t *matFile, const char* varName){
	matvar_t *tofMat = Mat_VarRead(matFile, varName);
	if(tofMat == NULL){
		throw tpat_exception("tpat_arc_data::readTOFFromMat: Could not read data vector");
	}else{
		int numSteps = tofMat->dims[0];

		if(segs.size() == 0)
			throw tpat_exception("tpat_arc_data::readTOFFromMat: Node vector has not been initialized");

		if(numSteps != (int)segs.size())
			throw tpat_exception("tpat_arc_data::readTOFFromMat: Epoch vector has different size than the initialized segment evctor");

		if(tofMat->dims[1] != 1)
			throw tpat_exception("tpat_arc_data::readTOFFromMat: Incompatible data file: Epoch vector has more than one column");

		if(tofMat->class_type == MAT_C_DOUBLE && tofMat->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(tofMat->data);

			if(data != NULL){
				for(int i = 0; i < numSteps; i++){
					segs[i].setTOF(data[i]);
				}
			}
		}else{
			throw tpat_exception("tpat_arc_data::readTOFFromMat: Incompatible data file: unsupported data type or class");
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
void tpat_arc_data::readExtraParamFromMat(mat_t *matFile, int varIx, const char *varName){
	if(varIx > numExtraParam || varIx < 0)
		throw tpat_exception("tpat_arc_data::readExtraParamFromMat: Could not read extra parameter; index out of bounds");

	// Get starting index of this extra param within a arc step's extra parameter vector
	int ix0 = 0;
	for(int i = 0; i < varIx; i++){ ix0 += extraParamRowSize[i]; }

	matvar_t *matvar = Mat_VarRead(matFile, varName);
	if(matvar == NULL){
		throw tpat_exception("tpat_arc_data::readExtraParamFromMat: Could not read data vector");
	}else{
		int numSteps = matvar->dims[0];
		
		if(nodes.size() == 0){
			throw tpat_exception("tpat_arc_data::readExtraParamFromMat: Step vector has not been initialized!");
		}

		if(matvar->dims[1] != ((size_t)extraParamRowSize[varIx])){
			char message[64];
			sprintf(message, "tpat_arc_data::readExtraParamFromMat: Incompatible data file: %s width is not %d", varName, extraParamRowSize[varIx]);
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
			throw tpat_exception("tpat_arc_data::readExtraParamFromMat: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//===============================================

/**
 *	@brief Update the constraints for every node so that their node numberes
 * 	match the node/step they belong to.
 */
void tpat_arc_data::updateCons(){
	for(size_t n = 0; n < nodes.size(); n++){
		std::vector<tpat_constraint> nodeCons = nodes[n].getConstraints();
		for(size_t c = 0; c < nodeCons.size(); c++){
			nodeCons[c].setNode(n);
		}
		nodes[n].setConstraints(nodeCons);
	}
}//====================================================




