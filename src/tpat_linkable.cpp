/**
 *  @file tpat_linkable.cpp
 *	@brief A basic object that encapsulates the behavior of objects that link together
 *
 *	@author Andrew Cox
 *	@version May 1, 2016
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
#include "tpat_linkable.hpp"

#include "tpat_exceptions.hpp"

#include <algorithm>

const int TPAT_Linkable::INVALID_ID = -1;
const int TPAT_Linkable::NUM_LINKS = 2;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Default constructor
 */
TPAT_Linkable::TPAT_Linkable(){}

/**
 *  @brief Copy constructor
 * 
 *  @param obj reference to a linkable object
 */
TPAT_Linkable::TPAT_Linkable(const TPAT_Linkable &obj){
	copyMe(obj);
}//====================================================

/**
 *  @brief Destructor
 */
TPAT_Linkable::~TPAT_Linkable(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  @brief Assignment operator
 * 
 *  @param obj reference to a linkable object
 *  @return reference to a linkable object, now equivalent to the input reference
 */
TPAT_Linkable& TPAT_Linkable::operator =(const TPAT_Linkable &obj){
	copyMe(obj);
	return *this;
}//====================================================

/**
 *  @brief Equality comparison
 * 
 *  @param lhs 
 *  @param rhs 
 * 
 *  @return whether or not the two objects are identical
 */
bool operator ==(const TPAT_Linkable &lhs, const TPAT_Linkable &rhs){
	for(int i = 0; i < TPAT_Linkable::NUM_LINKS; i++){
		if(lhs.links[i] != rhs.links[i])
			return false;
	}

	return lhs.ID == rhs.ID;
}//====================================================

/**
 *  @brief Inequality comparison
 * 
 *  @param lhs 
 *  @param rhs 
 * 
 *  @return whether or not the two objects are nonidentical
 */
bool operator !=(const TPAT_Linkable &lhs, const TPAT_Linkable &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Add a link from this object to another with the specified ID
 *  @details Only two links are allowed; if both slots are already filled,
 *  an exception is thrown
 * 
 *  @param id the ID of the object that is linked to this object
 *  @throws TPAT_Exception if two links are already present
 */
void TPAT_Linkable::addLink(int id){
	int ix = getLinkIx(INVALID_ID);	// Try to find an empty slot
	if(ix != INVALID_ID){
		setLink(ix, id);	// This function will also check for duplicate ID in the link slots
	}else{
		throw TPAT_Exception("linkable::addLink: Two links already present, cannot add a new one");
	}
}//===========================================

/**
 *  @brief Remove all links from this object to other objects
 *  @details This function resets the array of links
 */
void TPAT_Linkable::clearLinks(){
	for(int i = 0; i < NUM_LINKS; i++)
		links[i] = INVALID_ID;
}//===========================================

/**
 *  @brief Retrieve the ID of this object
 *  @return the ID of this object
 */
int TPAT_Linkable::getID() const { return ID; }

/**
 *  @brief Locate a particular ID within the link slots associated with this object
 * 
 *  @param id the ID to search for
 *  @return the index of the ID within the links[] storage array. A value of 
 *  TPAT_Linkable::INVALID_ID is returned if the ID is not located
 */
int TPAT_Linkable::getLinkIx(int id) const{	
	const int *it = std::find(links, links+NUM_LINKS, id);
	
	if(it == links+NUM_LINKS){
		return INVALID_ID;
	}
	else{
		return it - links;
	}
}//============================================

/**
 *  @brief Retrieve one of the links
 * 
 *  @param ix The index of the link, must be in the integer set [0, linkable::NUM_LINKS)
 *  @return the ID stored at the specified index. If the index is out of 
 *  range, an exception is thrown
 *  @throws TPAT_Exception if <tt>ix</tt> is out of bounds
 */
int TPAT_Linkable::getLink(int ix) const {
	if(ix >= 0 && ix < NUM_LINKS){
		return links[ix];
	}else{
		throw TPAT_Exception("linkable::getLink: Index out of bounds");
	}
}//===========================================

/**
 *  @brief Determine whether an object with the specified ID is linked
 *  to this object
 * 
 *  @param id the ID of the object
 *  @return whether or not the object with the specified ID is linked
 *  to this object
 */
bool TPAT_Linkable::isLinkedTo(int id) const {
	return std::find(links, links + NUM_LINKS, id) != links+NUM_LINKS;
}//===========================================

/**
 *  @brief Remove a link to an object with the specified ID.
 *  If the object with the specified ID is not linked to this
 *  object, no further action is taken
 * 
 *  @param id the ID of the object to "unlink"
 */
void TPAT_Linkable::removeLink(int id){
	int *it = std::find(links, links + NUM_LINKS, id);
	if(it != links + NUM_LINKS)
		links[it-links] = INVALID_ID;	
}//===========================================

/**
 *  @brief Set the ID of this object
 *  @details Note: This should NOT be done after a segment is added to an arcset
 * 
 *  @param id the desired ID
 */
void TPAT_Linkable::setID(int id){ ID = id; } 

/**
 *  @brief Set the ID stored in the specified link
 * 
 *  @param ix the index of the link within the storage array; must be
 *  within the integer set [0, linkable::NUM_LINKS). Invalid indices will
 *  result in an exception
 *  @param id the ID to store; you cannot store the same ID twice; an exception
 *  will be thrown if you try
 *  @throws TPAT_Exception if <tt>ix</tt> is out of bounds
 *  @throws TPAT_Exception if a link already exists to <tt>id</tt>
 */
void TPAT_Linkable::setLink(int ix, int id) {
	if(ix >= 0 && ix < NUM_LINKS){
		links[ix] = id;
		if(links[(ix+1) % NUM_LINKS] == id && id != INVALID_ID){
			throw TPAT_Exception("linkable::setLink: Cannot link both slots to same ID");
		}
	}else{
		throw TPAT_Exception("linkable::setLink: Index out of bounds");
	}
}//===========================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  @brief Copy all aspects of the reference object into this one
 * 
 *  @param L reference to a linkable object
 */
void TPAT_Linkable::copyMe(const TPAT_Linkable &L){
	ID = L.ID;
	std::copy(L.links, L.links + NUM_LINKS, links);
}//=====================================================