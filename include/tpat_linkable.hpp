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

#ifndef H_TPAT_LINKABLE
#define H_TPAT_LINKABLE

#include "tpat.hpp"
// Forward Declarations


/**
 *	@brief A basic object that encapsulates the behavior of nodes and segments that hold
 *	information about how they link together
 *
 *	@author Andrew Cox
 *	@version May 1, 2016
 *	@copyright GNU GPL v3.0
 */
class TPAT_Linkable : public TPAT{
public:
	static const int INVALID_ID;	//!< Reserved ID value
	static const int NUM_LINKS;		//!< Number of links stored by a linkable object

	// *structors
	TPAT_Linkable();
	TPAT_Linkable(const TPAT_Linkable&);
	virtual ~TPAT_Linkable();
	
	// Operators
	TPAT_Linkable& operator =(const TPAT_Linkable&);
	friend bool operator ==(const TPAT_Linkable&, const TPAT_Linkable&);
	friend bool operator !=(const TPAT_Linkable&, const TPAT_Linkable&);
	
	// Set and Get Functions
	void addLink(int);
	void clearLinks();
	int getID() const;
	int getLink(int) const;
	int getLinkIx(int) const;
	bool isLinkedTo(int) const;
	void removeLink(int);
	void setID(int); 
	void setLink(int, int);
	
protected:
	int ID = INVALID_ID;	//!< The ID associated with this object
	
	/**
	 * An array that stores up to two links in the form of the IDs of 
	 * the objects that this one is linked to.
	 */
	int links[2] = {INVALID_ID, INVALID_ID};

	void copyMe(const TPAT_Linkable&);
};

#endif