/**
 *  @file Linkable.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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
#pragma once

#include "Core.hpp"
// Forward Declarations

namespace astrohelion{
	/**
	 *	@brief A basic object that encapsulates the behavior of nodes and segments that hold
	 *	information about how they link together
	 *
	 *	@author Andrew Cox
	 *	@version May 1, 2016
	 *	@copyright GNU GPL v3.0
	 */
	class Linkable : public Core{
	public:
		static const int INVALID_ID;	//!< Reserved ID value
		static const int NUM_LINKS;		//!< Number of links stored by a linkable object

		/**
		 *  \name *structors
		 *  \{
		 */
		Linkable();
		Linkable(const Linkable&);
		virtual ~Linkable();
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		Linkable& operator =(const Linkable&);
		friend bool operator ==(const Linkable&, const Linkable&);
		friend bool operator !=(const Linkable&, const Linkable&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		int getID() const;
		int getLink(int) const;
		int getLinkIx(int) const;
		bool isLinkedTo(int) const;
		
		virtual void setID(int); 
		virtual void setLink(int, int);
		//\}
		
		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		virtual void addLink(int);
		void clearLinks();
		virtual void removeLink(int);
		//\}
	protected:
		int ID = INVALID_ID;	//!< The ID associated with this object
		
		/**
		 * An array that stores up to two links in the form of the IDs of 
		 * the objects that this one is linked to.
		 */
		int links[2] = {INVALID_ID, INVALID_ID};

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void copyMe(const Linkable&);
		//\}
	};
}