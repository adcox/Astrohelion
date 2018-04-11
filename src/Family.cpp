/**
 *  @file Family.cpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version October 5, 2017
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

#include "Family.hpp"

#include <cmath>

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors and Desctructor
//-----------------------------------------------------------------------------

/**
 * @brief Default constructor, requires system data pointer
 * 
 * @param pSys pointer to system data object common to all family members
 */
Family::Family(const SysData *pSys) : pSysData(pSys){}

/**
 * @brief Construct a Family object from another object
 * 
 * @param f reference to a Family object
 */
Family::Family(const Family &f) : pSysData(f.pSysData){
	copyMe(f);
}//====================================================

/**
 * @brief Default destructor; no functionality at this time
 */
Family::~Family() {}

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

/**
 * @brief Copy constructor
 * 
 * @param f Reference to a Family object
 * @return a reference to this Family object, now equal to the input `f`
 */
Family& Family::operator= (const Family &f){
	copyMe(f);
	return *this;
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 * @brief Retrieve the tolerance used when finding matches between user-
 * specified values and stored data
 * 
 * @return the numerical tolerance for detecting matches. A difference less
 * than this magnitude is considered zero.
 * @see setMatchTol()
 */
double Family::getMatchTol() const { return matchTol; }

/**
 * @brief Retrieve the name of the family
 * @return the name of the family
 * @see setName()
 */
std::string Family::getName() const{ return name; }

/**
 * @brief Retrieve the sortType associated with this family
 * @return a value of FamSort_tp that dictates how the family members are sorted
 * @see setSortType()
 */
FamSort_tp Family::getSortType() const { return sortType; }

/**
 * @brief Retrieve a string that describes the sort type associated with this 
 * family
 * @return a string that describes the sort type
 * @see getSortType()
 */
const char* Family::getSortTypeStr() const { return sortTypeToStr(sortType); }

/**
 * @brief Retrieve the system data pointer associated with this family
 * @return the system data pointer associated with this family
 */
const SysData* Family::getSysData() const { return pSysData; }

/**
 * @brief Set the match tolerance for the family
 * 
 * @param tol a tolerance used to determine which numerical values are 
 * practically zero
 * @see getMatchTol()
 */
void Family::setMatchTol(double tol){ matchTol = tol; }

/**
 * @brief Set the name of the family
 * 
 * @param name string describing the family
 * @see getName()
 */
void Family::setName(std::string name){ this->name = name; }

/**
 * @brief Set the family sort type
 * 
 * @param tp describes how family members are sorted when sortMembers() is
 * called
 * @see getSortType()
 */
void Family::setSortType(FamSort_tp tp){ sortType = tp; }

//-----------------------------------------------------------------------------
//      Family Operations
//-----------------------------------------------------------------------------

/**
 *	@brief Locate places in a data set where a specific value probably exists
 *
 *	This algorithm will locate both exact matches (within a tolerance) and intersections,
 *	assuming the data is continuous. If an intersection is found, the index of the point
 *	before the intersection is returned.
 *
 *	@param value the value to search for
 *	@param data a pointer to a data set to search in
 *	@return a vector of integers representing the indices of matches
 */
std::vector<unsigned int> Family::findMatches(double value, 
	std::vector<double> *data) const{

	double numBins = data->size() > 500 ? 100 : (data->size()/5.0);
	int binSize = std::floor((data->size())/numBins);

	std::vector<unsigned int> matches;
	double minDiff = 0;
	double diff = 0;
	unsigned int minDiffIx = 0;
	for(unsigned int n = 0; n < data->size(); n++){
		diff = std::abs(data->at(n) - value);

		// Search for acceptable minema in bins
		if(n % binSize == 0){
			// reset
			minDiff = diff;
			minDiffIx = n;

			if(minDiff < matchTol)
				matches.push_back(minDiffIx);
		}else{
			if(diff < minDiff){
				minDiff = diff;
				minDiffIx = n;
			}
		}

		// Search for intersections
		if(n > 0 && (data->at(n) - value)*(data->at(n-1) - value) < 0){
			matches.push_back(n-1);
		}
	}

	return matches;
}//=====================================================

//-----------------------------------------------------------------------------
//      Utility
//-----------------------------------------------------------------------------

/**
 * @brief Convert the numerical type into a human-readable string
 * 
 * @param tp family sort type
 * @return human-readable string version of the numerical type
 */
const char* Family::sortTypeToStr(FamSort_tp tp){
	switch(tp){
		case FamSort_tp::SORT_X: return "SORT_X"; break;
		case FamSort_tp::SORT_Y: return "SORT_Y"; break;
		case FamSort_tp::SORT_Z: return "SORT_Z"; break;
		case FamSort_tp::SORT_VX: return "SORT_VX"; break;
		case FamSort_tp::SORT_VY: return "SORT_VY"; break;
		case FamSort_tp::SORT_VZ: return "SORT_VZ"; break;
		case FamSort_tp::SORT_JC: return "SORT_JC"; break;
		case FamSort_tp::SORT_TOF: return "SORT_TOF"; break;
		case FamSort_tp::SORT_NONE: return "NO SORTING"; break;
		default: return "Unrecognized Type!";
	}
}//====================================================

/**
 * @brief Copy all member data between the input Family object and this family
 * 
 * @param f reference to another Family object
 */
void Family::copyMe(const Family &f){
	name = f.name;
	sortType = f.sortType;
	matchTol = f.matchTol;
}//====================================================

}// End of astrohelion namespace