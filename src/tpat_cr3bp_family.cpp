/**
 *	@file tpat_cr3bp_family.cpp
 *	@brief Data object for a CR3BP family
 */
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

#include "tpat_cr3bp_family.hpp"

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_cr3bp_traj.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Create an empty family for the specified system
 */
tpat_cr3bp_family::tpat_cr3bp_family(tpat_cr3bp_sys_data data){
	sysData = data;
}//====================================================

/**
 *	@brief load a family from a file
 */
tpat_cr3bp_family::tpat_cr3bp_family(std::string filepath){
	throw tpat_exception("Not yet implemented!");
}//====================================================

/**
 *	@brief Copy Constructor
 */
tpat_cr3bp_family::tpat_cr3bp_family(const tpat_cr3bp_family& fam){
	copyMe(fam);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_cr3bp_family::~tpat_cr3bp_family(){}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 */
tpat_cr3bp_family& tpat_cr3bp_family::operator= (const tpat_cr3bp_family& fam){
	copyMe(fam);
	return *this;
}


//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a member to the family
 *	@param mem a new family member; NOTE: <tt>mem</tt> should represent
 *	a trajectory that exists in the same system as the other family members.
 */
void tpat_cr3bp_family::addMember(tpat_cr3bp_family_member mem){
	members.push_back(mem);
}//====================================================

/**
 *	@brief Retrieve a family member by its index
 *	@param ix the index of the member. If the index is less than zero,
 *	it will count backwards from the end.
 */
tpat_cr3bp_family_member tpat_cr3bp_family::getMember(int ix) const{
	if(ix < 0)
		ix += members.size();

	return members.at(ix);
}//====================================================

/**
 *	@brief Locate all members with the specified value of Jacobi Constant
 *	@param jc the desired value for Jacobi
 *	@return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<tpat_cr3bp_family_member> tpat_cr3bp_family::getMemberByJacobi(double jc) const{
	// Get an array of all the jacobi values
	std::vector<double> allJC;
	for(int n = 0; n < ((int)members.size()); n++){
		allJC.push_back(members[n].getJC());
	}

	tpat_constraint jacobiCon(tpat_constraint::JC, 0, &jc, 1);

	return getMatchingMember(jc, &allJC, jacobiCon);
}//==============================================

/**
 *	@brief Locate all members with the specified time-of-flight
 *	@param tof the desired value for time-of-flight
 *	@return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<tpat_cr3bp_family_member> tpat_cr3bp_family::getMemberByTOF(double tof) const{
	std::vector<double> allTOF;
	for(int n = 0; n < ((int)members.size()); n++){
		allTOF.push_back(members[n].getTOF());
	}

	tpat_constraint tofCon(tpat_constraint::TOF, 0, &tof, 1);

	return getMatchingMember(tof, &allTOF, tofCon);
}//==============================================

/**
 *	@brief Locate all members with the specified value for the specified state variable in the IC
 *	@param value the desired value
 *	@param ix the index of the state variable [0, 5]
 *	@return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<tpat_cr3bp_family_member> tpat_cr3bp_family::getMemberByStateVar(double value, int ix) const{
	if(ix < 0 || ix > 5){
		throw tpat_exception("Invalid state index; out of range");
	}

	std::vector<double> allVals;
	for(int n = 0; n < ((int)members.size()); n++){
		allVals.push_back(members[n].getIC()[ix]);
	}

	double conData[] = {NAN,NAN,NAN,NAN,NAN,NAN};
	conData[ix] = value;
	tpat_constraint stateCon(tpat_constraint::STATE, 0, conData, 6);

	return getMatchingMember(value, &allVals, stateCon);
}//==============================================

/**
 *	@brief Retrieve the name of this family
 *	@return a descriptive name
 */
std::string tpat_cr3bp_family::getName() const { return name; }

/**
 *	@brief Determine which variable best naturally describes the flow of the family
 *	@return the sorting variable
 */
tpat_cr3bp_family::sortVar_t tpat_cr3bp_family::getSortType() const { return sortType; }

/**
 *	@brief Retrieve the system data for this family
 *	@return the system data describing the system all family members exist in
 */
tpat_cr3bp_sys_data tpat_cr3bp_family::getSysData() const { return sysData; }

/**
 *	@brief Set the name of this family
 *	@param n a descriptive name
 */
void tpat_cr3bp_family::setName(std::string n){ name = n; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy the family from another family
 *	@param fam a different family
 */
void tpat_cr3bp_family::copyMe(const tpat_cr3bp_family& fam){
	members = fam.members;
	sysData = fam.sysData;
	name = fam.name;
	sortType = fam.sortType;
}//====================================================

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
std::vector<int> tpat_cr3bp_family::findMatches(double value, std::vector<double> *data) const{
	std::vector<int> matches;
	for(int n = 0; n < ((int)data->size()); n++){
		if(std::abs(data->at(n) - value) < matchTol){
			matches.push_back(n);
		}else if(n > 0 && (data->at(n) - value)*(data->at(n-1) - value) < 0){
			matches.push_back(n-1);
		}
	}

	return matches;
}//=====================================================

/**
 *	@brief Locate a family member with a specific attribute
 *	
 *	This function locates a family member or set of members that have a specific value
 *	for one of the variables of interest (e.g. coordinates, Jacobi, TOF). Exact matches
 *	and interpolated matches are returned; interpolated matches are computed using a
 *	differential corrections algorithm.
 *
 *	@param value the value the family member should have
 *	@param dataSet a pointer to a vector containing the set of values to search for matches
 *	in. For example, if the <tt>value</tt> I pass in contains a specific TOF, then 
 *	<tt>dataSet</tt> points to a vector containing the TOFs for the entire family, sorted
 *	according to this family's <tt>sortType</tt>.
 *	@param matchCon a constraint that can be applied in a corrections scheme that will
 *	ensure the corrected trajectory has the desired value for the variable of interest.
 *
 *	@return a vector of matches. If no matches are returned, the vector will be empty.
 */
std::vector<tpat_cr3bp_family_member> tpat_cr3bp_family::getMatchingMember(double value, std::vector<double> *dataSet,
	tpat_constraint matchCon) const{
	// Locate possible candidates
	std::vector<int> matches = findMatches(value, dataSet);
	if(matches.size() == 0){
		throw tpat_exception("Could not locate any matches. The family either has too few members to facilitate an accurate search or the desired trajectory does not exist.");
	}

	std::vector<tpat_cr3bp_family_member> matchMembers;
	for(int n = 0; n < ((int)matches.size()); n++){
		// Check to see if they are "close enough"
		if(std::abs(dataSet->at(matches[n]) - value) < matchTol){
			matchMembers.push_back(members[matches[n]]);
		}else{	// If not, employ corrections
			
			// Create a nodeset and a constraint to make the orbit periodic
			tpat_cr3bp_nodeset memberSet(members[n].getIC(), sysData, members[n].getTOF(), numNodes);
			double conData[] = {0,0,0,0,0,NAN};
			tpat_constraint periodicCon(tpat_constraint::MATCH_CUST, numNodes-1, conData, 6);

			memberSet.addConstraint(periodicCon);
			memberSet.addConstraint(matchCon);

			// Correct the nodeset while constraining the orbit to have the desired characterstic
			tpat_correction_engine corrector;
			try{
				corrector.correct_cr3bp(&memberSet);
				tpat_cr3bp_nodeset newNodes = corrector.getCR3BPOutput();

				tpat_cr3bp_traj newTraj = tpat_cr3bp_traj::fromNodeset(newNodes);
				tpat_cr3bp_family_member newMember(newTraj);
				matchMembers.push_back(newMember);
			}catch(tpat_diverge &e){}
		}
	}

	return matchMembers;
}//==============================================

/**
 *	@brief Populate an array with a single coordinate from each family member
 *	@param ix the index of the coordinate within the IC vector
 *	@param array a pointer to the array we want to populate
 */
void tpat_cr3bp_family::getCoord(int ix, std::vector<double> *array) const{
	for(int i = 0; i < ((int)members.size()); i++){
		array->push_back(members[i].getIC()[ix]);
	}
}//====================================================

/**
 *	@brief Sort the family members by the specified sort variable (in ascending order)
 *	
 *	The sorting variable is specified by <tt>sortType</tt>; this is the variable
 *	that best describes the natural progression of the family. For example,
 *	Lyapunov orbits can be evolved naturally by varying the x-coordinate of the IC.
 *
 *	The family must be sorted (or be loaded from an already sorted family file) before
 *	you can retrieve family members. The process will run without sorting, but the results
 *	will likely be wonky.
 */
void tpat_cr3bp_family::sort(){
	// Create an array containing the independent variable from each family member
	std::vector<double> dataToSort;
	switch(sortType){
		case SORT_X:
			getCoord(0, &dataToSort); break;
		case SORT_Y:
			getCoord(1, &dataToSort); break;
		case SORT_Z:
			getCoord(2, &dataToSort); break;
		case SORT_VX:
			getCoord(3, &dataToSort); break;
		case SORT_VY:
			getCoord(4, &dataToSort); break;
		case SORT_VZ:
			getCoord(5, &dataToSort); break;
		case SORT_JC:
			for(int i = 0; i < ((int)members.size()); i++){
				dataToSort.push_back(members[i].getJC());
			}
			break;
		case SORT_TOF:
			for(int i = 0; i < ((int)members.size()); i++){
				dataToSort.push_back(members[i].getTOF());
			}
			break;
		default:
			throw tpat_exception("Unrecognized sorting type");
			break;
	}

	// Sort the data, retrieve the indices of the now-sorted elements
	std::vector<int> indices = getSortedInd(dataToSort);

	// Use those indices to sort the family members
	std::vector<tpat_cr3bp_family_member> sortedMembers;
	for(int n = 0; n < ((int)indices.size()); n++){
		sortedMembers.push_back(members[indices[n]]);
	}

	// Reassign the members array
	members = sortedMembers;
}//===================================================

