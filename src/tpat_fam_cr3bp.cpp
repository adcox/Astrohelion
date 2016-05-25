/**
 *	@file tpat_fam_cr3bp.cpp
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
 
#include "tpat_fam_cr3bp.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

#include <algorithm>
#include <cstring>
//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Create an empty family for the specified system
 *	@param data a CR3BP system data object
 */
TPAT_Fam_CR3BP::TPAT_Fam_CR3BP(TPAT_Sys_Data_CR3BP data){
	sysData = data;
}//====================================================

/**
 *	@brief load a family from a file
 *	@param filepath an aboslute or relative filepath to the data file
 *	@throws TPAT_Exception if the file cannot be opened
 */
TPAT_Fam_CR3BP::TPAT_Fam_CR3BP(const char* filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw TPAT_Exception("TPAT_Fam_CR3BP: Could not load family from file");
	}

	loadMemberData(matfp);
	loadEigVals(matfp);
	name = readStringFromMat(matfp, NAME_VAR_NAME, MAT_T_UINT8, MAT_C_CHAR);
	double type = readDoubleFromMat(matfp, SORT_TYPE_VAR_NAME);
	sortType = static_cast<TPAT_SortFam_Tp>((int)type);
	// sysData.readFromMat(matfp);

	Mat_Close(matfp);
	sysData = TPAT_Sys_Data_CR3BP(filepath);
}//====================================================

/**
 *	@brief Copy Constructor
 *	@brief fam a family reference
 */
TPAT_Fam_CR3BP::TPAT_Fam_CR3BP(const TPAT_Fam_CR3BP& fam){
	copyMe(fam);
}//====================================================

/**
 *	@brief Destructor
 */
TPAT_Fam_CR3BP::~TPAT_Fam_CR3BP(){}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param fam a family reference
 */
TPAT_Fam_CR3BP& TPAT_Fam_CR3BP::operator= (const TPAT_Fam_CR3BP& fam){
	copyMe(fam);
	return *this;
}//====================================================


//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a member to the family
 *	@param mem a new family member; NOTE: <tt>mem</tt> should represent
 *	a trajectory that exists in the same system as the other family members.
 */
void TPAT_Fam_CR3BP::addMember(TPAT_FamMember_CR3BP mem){
	members.push_back(mem);
}//====================================================

/**
 *	@brief Retrieve a family member by its index
 *	@param ix the index of the member. If the index is less than zero,
 *	it will count backwards from the end.
 */
TPAT_FamMember_CR3BP TPAT_Fam_CR3BP::getMember(int ix) const{
	if(ix < 0)
		ix += members.size();

	return members.at(ix);
}//====================================================

/**
 *	@brief Locate all members with the specified value of Jacobi Constant
 *	@param jc the desired value for Jacobi
 *	@return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<TPAT_FamMember_CR3BP> TPAT_Fam_CR3BP::getMemberByJacobi(double jc) const{
	// Get an array of all the jacobi values
	std::vector<double> allJC;
	for(int n = 0; n < ((int)members.size()); n++){
		allJC.push_back(members[n].getJacobi());
	}

	TPAT_Constraint jacobiCon(TPAT_Constraint_Tp::JC, 0, &jc, 1);

	return getMatchingMember(jc, &allJC, jacobiCon);
}//==============================================

/**
 *	@brief Locate all members with the specified time-of-flight
 *	@param tof the desired value for time-of-flight
 *	@return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<TPAT_FamMember_CR3BP> TPAT_Fam_CR3BP::getMemberByTOF(double tof) const{
	std::vector<double> allTOF;
	for(size_t n = 0; n < members.size(); n++){
		allTOF.push_back(members[n].getTOF());
	}

	TPAT_Constraint tofCon(TPAT_Constraint_Tp::TOF, 0, &tof, 1);

	return getMatchingMember(tof, &allTOF, tofCon);
}//==============================================

/**
 *	@brief Locate all members with the specified value for the specified state variable in the IC
 *	@param value the desired value
 *	@param ix the index of the state variable [0, 5]
 *	@return a vector of matching family members (some may be interpolated/corrected)
 *	@throws TPAT_Exception if <tt>ix</tt> is out of range
 */
std::vector<TPAT_FamMember_CR3BP> TPAT_Fam_CR3BP::getMemberByStateVar(double value, int ix) const{
	if(ix < 0 || ix > 5){
		throw TPAT_Exception("TPAT_Fam_CR3BP::getMemberByStateVar: Invalid state index; out of range");
	}

	std::vector<double> allVals;
	for(size_t n = 0; n < members.size(); n++){
		allVals.push_back(members[n].getIC()[ix]);
	}

	double conData[] = {NAN,NAN,NAN,NAN,NAN,NAN};
	conData[ix] = value;
	TPAT_Constraint stateCon(TPAT_Constraint_Tp::STATE, 0, conData, 6);

	return getMatchingMember(value, &allVals, stateCon);
}//==============================================

/**
 *	@brief Retrieve the name of this family
 *	@return a descriptive name
 */
std::string TPAT_Fam_CR3BP::getName() const { return name; }

/**
 *	@brief Retrieve number of members in this family
 *	@return number of members in this family
 */
int TPAT_Fam_CR3BP::getNumMembers() const { return members.size(); }

/**
 *	@brief Determine which variable best naturally describes the flow of the family
 *	@return the sorting variable
 */
TPAT_SortFam_Tp TPAT_Fam_CR3BP::getSortType() const { return sortType; }

/**
 *	@brief Retrieve a string describing the sort type in human-readable format
 *	@return a string describing the sort type in human-readable format
 */
const char* TPAT_Fam_CR3BP::getSortTypeStr() const{
	switch(sortType){
		case TPAT_SortFam_Tp::SORT_X: return "SORT_X"; break;
		case TPAT_SortFam_Tp::SORT_Y: return "SORT_Y"; break;
		case TPAT_SortFam_Tp::SORT_Z: return "SORT_Z"; break;
		case TPAT_SortFam_Tp::SORT_VX: return "SORT_VX"; break;
		case TPAT_SortFam_Tp::SORT_VY: return "SORT_VY"; break;
		case TPAT_SortFam_Tp::SORT_VZ: return "SORT_VZ"; break;
		case TPAT_SortFam_Tp::SORT_JC: return "SORT_JC"; break;
		case TPAT_SortFam_Tp::SORT_TOF: return "SORT_TOF"; break;
		case TPAT_SortFam_Tp::SORT_NONE: return "NO SORTING"; break;
		default: return "Unrecognized Type!";
	}
}//===========================================

/**
 *	@brief Retrieve the system data for this family
 *	@return the system data describing the system all family members exist in
 */
TPAT_Sys_Data_CR3BP TPAT_Fam_CR3BP::getSysData() const { return sysData; }

/**
 *	@brief Retrieve a pointer to the system data object for this family
 *	@return a pointer to the system data object for this family
 */
TPAT_Sys_Data_CR3BP* TPAT_Fam_CR3BP::getSysDataPtr() { return &sysData; }

/**
 *	@brief Set the name of this family
 *	@param n a descriptive name; must be less than 128 characters
 *	if you want to save to a matlab file
 */
void TPAT_Fam_CR3BP::setName(std::string n){ name = n; }

/**
 *	@brief Set the sort type for this family
 *	@param type the sort type
 */
void TPAT_Fam_CR3BP::setSortType(TPAT_SortFam_Tp type){ sortType = type; }
//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Copy the family from another family
 *	@param fam a different family
 */
void TPAT_Fam_CR3BP::copyMe(const TPAT_Fam_CR3BP& fam){
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
std::vector<int> TPAT_Fam_CR3BP::findMatches(double value, std::vector<double> *data) const{
	double numBins = data->size() > 500 ? 100 : (data->size()/5.0);
	int binSize = floor((data->size())/numBins);

	std::vector<int> matches;
	double minDiff = 0;
	double diff = 0;
	int minDiffIx = 0;
	for(size_t n = 0; n < data->size(); n++){
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
std::vector<TPAT_FamMember_CR3BP> TPAT_Fam_CR3BP::getMatchingMember(double value, std::vector<double> *dataSet,
	TPAT_Constraint matchCon) const{
	// Locate possible candidates
	std::vector<int> matches = findMatches(value, dataSet);
	std::vector<TPAT_FamMember_CR3BP> matchMembers;
	if(matches.size() == 0){
		printErr("Could not locate any matches. The family either has too few members to facilitate an accurate search or the desired trajectory does not exist.\n");
		return matchMembers;	// empty set
	}else{
		printColor(GREEN, "Located %zu matches; applying corrections\n", matches.size());
	}

	TPAT_Sys_Data_CR3BP tempSys = sysData;
	
	for(int n = 0; n < ((int)matches.size()); n++){
		int idx = matches[n];

		// Check to see if they are "close enough"
		if(std::abs(dataSet->at(idx) - value) < matchTol){
			matchMembers.push_back(members[idx]);
			printf("  Candidate %d is close enough; no corrections required\n", n);
		}else{	// If not, employ corrections
			printf("  Correcting candidate %d...\n", n);
			// Create a nodeset and a constraint to make the orbit periodic
			double tof = members[idx].getTOF();
			int numNodes = tof > 2 ? floor(tof) : 2;

			TPAT_Nodeset_CR3BP memberSet(members[idx].getIC(), &tempSys, tof, numNodes);
			double end = numNodes-1;
			double conData[] = {end,end,end,end,end,end};
			TPAT_Constraint periodicCon(TPAT_Constraint_Tp::MATCH_CUST, 0, conData, 6);

			memberSet.addConstraint(periodicCon);
			memberSet.addConstraint(matchCon);

			// Correct the nodeset while constraining the orbit to have the desired characterstic
			TPAT_Correction_Engine corrector;
			corrector.setTol(1e-11);
			try{
				TPAT_Nodeset_CR3BP newNodes(&tempSys);
				corrector.multShoot(&memberSet, &newNodes);
				
				TPAT_Traj_CR3BP newTraj = TPAT_Traj_CR3BP::fromNodeset(newNodes);
				TPAT_FamMember_CR3BP newMember(newTraj);
				matchMembers.push_back(newMember);
			}catch(TPAT_Diverge &e){}
		}
	}

	return matchMembers;
}//==============================================

/**
 *	@brief Populate an array with a single coordinate from each family member
 *	@param ix the index of the coordinate within the IC vector
 *	@param array a pointer to the array we want to populate
 */
void TPAT_Fam_CR3BP::getCoord(int ix, std::vector<double> *array) const{
	for(size_t i = 0; i < members.size(); i++){
		std::vector<double> ic = members[i].getIC();
		array->push_back(ic[ix]);
	}
}//====================================================

/**
 *	@brief Locate all bifurcations in the family by analyzing the 
 *	eigenvalues
 *
 *	Eigenvalues MUST be sorted, or this will yield completely bogus
 *	results
 */
std::vector<int> TPAT_Fam_CR3BP::findBifurcations(){
	// TODO: Incorporate much more advanced ways to compute this
	double okErr = 1e-3;
	cdouble one(1,0);

	TPAT_EigValSet_Tp setTypes[3];

	// Determine what "type" of eigenvalue pair we have
	for(int set = 0; set < 3; set++){
		double sumImag = 0;
		// double sumReal = 0;
		double sumDistFromOne = 0;

		for(size_t m = 0; m < members.size(); m++){
			std::vector<cdouble> eigs = members[m].getEigVals();
			sumImag += (std::abs(std::imag(eigs[set*2])) + std::abs(std::imag(eigs[set*2+1])))/2;
			sumDistFromOne += (std::abs(eigs[set*2] - one) + std::abs(eigs[set*2+1] - one))/2;
		}

		double count = members.size();
		double meanImag = sumImag/count;
		// double meanReal = sumReal/count;
		double meanDistFromOne = sumDistFromOne/count;

		if(meanImag > okErr){
			// significant imaginary parts
			setTypes[set] = TPAT_EigValSet_Tp::EIGSET_COMP_CONJ;
		}else{
			if(meanDistFromOne < okErr){
				setTypes[set] = TPAT_EigValSet_Tp::EIGSET_ONES;
			}else{
				setTypes[set] = TPAT_EigValSet_Tp::EIGSET_REAL_RECIP;
			}
		}
	}

	// Find bifurcations
	std::vector<int> bifs;
	for(size_t m = 1; m < members.size(); m++){
		for(int set = 0; set < 3; set++){
			std::vector<cdouble> eigs = members[m].getEigVals();
			std::vector<cdouble> prevEigs = members[m-1].getEigVals();

			switch(setTypes[set]){
				case TPAT_EigValSet_Tp::EIGSET_REAL_RECIP:
				{
					// Compute distance of each eigenvalues from +/- 1
					double d1 = 1 - std::abs(std::real(eigs[set*2]));
					double d2 = 1 - std::abs(std::real(eigs[set*2+1]));
					double prev_d1 = 1 - std::abs(std::real(prevEigs[set*2]));
					double prev_d2 = 1 - std::abs(std::real(prevEigs[set*2+1]));

					// Check to make sure magnitude of differences is significant
					if( std::abs(d1) > okErr && std::abs(d2) > okErr &&
						std::abs(prev_d1) > okErr && std::abs(prev_d2) > okErr){
						if(d1*prev_d1 < 0 && d2*prev_d2 < 0){
							// Sign changed; bifurcation!
							// printf("Located a bifurcation!\n");
							// printf("  Member %03zu - %03zu\n", m-1, m);
							bifs.push_back(m-1);
						}
					}
					break;
				}
				case TPAT_EigValSet_Tp::EIGSET_COMP_CONJ:
				{
					double meanPrevImag = (std::abs(std::imag(prevEigs[set*2])) + std::abs(std::imag(prevEigs[set*2+1])))/2;
					double meanImag = (std::abs(std::imag(eigs[set*2])) + std::abs(std::imag(eigs[set*2+1])))/2;

					// Check to see if we moved from complex to real or vice versa
					if( (meanPrevImag > okErr) != (meanImag > okErr) ){
						// printf("Located a bifurcation!\n");
						// printf("  Member %03zu - %03zu\n", m-1, m);
						bifs.push_back(m-1);
					}
				}
				default: break;
			}
		}
	}

	return bifs;
}//====================================================

/**
 *	@brief Attempt to load data from the specified matlab file
 *	@param matFile a pointer to the opened matlab file
 *	@throws TPAT_Exception if the variable cannot be loaded
 */	
void TPAT_Fam_CR3BP::loadMemberData(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, DATA_VAR_NAME);
	if(matvar == NULL){
		throw TPAT_Exception("Could not read member data into family");
	}else{
		int numMembers = matvar->dims[0];
		if(matvar->dims[1] != DATA_WIDTH){
			throw TPAT_Exception("Incompatible data file: Data widths are different.");
		}

		if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(matvar->data);

			if(data != NULL){
				for(int i = 0; i < numMembers; i++){
					double state[] = {0,0,0,0,0,0};
					state[0] = data[0*numMembers + i];
					state[1] = data[1*numMembers + i];
					state[2] = data[2*numMembers + i];
					state[3] = data[3*numMembers + i];
					state[4] = data[4*numMembers + i];
					state[5] = data[5*numMembers + i];
					TPAT_FamMember_CR3BP mem(state,
						data[6*numMembers + i], data[7*numMembers + i],
						data[8*numMembers + i], data[9*numMembers + i],
						data[10*numMembers + i]);

					members.push_back(mem);
				}
			}
		}else{
			throw TPAT_Exception("TPAT_Fam_CR3BP::loadMemberData: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//==============================================

/**
 *	@brief Load eigenvalues from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the eigenvalues
 *	@param matFile a pointer to the data file in question
 *	@throws TPAT_Exception if the variable cannot be loaded
 */
void TPAT_Fam_CR3BP::loadEigVals(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, EIG_VAR_NAME);
	if(matvar == NULL){
		throw TPAT_Exception("Could not read eigenvalues into family");
	}else{
		int numMembers = matvar->dims[0];
		if(matvar->dims[1] != 6){
			throw TPAT_Exception("Incompatible data file: Data widths are different.");
		}

		if((int)(members.size()) != numMembers){
			throw TPAT_Exception("TPAT_Fam_CR3BP::loadEigVals: # eigenvalues is not same as number of members");
		}

		if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = static_cast<mat_complex_split_t *>(matvar->data);

			if(splitVals != NULL){
				// splitVals holds two void pointers to the real and imaginary parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				// Read data from column-major order matrix, store in row-major order vector
				for(int i = 0; i < numMembers; i++){
					std::vector<cdouble> vals;
					for(int j = 0; j < 6; j++){
						cdouble temp(realParts[j*numMembers + i], imagParts[j*numMembers + i]);
						vals.push_back(temp);
					}

					members[i].setEigVals(vals);
				}
			}
		}else{
			throw TPAT_Exception("TPAT_Fam_CR3BP::loadEigVals: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//=============================================

/**
 *	@brief Sort all members' eigenvalues so they are in the same order.
 *
 *	This is necessary before bifurcations can be accurately located.
 */
void TPAT_Fam_CR3BP::sortEigs(){
	std::vector<cdouble> allEigs;
	// Create a vector with all the eigenvalues
	for(size_t m = 0; m < members.size(); m++){
		std::vector<cdouble> vals = members[m].getEigVals();
		allEigs.insert(allEigs.end(), vals.begin(), vals.end());
	}

	// Sort eigenvalues
	std::vector<int> sortedIxs;
	std::vector<cdouble> sortedEigs = sortEig(allEigs, &sortedIxs);

	// Update all eigenvalues in the family members
	for(size_t m = 0; m < members.size(); m++){
		std::vector<cdouble> sortedVals(sortedEigs.begin() + m*6, sortedEigs.begin() + (m+1)*6);
		members[m].setEigVals(sortedVals);
	}
}//=================================================

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
 *	
 *	@throws TPAT_Exception if the sorting type is not recognized
 */
void TPAT_Fam_CR3BP::sortMembers(){
	// Don't do any sorting if the sort type is NONE
	if(sortType == TPAT_SortFam_Tp::SORT_NONE)
		return;

	// Create an array containing the independent variable from each family member
	std::vector<double> dataToSort;
	switch(sortType){
		case TPAT_SortFam_Tp::SORT_X:
			getCoord(0, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_Y:
			getCoord(1, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_Z:
			getCoord(2, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_VX:
			getCoord(3, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_VY:
			getCoord(4, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_VZ:
			getCoord(5, &dataToSort); break;
		case TPAT_SortFam_Tp::SORT_JC:
			for(int i = 0; i < ((int)members.size()); i++){
				dataToSort.push_back(members[i].getJacobi());
			}
			break;
		case TPAT_SortFam_Tp::SORT_TOF:
			for(int i = 0; i < ((int)members.size()); i++){
				dataToSort.push_back(members[i].getTOF());
			}
			break;
		default:
			throw TPAT_Exception("Unrecognized sorting type");
			break;
	}

	// Sort the data, retrieve the indices of the now-sorted elements
	std::vector<int> indices = TPAT_Util::getSortedInd(dataToSort);

	// Use those indices to sort the family members
	std::vector<TPAT_FamMember_CR3BP> sortedMembers;
	for(int n = 0; n < ((int)indices.size()); n++){
		sortedMembers.push_back(members[indices[n]]);
	}

	// Reassign the members array
	members = sortedMembers;
}//===================================================

/**
 *	@brief Save the family to a mat-file
 *	@param filename a path the the mat-file in question
 */
void TPAT_Fam_CR3BP::saveToMat(const char *filename){
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
		printErr("TPAT_Fam_CR3BP::saveToMat: Error creating MAT file\n");
	}else{
		// save things
		saveMembers(matfp);
		saveMiscData(matfp);
		saveEigVals(matfp);
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Save ICs, TOFs, and JCs for each member in a matrix
 *	@param matFile a pointer to the destination matlab file
 */
void TPAT_Fam_CR3BP::saveMembers(mat_t *matFile){
	if(members.size() > 0){
		// Create a vector with member data stored in column-major order
		std::vector<double> allData(members.size()*DATA_WIDTH);
		for(size_t i = 0; i < members.size(); i++){
			for(size_t j = 0; j < DATA_WIDTH; j++){
				if(j < 6)
					allData[j*members.size() + i] = members[i].getIC()[j];
				else if(j == 6)
					allData[j*members.size() + i] = members[i].getTOF();
				else if(j == 7)
					allData[j*members.size() + i] = members[i].getJacobi();
				else if(j == 8)
					allData[j*members.size() + i] = members[i].getXAmplitude();
				else if(j == 9)
					allData[j*members.size() + i] = members[i].getYAmplitude();
				else
					allData[j*members.size() + i] = members[i].getZAmplitude();
			}
		}

		size_t dims[2] = {members.size(), static_cast<size_t>(DATA_WIDTH)};
		matvar_t *matvar = Mat_VarCreate(DATA_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(allData[0]), MAT_F_DONT_COPY_DATA);
		saveVar(matFile, matvar, DATA_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	@brief Save eigenvalue data to a mat file
 *	@param matFile a pointer to the mat file in question
 *	@throws TPAT_Exception if a family member does not have six eigenvalues
 */
void TPAT_Fam_CR3BP::saveEigVals(mat_t *matFile){
	if(members.size() > 0){
		// Separate all eigenvalues into real and complex parts
		std::vector<double> realParts(members.size()*6);
		std::vector<double> imagParts(members.size()*6);
		for(size_t i = 0; i < members.size(); i++){
			std::vector<cdouble> vals = members[i].getEigVals();
				if(vals.size() != 6)
					throw TPAT_Exception("TPAT_Fam_CR3BP::saveEigVals: family member does not have 6 eigenvalues!");
			for(int j = 0; j < 6; j++){
				realParts[j*members.size() + i] = std::real(vals[j]);
				imagParts[j*members.size() + i] = std::imag(vals[j]);
			}
		}

		// create a special variable for them
		mat_complex_split_t splitVals = {&(realParts[0]), &(imagParts[0])};

		size_t dims[2] = {members.size(), 6};
		matvar_t *matvar = Mat_VarCreate(EIG_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &splitVals, MAT_F_COMPLEX);
		saveVar(matFile, matvar, EIG_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	@brief Save other useful information to a matlab file
 *	@param matFile the destination matlab file
 */
void TPAT_Fam_CR3BP::saveMiscData(mat_t *matFile){
	// sortType
	int type = static_cast<int>(sortType);
	size_t dims[2] = {1,1};
	matvar_t *typeVar = Mat_VarCreate(SORT_TYPE_VAR_NAME, MAT_C_INT32, MAT_T_INT32, 2, dims, &type, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, typeVar, SORT_TYPE_VAR_NAME, MAT_COMPRESSION_NONE);

	// name
	char name_str[128];
	std::strcpy(name_str, name.c_str());
	dims[1] = name.length();
	matvar_t *nameVar = Mat_VarCreate(NAME_VAR_NAME, MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(name_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, nameVar, NAME_VAR_NAME, MAT_COMPRESSION_NONE);
}//====================================================

