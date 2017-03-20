/**
 *	\file Fam_cr3bp.cpp
 *	\brief Data object for a CR3BP family
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
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
 
#include "Fam_cr3bp.hpp"

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "Common.hpp"
#include "Constraint.hpp"
#include "CorrectionEngine.hpp"
#include "MultShootData.hpp"
#include "Nodeset_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <cstring>


namespace astrohelion{
//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	\brief Create an empty family for the specified system
 *	\param data a CR3BP system data object
 */
Fam_cr3bp::Fam_cr3bp(SysData_cr3bp data){
	sysData = data;
}//====================================================

/**
 *	\brief load a family from a file
 *	\param filepath an aboslute or relative filepath to the data file
 *	\throws Exception if the file cannot be opened
 */
Fam_cr3bp::Fam_cr3bp(const char* filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Fam_cr3bp: Could not load family from file");
	}

	loadMemberData(matfp);
	loadEigVals(matfp);
	loadEigVecs(matfp);
	loadSTMs(matfp);
	name = astrohelion::readStringFromMat(matfp, NAME_VAR_NAME, MAT_T_UINT8, MAT_C_CHAR);
	int type = static_cast<int>(astrohelion::readDoubleFromMat(matfp, SORT_TYPE_VAR_NAME));
	sortType = static_cast<FamSort_tp>(type);
	// sysData.readFromMat(matfp);

	Mat_Close(matfp);
	sysData = SysData_cr3bp(filepath);
}//====================================================

/**
 *	\brief Copy Constructor
 *	\brief fam a family reference
 */
Fam_cr3bp::Fam_cr3bp(const Fam_cr3bp& fam){
	copyMe(fam);
}//====================================================

/**
 *	\brief Destructor
 */
Fam_cr3bp::~Fam_cr3bp(){}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	\brief Assignment operator
 *	\param fam a family reference
 */
Fam_cr3bp& Fam_cr3bp::operator= (const Fam_cr3bp& fam){
	copyMe(fam);
	return *this;
}//====================================================


//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	\brief Add a member to the family
 *	\param mem a new family member; NOTE: <tt>mem</tt> should represent
 *	a trajectory that exists in the same system as the other family members.
 */
void Fam_cr3bp::addMember(FamMember_cr3bp mem){
	members.push_back(mem);
}//====================================================

/**
 *	\brief Retrieve a family member by its index
 *	\param ix the index of the member. If the index is less than zero,
 *	it will count backwards from the end.
 */
FamMember_cr3bp Fam_cr3bp::getMember(int ix) const{
	if(ix < 0)
		ix += members.size();

	return members.at(ix);
}//====================================================

/**
 *	\brief Locate all members with the specified value of Jacobi Constant
 *	\param jc the desired value for Jacobi
 *	\return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<FamMember_cr3bp> Fam_cr3bp::getMemberByJacobi(double jc) const{
	// Get an array of all the jacobi values
	std::vector<double> allJC;
	for(unsigned int n = 0; n < members.size(); n++){
		allJC.push_back(members[n].getJacobi());
	}

	Constraint jacobiCon(Constraint_tp::JC, 0, &jc, 1);

	return getMatchingMember(jc, &allJC, jacobiCon);
}//==============================================

/**
 *	\brief Locate all members with the specified time-of-flight
 *	\param tof the desired value for time-of-flight
 *	\return a vector of matching family members (some may be interpolated/corrected)
 */
std::vector<FamMember_cr3bp> Fam_cr3bp::getMemberByTOF(double tof) const{
	std::vector<double> allTOF;
	for(unsigned int n = 0; n < members.size(); n++){
		allTOF.push_back(members[n].getTOF());
	}

	Constraint tofCon(Constraint_tp::TOF, 0, &tof, 1);

	return getMatchingMember(tof, &allTOF, tofCon);
}//==============================================

/**
 *	\brief Locate all members with the specified value for the specified state variable in the IC
 *	\param value the desired value
 *	\param ix the index of the state variable [0, 5]
 *	\return a vector of matching family members (some may be interpolated/corrected)
 *	\throws Exception if <tt>ix</tt> is out of range
 */
std::vector<FamMember_cr3bp> Fam_cr3bp::getMemberByStateVar(double value, int ix) const{
	if(ix < 0 || ix > 5){
		throw Exception("Fam_cr3bp::getMemberByStateVar: Invalid state index; out of range");
	}

	std::vector<double> allVals;
	for(unsigned int n = 0; n < members.size(); n++){
		allVals.push_back(members[n].getIC()[ix]);
	}

	double conData[] = {NAN,NAN,NAN,NAN,NAN,NAN};
	conData[ix] = value;
	Constraint stateCon(Constraint_tp::STATE, 0, conData, 6);

	return getMatchingMember(value, &allVals, stateCon);
}//==============================================

/**
 *	\brief Retrieve the name of this family
 *	\return a descriptive name
 */
std::string Fam_cr3bp::getName() const { return name; }

/**
 *	\brief Retrieve number of members in this family
 *	\return number of members in this family
 */
int Fam_cr3bp::getNumMembers() const { return members.size(); }

/**
 *	\brief Determine which variable best naturally describes the flow of the family
 *	\return the sorting variable
 */
FamSort_tp Fam_cr3bp::getSortType() const { return sortType; }

/**
 *	\brief Retrieve a string describing the sort type in human-readable format
 *	\return a string describing the sort type in human-readable format
 */
const char* Fam_cr3bp::getSortTypeStr() const{
	switch(sortType){
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
}//===========================================

/**
 *	\brief Retrieve the system data for this family
 *	\return the system data describing the system all family members exist in
 */
SysData_cr3bp Fam_cr3bp::getSysData() const { return sysData; }

/**
 *	\brief Retrieve a pointer to the system data object for this family
 *	\return a pointer to the system data object for this family
 */
SysData_cr3bp* Fam_cr3bp::getSysDataPtr() { return &sysData; }

/**
 *	\brief Set the name of this family
 *	\param n a descriptive name; must be less than 128 characters
 *	if you want to save to a matlab file
 */
void Fam_cr3bp::setName(std::string n){ name = n; }

/**
 *	\brief Set the sort type for this family
 *	\param type the sort type
 */
void Fam_cr3bp::setSortType(FamSort_tp type){ sortType = type; }
//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	\brief Copy the family from another family
 *	\param fam a different family
 */
void Fam_cr3bp::copyMe(const Fam_cr3bp& fam){
	members = fam.members;
	sysData = fam.sysData;
	name = fam.name;
	sortType = fam.sortType;
}//====================================================

/**
 *	\brief Locate places in a data set where a specific value probably exists
 *
 *	This algorithm will locate both exact matches (within a tolerance) and intersections,
 *	assuming the data is continuous. If an intersection is found, the index of the point
 *	before the intersection is returned.
 *
 *	\param value the value to search for
 *	\param data a pointer to a data set to search in
 *	\return a vector of integers representing the indices of matches
 */
std::vector<int> Fam_cr3bp::findMatches(double value, std::vector<double> *data) const{
	double numBins = data->size() > 500 ? 100 : (data->size()/5.0);
	int binSize = floor((data->size())/numBins);

	std::vector<int> matches;
	double minDiff = 0;
	double diff = 0;
	int minDiffIx = 0;
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

/**
 *	\brief Locate a family member with a specific attribute
 *	
 *	This function locates a family member or set of members that have a specific value
 *	for one of the variables of interest (e.g. coordinates, Jacobi, TOF). Exact matches
 *	and interpolated matches are returned; interpolated matches are computed using a
 *	differential corrections algorithm.
 *
 *	\param value the value the family member should have
 *	\param dataSet a pointer to a vector containing the set of values to search for matches
 *	in. For example, if the <tt>value</tt> I pass in contains a specific TOF, then 
 *	<tt>dataSet</tt> points to a vector containing the TOFs for the entire family, sorted
 *	according to this family's <tt>sortType</tt>.
 *	\param matchCon a constraint that can be applied in a corrections scheme that will
 *	ensure the corrected trajectory has the desired value for the variable of interest.
 *
 *	\return a vector of matches. If no matches are returned, the vector will be empty.
 */
std::vector<FamMember_cr3bp> Fam_cr3bp::getMatchingMember(double value, std::vector<double> *dataSet,
	Constraint matchCon) const{
	// Locate possible candidates
	std::vector<int> matches = findMatches(value, dataSet);
	std::vector<FamMember_cr3bp> matchMembers;
	if(matches.size() == 0){
		astrohelion::printErr("Could not locate any matches. The family either has too few members to facilitate an accurate search or the desired trajectory does not exist.\n");
		return matchMembers;	// empty set
	}else{
		astrohelion::printColor(GREEN, "Located %zu matches; applying corrections\n", matches.size());
	}

	SysData_cr3bp tempSys = sysData;
	
	for(unsigned int n = 0; n < matches.size(); n++){
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

			Nodeset_cr3bp memberSet(&tempSys, members[idx].getIC(), tof, numNodes);
			double end = numNodes-1;
			double conData[] = {end,end,end,end,end,end};
			Constraint periodicCon(Constraint_tp::MATCH_CUST, 0, conData, 6);

			memberSet.addConstraint(periodicCon);
			memberSet.addConstraint(matchCon);

			// Correct the nodeset while constraining the orbit to have the desired characterstic
			CorrectionEngine corrector;
			corrector.setTol(1e-11);
			try{
				Nodeset_cr3bp newNodes(&tempSys);
				corrector.multShoot(&memberSet, &newNodes);
				
				Traj_cr3bp newTraj = Traj_cr3bp::fromNodeset(newNodes);
				FamMember_cr3bp newMember(newTraj);
				matchMembers.push_back(newMember);
			}catch(DivergeException &e){}
		}
	}

	return matchMembers;
}//==============================================

/**
 *	\brief Populate an array with a single coordinate from each family member
 *	\param ix the index of the coordinate within the IC vector
 *	\param array a pointer to the array we want to populate
 */
void Fam_cr3bp::getCoord(int ix, std::vector<double> *array) const{
	for(unsigned int i = 0; i < members.size(); i++){
		std::vector<double> ic = members[i].getIC();
		array->push_back(ic[ix]);
	}
}//====================================================

/**
 *	\brief Locate all bifurcations in the family by analyzing the 
 *	eigenvalues
 *
 *	Eigenvalues MUST be sorted, or this will yield completely bogus
 *	results
 */
std::vector<int> Fam_cr3bp::findBifurcations(){
	double okErr = 1e-6;
	std::vector<int> bifs;
	double prevStab[3] = {0};
	for(unsigned int m = 0; m < members.size(); m++){
		double stab[3] = {0};
		std::vector<cdouble> eigs = members[m].getEigVals();
		for(unsigned int i = 0; i < 3; i++){
			stab[i] = 1.0 - std::abs(0.5*(eigs[2*i] + eigs[2*i+1]));
		}

		if(m == 0)
			std::copy(stab, stab+3, prevStab);
		else{
			// Look for changes in sign (stability index crosses 1)
			for(unsigned int i = 0; i < 3; i++){
				if(stab[i]*prevStab[i] < 0 && std::abs(stab[i] - prevStab[i]) > okErr){
					// Bifurcation!
					bifs.push_back(m-1);
					break;
				}
			}

			std::copy(stab, stab+3, prevStab);
		}
	}

	// // TODO: Incorporate much more advanced ways to compute this
	// double okErr = 1e-3;
	// cdouble one(1,0);

	// EigValSet_tp setTypes[3];

	// // Determine what "type" of eigenvalue pair we have
	// for(int set = 0; set < 3; set++){
	// 	double sumImag = 0;
	// 	// double sumReal = 0;
	// 	double sumDistFromOne = 0;

	// 	for(unsigned int m = 0; m < members.size(); m++){
	// 		std::vector<cdouble> eigs = members[m].getEigVals();
	// 		sumImag += (std::abs(std::imag(eigs[set*2])) + std::abs(std::imag(eigs[set*2+1])))/2;
	// 		sumDistFromOne += (std::abs(eigs[set*2] - one) + std::abs(eigs[set*2+1] - one))/2;
	// 	}

	// 	double count = members.size();
	// 	double meanImag = sumImag/count;
	// 	// double meanReal = sumReal/count;
	// 	double meanDistFromOne = sumDistFromOne/count;

	// 	if(meanImag > okErr){
	// 		// significant imaginary parts
	// 		setTypes[set] = EigValSet_tp::EIGSET_COMP_CONJ;
	// 	}else{
	// 		if(meanDistFromOne < okErr){
	// 			setTypes[set] = EigValSet_tp::EIGSET_ONES;
	// 		}else{
	// 			setTypes[set] = EigValSet_tp::EIGSET_REAL_RECIP;
	// 		}
	// 	}
	// }

	// // Find bifurcations
	// std::vector<int> bifs;
	// for(unsigned int m = 1; m < members.size(); m++){
	// 	for(int set = 0; set < 3; set++){
	// 		std::vector<cdouble> eigs = members[m].getEigVals();
	// 		std::vector<cdouble> prevEigs = members[m-1].getEigVals();

	// 		switch(setTypes[set]){
	// 			case EigValSet_tp::EIGSET_REAL_RECIP:
	// 			{
	// 				// Compute distance of each eigenvalues from +/- 1
	// 				double d1 = 1 - std::abs(std::real(eigs[set*2]));
	// 				double d2 = 1 - std::abs(std::real(eigs[set*2+1]));
	// 				double prev_d1 = 1 - std::abs(std::real(prevEigs[set*2]));
	// 				double prev_d2 = 1 - std::abs(std::real(prevEigs[set*2+1]));

	// 				// Check to make sure magnitude of differences is significant
	// 				if( std::abs(d1) > okErr && std::abs(d2) > okErr &&
	// 					std::abs(prev_d1) > okErr && std::abs(prev_d2) > okErr){
	// 					if(d1*prev_d1 < 0 && d2*prev_d2 < 0){
	// 						// Sign changed; bifurcation!
	// 						// printf("Located a bifurcation!\n");
	// 						// printf("  Member %03zu - %03zu\n", m-1, m);
	// 						bifs.push_back(m-1);
	// 					}
	// 				}
	// 				break;
	// 			}
	// 			case EigValSet_tp::EIGSET_COMP_CONJ:
	// 			{
	// 				double meanPrevImag = (std::abs(std::imag(prevEigs[set*2])) + std::abs(std::imag(prevEigs[set*2+1])))/2;
	// 				double meanImag = (std::abs(std::imag(eigs[set*2])) + std::abs(std::imag(eigs[set*2+1])))/2;

	// 				// Check to see if we moved from complex to real or vice versa
	// 				if( (meanPrevImag > okErr) != (meanImag > okErr) ){
	// 					// printf("Located a bifurcation!\n");
	// 					// printf("  Member %03zu - %03zu\n", m-1, m);
	// 					bifs.push_back(m-1);
	// 				}
	// 			}
	// 			default: break;
	// 		}
	// 	}
	// }

	return bifs;
}//====================================================

/**
 *	\brief Attempt to load data from the specified matlab file
 *	\param matFile a pointer to the opened matlab file
 *	\throws Exception if the variable cannot be loaded
 */	
void Fam_cr3bp::loadMemberData(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, DATA_VAR_NAME);
	if(matvar == NULL){
		throw Exception("Could not read member data into family");
	}else{
		int numMembers = matvar->dims[0];
		if(matvar->dims[1] != DATA_WIDTH){
			throw Exception("Incompatible data file: Data widths are different.");
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
					FamMember_cr3bp mem(state,
						data[6*numMembers + i], data[7*numMembers + i],
						data[8*numMembers + i], data[9*numMembers + i],
						data[10*numMembers + i]);

					members.push_back(mem);
				}
			}
		}else{
			throw Exception("Fam_cr3bp::loadMemberData: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//==============================================

/**
 *	\brief Load eigenvalues from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the eigenvalues
 *	\param matFile a pointer to the data file in question
 *	\throws Exception if the variable cannot be loaded
 */
void Fam_cr3bp::loadEigVals(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, EIG_VAR_NAME);
	if(matvar == NULL){
		throw Exception("Could not read eigenvalues into family");
	}else{
		unsigned int numMembers = matvar->dims[0];
		if(matvar->dims[1] != 6){
			throw Exception("Incompatible data file: Data widths are different.");
		}

		if(members.size() != numMembers){
			throw Exception("Fam_cr3bp::loadEigVals: # eigenvalues is not same as number of members");
		}

		if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = static_cast<mat_complex_split_t *>(matvar->data);

			if(splitVals != NULL){
				// splitVals holds two void pointers to the real and imaginary parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				// Read data from column-major order matrix, store in row-major order vector
				for(unsigned int i = 0; i < numMembers; i++){
					std::vector<cdouble> vals;
					for(int j = 0; j < 6; j++){
						cdouble temp(realParts[j*numMembers + i], imagParts[j*numMembers + i]);
						vals.push_back(temp);
					}

					members[i].setEigVals(vals);
				}
			}
		}else{
			throw Exception("Fam_cr3bp::loadEigVals: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//=============================================

/**
 *	\brief Load eigenvectors from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the eigenvectors
 *	\param pMatFile a pointer to the data file in question
 *	\throws Exception if the variable cannot be loaded
 */
void Fam_cr3bp::loadEigVecs(mat_t* pMatFile){
	matvar_t *pMatvar = Mat_VarRead(pMatFile, EIGVEC_VAR_NAME);
	if(pMatvar == NULL){
		throw Exception("Fam_cr3bp::loadEigVecs: Could not read data vector");
	}else{
		unsigned int numSteps = pMatvar->dims[2];

		if(members.size() == 0){
			throw Exception("Fam_cr3bp::loadEigVecs: Member vector has not been initialized!");
		}

		if(numSteps != members.size() ){
			throw Exception("Fam_cr3bp::loadEigVecs: Eigenvector vector does not have the same number of elements as the member vector");
		}

		if(pMatvar->dims[0] != 6 || pMatvar->dims[1] != 6){
			throw Exception("Fam_cr3bp::loadEigVecs: Incompatible data file: Eigenvector matrix is not 6x6.");
		}

		if(pMatvar->class_type == MAT_C_DOUBLE && pMatvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = static_cast<mat_complex_split_t *>(pMatvar->data);

			if(splitVals != NULL){
				// splitVals holds two void pointers to the real and imaginary parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				for(unsigned int i = 0; i < numSteps; i++){
					cdouble vecData[36];
					for(unsigned int j = 0; j < 36; j++){
						vecData[j] = cdouble(realParts[j*numSteps + i], imagParts[j*numSteps + i]);
					}
					MatrixXRcd vecMat = Eigen::Map<MatrixXRcd>(vecData, 6, 6);
					members[i].setEigVecs(vecMat.transpose());
				}
			}
		}else{
			throw Exception("Fam_cr3bp::loadEigVecs: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(pMatvar);
}//=============================================

/**
 *	\brief Load state transition matrices from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the STMs
 *	\param pMatFile a pointer to the data file in question
 *	\throws Exception if the variable cannot be loaded
 */
void Fam_cr3bp::loadSTMs(mat_t* pMatFile){
	matvar_t *pAllSTM = Mat_VarRead(pMatFile, VARNAME_STM);
	if(pAllSTM == NULL){
		throw Exception("Fam_cr3bp::loadSTMs: Could not read data vector");
	}else{
		unsigned int numSteps = pAllSTM->dims[2];

		if(members.size() == 0){
			throw Exception("Fam_cr3bp::loadSTMs: Member vector has not been initialized!");
		}

		if(numSteps != members.size() ){
			throw Exception("Fam_cr3bp::loadSTMs: STM vector does not have the same number of elements as the member vector");
		}

		if(pAllSTM->dims[0] != 6 || pAllSTM->dims[1] != 6){
			throw Exception("Fam_cr3bp::loadSTMs: Incompatible data file: STM is not 6x6.");
		}

		if(pAllSTM->class_type == MAT_C_DOUBLE && pAllSTM->data_type == MAT_T_DOUBLE){
			double *data = static_cast<double *>(pAllSTM->data);

			if(data != NULL){
				for(unsigned int i = 0; i < numSteps; i++){
					double stmEl[36];
					for(unsigned int j = 0; j < 36; j++){
						stmEl[j] = data[36*i + j];
					}

					MatrixXRd P = Eigen::Map<MatrixXRd>(stmEl, 6, 6);
					members[i].setSTM(P.transpose());
				}
			}
		}else{
			throw Exception("Fam_cr3bp::loadSTMs: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(pAllSTM);
}//=============================================

/**
 *	\brief Sort all members' eigenvalues so they are in the same order.
 *	\details This is necessary before bifurcations can be accurately located.
 */
void Fam_cr3bp::sortEigs(){
	std::vector<cdouble> allVals;
	std::vector<MatrixXRcd> allVecs;
	// Create a vector with all the eigenvalues
	for(unsigned int m = 0; m < members.size(); m++){
		std::vector<cdouble> vals = members[m].getEigVals();
		allVals.insert(allVals.end(), vals.begin(), vals.end());
		allVecs.insert(allVecs.end(), members[m].getEigVecs());
	}

	// Sort eigenvalues
	std::vector<unsigned int> sortedIxs = sortEig(allVals, allVecs);
	// std::vector<cdouble> sortedEigs = sortEig(allVals, &sortedIxs);

	// Update all eigenvalues in the family members
	std::vector<cdouble> sortedVals;
	MatrixXRcd sortedVecs = MatrixXRcd::Zero(6,6);
	for(unsigned int m = 0; m < members.size(); m++){
		// Construct new eigenvalue and eigenvector objects with the new index order
		for(unsigned int c = 0; c < 6; c++){
			unsigned int ix = sortedIxs[m*6+c];
			sortedVals.push_back(allVals[m*6+ix]);

			sortedVecs.col(c) = allVecs[m].col(ix);
		}
		
		// Update member with sorted vectors and values
		members[m].setEigVals(sortedVals);
		members[m].setEigVecs(sortedVecs);

		// Reset storage variables
		sortedVals.clear();
		sortedVecs = MatrixXRcd::Zero(6,6);
	}
}//=================================================

/**
 *	\brief Sort the family members by the specified sort variable (in ascending order)
 *	\details The sorting variable is specified by <tt>sortType</tt>; this is the variable
 *	that best describes the natural progression of the family. For example,
 *	Lyapunov orbits can be evolved naturally by varying the x-coordinate of the IC.
 *
 *	The family must be sorted (or be loaded from an already sorted family file) before
 *	you can retrieve family members. The process will run without sorting, but the results
 *	will likely be wonky.
 *	
 *	\throws Exception if the sorting type is not recognized
 */
void Fam_cr3bp::sortMembers(){
	// Don't do any sorting if the sort type is NONE
	if(sortType == FamSort_tp::SORT_NONE)
		return;

	// Create an array containing the independent variable from each family member
	std::vector<double> dataToSort;
	switch(sortType){
		case FamSort_tp::SORT_X:
			getCoord(0, &dataToSort); break;
		case FamSort_tp::SORT_Y:
			getCoord(1, &dataToSort); break;
		case FamSort_tp::SORT_Z:
			getCoord(2, &dataToSort); break;
		case FamSort_tp::SORT_VX:
			getCoord(3, &dataToSort); break;
		case FamSort_tp::SORT_VY:
			getCoord(4, &dataToSort); break;
		case FamSort_tp::SORT_VZ:
			getCoord(5, &dataToSort); break;
		case FamSort_tp::SORT_JC:
			for(unsigned int i = 0; i < members.size(); i++){
				dataToSort.push_back(members[i].getJacobi());
			}
			break;
		case FamSort_tp::SORT_TOF:
			for(unsigned int i = 0; i < members.size(); i++){
				dataToSort.push_back(members[i].getTOF());
			}
			break;
		default:
			throw Exception("Unrecognized sorting type");
			break;
	}

	// Sort the data, retrieve the indices of the now-sorted elements
	std::vector<int> indices = astrohelion::getSortedInd(dataToSort);

	// Use those indices to sort the family members
	std::vector<FamMember_cr3bp> sortedMembers;
	for(unsigned int n = 0; n < indices.size(); n++){
		sortedMembers.push_back(members[indices[n]]);
	}

	// Reassign the members array
	members = sortedMembers;
}//===================================================

/**
 *	\brief Save the family to a mat-file
 *	\param filename a path the the mat-file in question
 */
void Fam_cr3bp::saveToMat(const char *filename){
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
		astrohelion::printErr("Fam_cr3bp::saveToMat: Error creating MAT file\n");
	}else{
		// save things
		saveMembers(matfp);
		saveMiscData(matfp);
		saveEigVals(matfp);
		saveEigVecs(matfp);
		saveSTMs(matfp);
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	\brief Save ICs, TOFs, and JCs for each member in a matrix
 *	\param pMatFile a pointer to the destination matlab file
 */
void Fam_cr3bp::saveMembers(mat_t *pMatFile){
	if(members.size() > 0){
		// Create a vector with member data stored in column-major order
		std::vector<double> allData(members.size()*DATA_WIDTH);
		for(unsigned int i = 0; i < members.size(); i++){
			for(unsigned int j = 0; j < DATA_WIDTH; j++){
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
		astrohelion::saveVar(pMatFile, matvar, DATA_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save eigenvalue data to a mat file
 *	\param pMatFile a pointer to the mat file in question
 *	\throws Exception if a family member does not have six eigenvalues
 */
void Fam_cr3bp::saveEigVals(mat_t *pMatFile){
	if(members.size() > 0){
		// Separate all eigenvalues into real and complex parts
		std::vector<double> realParts(members.size()*6);
		std::vector<double> imagParts(members.size()*6);
		for(unsigned int i = 0; i < members.size(); i++){
			std::vector<cdouble> vals = members[i].getEigVals();
				if(vals.size() != 6)
					throw Exception("Fam_cr3bp::saveEigVals: family member does not have 6 eigenvalues!");
			for(unsigned int j = 0; j < 6; j++){
				realParts[j*members.size() + i] = std::real(vals[j]);
				imagParts[j*members.size() + i] = std::imag(vals[j]);
			}
		}

		// create a special variable for them
		mat_complex_split_t splitVals = {&(realParts[0]), &(imagParts[0])};

		size_t dims[2] = {members.size(), 6};
		matvar_t *matvar = Mat_VarCreate(EIG_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, matvar, EIG_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save eigenvector data to a mat file
 *	\param pMatFile a pointer to the mat file in question
 */
void Fam_cr3bp::saveEigVecs(mat_t *pMatFile){
	if(members.size() > 0){
		
		std::vector<double> allVec_real(members.size()*36);
		std::vector<double> allVec_imag(members.size()*36);

		for(unsigned int i = 0; i < members.size(); i++){
			// get the transpose of the STM matrix; we need to store it in column-major order
			// and it's currently in row-major order
			MatrixXRcd P = members[i].getEigVecs().transpose();
			
			// Retrieve the data from the matrix
			cdouble *matData = P.data();
			
			// Store that data in our huge vector
			for(unsigned int j = 0; j < 36; j++){
				allVec_real[i*36 + j] = std::real(matData[j]);
				allVec_imag[i*36 + j] = std::imag(matData[j]);
			}
		}
		mat_complex_split_t splitVals = {&(allVec_real[0]), &(allVec_imag[0])};

		size_t dims[3] = {6, 6, members.size()};
		matvar_t *pMatVar = Mat_VarCreate(EIGVEC_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, pMatVar, EIGVEC_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save STM data to a mat file
 *	\param pMatFile a pointer to the mat file in question
 */
void Fam_cr3bp::saveSTMs(mat_t *pMatFile){
	if(members.size() > 0){
		std::vector<double> allSTMEl(members.size()*36);

		for(unsigned int i = 0; i < members.size(); i++){
			// get the transpose of the STM matrix; we need to store it in column-major order
			// and it's currently in row-major order
			MatrixXRd P = members[i].getSTM().transpose();
			// Retrieve the data from the matrix
			double *matData = P.data();
			// Store that data in our huge vector
			std::copy(matData, matData+36, &(allSTMEl[0]) + i*36);
		}

		size_t dims[3] = {6, 6, members.size()};
		matvar_t *pMatVar = Mat_VarCreate(VARNAME_STM, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, &(allSTMEl[0]), MAT_F_DONT_COPY_DATA);
		astrohelion::saveVar(pMatFile, pMatVar, VARNAME_STM, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save other useful information to a matlab file
 *	\param pMatFile the destination matlab file
 */
void Fam_cr3bp::saveMiscData(mat_t *pMatFile){
	// sortType
	int type = static_cast<int>(sortType);
	size_t dims[2] = {1,1};
	matvar_t *typeVar = Mat_VarCreate(SORT_TYPE_VAR_NAME, MAT_C_INT32, MAT_T_INT32, 2, dims, &type, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(pMatFile, typeVar, SORT_TYPE_VAR_NAME, MAT_COMPRESSION_NONE);

	// name
	char name_str[128];
	std::strcpy(name_str, name.c_str());
	dims[1] = name.length();
	matvar_t *nameVar = Mat_VarCreate(NAME_VAR_NAME, MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(name_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(pMatFile, nameVar, NAME_VAR_NAME, MAT_COMPRESSION_NONE);
}//====================================================



}// END of Astrohelion namespace