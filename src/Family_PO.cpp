/**
 *  @file Family_PO.cpp
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

#include "Family_PO.hpp"

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "ControlLaw.hpp"
#include "DynamicsModel.hpp"
#include "Exceptions.hpp"
#include "MultShootEngine.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors and Desctructor
//-----------------------------------------------------------------------------

/**
 * @brief Default constructor
 * @param pSys System data associated with every member of the family
 */
Family_PO::Family_PO(const SysData *pSys) : Family(pSys) {}

/**
 * @brief Construct the family from another object
 * @param f reference to another family object
 */
Family_PO::Family_PO(const Family_PO &f): Family(f) { copyMe(f); }

/**
 * @brief Default destructor; does nothing
 */
Family_PO::~Family_PO() {}

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

/**
 * @brief Copy constructor
 * 
 * @param f Another object to copy
 * @return reference to this object after it is copied from the input, `f`
 */
Family_PO& Family_PO::operator= (const Family_PO &f){
	copyMe(f);
	return *this;
}//====================================================

/**
 *  @brief Concatenates two Family_PO objects
 *  @details Performs the operation `combo = lhs + rhs`. The data from 
 *  `lhs` is copied directly into `combo`, thus, the name, sort type,
 *  and match tolerance of `combo` are identical to those of `lhs`
 *  regardless of the values of these parameters in `rhs`. The members,
 *  eigenvalues, and eigenvectors from `rhs` are appended to those of 
 *  `lhs`.
 * 
 *  @param lhs a Family_PO object reference
 *  @param rhs a Family_PO object reference
 */
Family_PO operator+ (const Family_PO &lhs, const Family_PO &rhs){
	if(*(lhs.pSysData) != *(rhs.pSysData))
		throw Exception("Family_PO, operator+ : lhs and rhs use different "
			"system data structures");

	Family_PO combo(lhs);	// Copy contents of lhs

	// Append contents of rhs to lhs
	combo.members.insert(combo.members.end(), rhs.members.begin(), 
		rhs.members.end());
	combo.memberEigVals.insert(combo.memberEigVals.end(), 
		rhs.memberEigVals.begin(), rhs.memberEigVals.end());
	combo.memberEigVecs.insert(combo.memberEigVecs.end(), 
		rhs.memberEigVecs.begin(), rhs.memberEigVecs.end());

	return combo;
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Add a periodic orbit to the family
 *  @param arc reference to a periodic orbit; a copy is made and
 *  placed in a storage vector
 */
void Family_PO::addMember(const Arcset_periodic &arc){
	// Save the orbit data
	members.push_back(arc);
	
	// Compute the eigenvalue/vector information
	std::vector<cdouble> vals {};
	MatrixXRcd vecs;
	members.back().getEigData(&vals, &vecs);

	// Save the eigenvalue/vector information
	memberEigVals.insert(memberEigVals.end(), vals.begin(), vals.end());
	memberEigVecs.push_back(vecs);
}//====================================================

/**
 * @brief Retrieve the eigenvalues associated with a specific family member
 * 
 * @param ix Index of the family member. If ix < 0, it counts backward from the
 * end of the family.
 * @return a vector of complex eigenvalues associated with the specified
 * family member's monodromy matrix
 * 
 * @throws Exception if the eigenvalue and eigenvector arrays have not been
 * populated
 * @throws Exception if `ix` is out of bounds
 */
std::vector<cdouble> Family_PO::getEigVals(int ix) const{
	if(memberEigVecs.size() != members.size() || memberEigVals.size() == 0)
		throw Exception("Family_PO::getEigVals: eigenvalue/vector vectors"
			" have not been populated");

	while(ix < 0)
		ix += members.size();

	if(ix > static_cast<int>(members.size())){
		char msg[128];
		sprintf(msg, "Family_PO::getEigVals: index %d out of bounds, "
			"max = %zu\n", ix, members.size());
		throw Exception(msg);
	}

	const unsigned int nE = memberEigVecs[0].rows();
	if(memberEigVals.size() != nE*memberEigVecs.size()){
		throw Exception("Family_PO::getEigVals: eigenvalue vector does not "
			"have enough elements to supply eigenvalues for each family "
			"member");
	}

	std::vector<cdouble> vals(memberEigVals.begin() + ix*nE, 
		memberEigVals.begin() + (ix+1)*nE);
	return vals;
}//====================================================

/**
 * @brief Retrieve the eigenvectors associated with a specific family member
 * 
 * @param ix Index of the family member. If ix < 0, it counts backward from the
 * end of the family.
 * @return a matrix of complex eigenvectors associated with the specified
 * family member's monodromy matrix; each column is an eigenvector
 * 
 * @throws Exception if the eigenvector array has not been populated
 * @throws Exception if `ix` is out of bounds
 */
MatrixXRcd Family_PO::getEigVecs(int ix) const{
	if(memberEigVecs.size() != members.size())
		throw Exception("Family_PO::getEigVecs: eigenvector array is not "
			"same size as number of members");

	while(ix < 0)
		ix += members.size();

	if(ix > static_cast<int>(members.size())){
		char msg[128];
		sprintf(msg, "Family_PO::getEigVecs: index %d out of bounds, "
			"max = %zu\n", ix, members.size());
		throw Exception(msg);
	}

	return memberEigVecs[ix];
}//====================================================

/**
 *  @brief Retrieve a copy of a family member
 * 
 *  @param ix index of the family member within the storage vector
 *  @return a copy of a family member
 */
Arcset_periodic Family_PO::getMember(int ix) const{
	while(ix < 0)
		ix += members.size();

	if(ix > static_cast<int>(members.size())){
		char msg[128];
		sprintf(msg, "Family_PO::getMember: index %d out of bounds, "
			"max = %zu\n", ix, members.size());
		throw Exception(msg);
	}

	return members[ix];
}//====================================================

/**
 *  @brief Retrieve a reference to a family member
 * 
 *  @param ix index of the family member within the storage vector
 *  @return a reference to the family member
 */
Arcset_periodic& Family_PO::getMemberRef(int ix){
	while(ix < 0)
		ix += members.size();

	return members[ix];
}//====================================================

/**
 *  @brief Find a family member with a specific value of one of the initial states
 * 
 *  @param val the value of the initial state at the specified index, 
 *  nondimensional units
 *  @param stateIx index of the state variable (e.g., x = 0, y = 1, etc.)
 * 
 *  @return A vector of family members that match the criterion
 */
std::vector<Arcset_periodic> Family_PO::getMemberByState(double val, 
	unsigned int stateIx) const{

	const unsigned int core_dim = pSysData->getDynamicsModel()->getCoreStateSize();
	if(stateIx >= core_dim)
		throw Exception("Family_PO::getMemberByState: State index out of range");

	std::vector<double> allVals;
	getCoord(stateIx, &allVals);

	std::vector<double> conData(core_dim, NAN);
	conData[stateIx] = val;
	Constraint stateCon(Constraint_tp::STATE, 0, conData);

	return getMatchingMember(val, &allVals, stateCon);
}//====================================================

/**
 *  @brief Find a family member with the specified time-of-flight
 * 
 *  @param tof time-of-flight, nondimensional units
 *  @return a vector of family members with matching times-of-flight
 */
std::vector<Arcset_periodic> Family_PO::getMemberByTOF(double tof) const{
	std::vector<double> allTOF(members.size(), NAN);
	for(unsigned int n = 0; n < members.size(); n++){
		allTOF[n] = members[n].getTotalTOF();
	}

	Constraint tofCon(Constraint_tp::TOF_TOTAL, 0, &tof, 1);

	return getMatchingMember(tof, &allTOF, tofCon);
}//====================================================

/**
 *  @brief Retrieve the number of periodic orbits in the family
 *  @return the number of periodic orbits in the family
 */
unsigned int Family_PO::getNumMembers() const { return members.size(); }

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *	@brief Locate all bifurcations in the family by analyzing the 
 *	eigenvalues
 *	@details Eigenvalues MUST be sorted, or this will yield completely bogus
 *	results
 *	
 *	@return the indices of the family members located close to bifurcations
 */
std::vector<unsigned int> Family_PO::findBifurcations(){
	// Find the eigenvalues that are exactly (in theory) equal to one
	double onesErrTotal[3] = {0};
	unsigned int nE = memberEigVecs[0].rows();

	for(unsigned int m = 0; m < members.size(); m++){
		std::vector<cdouble> eigs(memberEigVals.begin() + m*nE, 
			memberEigVals.begin() + (m+1)*nE);
		// printf("Eigs[ %s, %s, %s, %s, %s, %s]\n", complexToStr(eigs[0]).c_str(),
		// 	complexToStr(eigs[1]).c_str(), complexToStr(eigs[2]).c_str(), complexToStr(eigs[3]).c_str(),
		// 	complexToStr(eigs[4]).c_str(), complexToStr(eigs[5]).c_str());

		for(unsigned int e = 0; e < 3; e++){
			onesErrTotal[e] += std::abs(1.0 - eigs[2*e]);
		}
	}
	unsigned int onesIx = 0;
	if(onesErrTotal[0] > onesErrTotal[1] || onesErrTotal[0] > onesErrTotal[1]){
		onesIx = onesErrTotal[1] < onesErrTotal[2] ? 1 : 2;
	}

	// printf("Ones at index %u\n", onesIx);

	double okErr = 1e-6;		// Tolerance to determine if an eigenvalue equals 1.0
	std::vector<unsigned int> bifs;
	double stab[3] = {0};		// Stability indices for current family member
	double prevStab[3] = {0};	// Stability indices from the previous family member
	for(unsigned int m = 0; m < members.size(); m++){
		
		// Compute the stability indices for this family member
		std::vector<cdouble> eigs(memberEigVals.begin() + m*nE, 
			memberEigVals.begin() + (m+1)*nE);
		for(unsigned int i = 0; i < 3; i++){
			stab[i] = 1.0 - std::abs(0.5*(eigs[2*i] + eigs[2*i+1]));

			// if(abs(1.0 - eigs[2*i]) < okErr && abs(1.0 - eigs[2*i+1]) < okErr){
			// 	// Eigenvalues are the two ones
			// 	onesErr = abs(stab[i]);
			// }
		}

		if(m == 0)
			std::copy(stab, stab+3, prevStab);
		else{
			// Look for changes in sign (stability index crosses 1)
			for(unsigned int i = 0; i < 3; i++){
				if(i != onesIx && stab[i]*prevStab[i] < 0 && 
					std::abs(stab[i] - prevStab[i]) > okErr){

					// Bifurcation!
					bifs.push_back(m-1);
					// printf("vvvv -Bifurcation- vvvv\n");
					break;
				}
			}

			std::copy(stab, stab+3, prevStab);
		}

		// printf("Member %03u: stability = [%.4f, %.4f, %.4f]; onesErr = %.2e\n", m,
		// 	stab[0], stab[1], stab[2], onesErr);
	}

	return bifs;
}//====================================================

/**
 *	@brief Sort all members' eigenvalues so they are in the same order.
 *	@details This is necessary before bifurcations can be accurately located.
 */
void Family_PO::sortEigs(){
	if(memberEigVecs.size() == 0)
		return;

	// Sort eigenvalues
	const unsigned int nE = memberEigVecs[0].rows();	// number of eigenvalues
	std::vector<unsigned int> sortedIxs = sortEig(memberEigVals, memberEigVecs);

	// Update all eigendata in the storage vectors
	std::vector<cdouble> sortedVals(nE, 0);
	MatrixXRcd sortedVecs = MatrixXRcd::Zero(nE, nE);
	unsigned int ix = 0, m = 0, c = 0;
	for(m = 0; m < members.size(); m++){
		// Construct new eigenvalue and eigenvector objects with the new index order
		for(c = 0; c < nE; c++){
			ix = sortedIxs[m*nE+c];
			sortedVals[c] = memberEigVals[m*nE+ix];

			sortedVecs.col(c) = memberEigVecs[m].col(ix);
		}
		
		// Update member with sorted vectors and values
		memberEigVecs[m] = sortedVecs;
		for(c = 0; c < nE; c++){
			memberEigVals[m*nE + c] = sortedVals[c];
		}

		// Reset storage variables
		sortedVals.assign(nE, 0);
		sortedVecs = MatrixXRcd::Zero(nE,nE);
	}
}//====================================================

/**
 *	@brief Sort the family members by the specified sort variable (in ascending order)
 *	@details The sorting variable is specified by `sortType`; this is the variable
 *	that best describes the natural progression of the family. For example,
 *	Lyapunov orbits can be evolved naturally by varying the x-coordinate of the IC.
 *
 *	The family must be sorted (or be loaded from an already sorted family file) 
 *	before you can retrieve family members. The process will run without sorting, 
 *	but the results will likely be wonky.
 *	
 *	@throws Exception if the sorting type is not recognized
 */
void Family_PO::sortMembers(){
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
			// for(unsigned int i = 0; i < members.size(); i++){
			// 	dataToSort.push_back(members[i].getJacobi());
			// }
			throw Exception("Family_PO::sortMembers: JC Sorting not "
				"implemented here!");
			break;
		case FamSort_tp::SORT_TOF:
			for(unsigned int i = 0; i < members.size(); i++){
				dataToSort.push_back(members[i].getTotalTOF());
			}
			break;
		default:
			throw Exception("Unrecognized sorting type");
			break;
	}

	// Sort the data, retrieve the indices of the now-sorted elements
	std::vector<int> indices = astrohelion::getSortedInd(dataToSort);

	// Use those indices to sort the family members
	std::vector<Arcset_periodic> sortedMembers;
	std::vector<MatrixXRcd> sortedEigVecs;
	std::vector<cdouble> sortedEigVals;
	unsigned int ix = 0;
	const unsigned int nE = memberEigVecs[0].rows();	// number of eigenvalues
	for(unsigned int n = 0; n < indices.size(); n++){
		ix = indices[n];
		sortedMembers.push_back(members[ix]);
		sortedEigVecs.push_back(memberEigVecs[ix]);
		sortedEigVals.insert(sortedEigVals.end(), memberEigVals.begin() + nE*ix, 
			memberEigVals.begin() + nE*(ix+1));
	}

	// Reassign the members vector
	members = sortedMembers;
	memberEigVecs = sortedEigVecs;
	memberEigVals = sortedEigVals;
}//====================================================

/**
 *	@brief Locate a family member with a specific attribute
 *	
 *	This function locates a family member or set of members that have a specific 
 *	value for one of the variables of interest (e.g. coordinates, Jacobi, TOF). 
 *	Exact matches and interpolated matches are returned; interpolated matches 
 *	are computed using a differential corrections algorithm.
 *
 *	@param val the value the family member should have
 *	@param data a pointer to a vector containing the set of values to search for 
 *	matches in. For example, if the `val` I pass in contains a specific TOF, then 
 *	`data` points to a vector containing the TOFs for the entire family, sorted
 *	according to this family's `sortType`.
 *	@param matchCon a constraint that can be applied in a corrections scheme 
 *	that will ensure the corrected trajectory has the desired value for the 
 *	variable of interest.
 *
 *	@return a vector of matches. If no matches are returned, the vector will be 
 *	empty.
 */
std::vector<Arcset_periodic> Family_PO::getMatchingMember(double val,
	std::vector<double> *data, Constraint matchCon) const{

	// Locate possible candidates
	std::vector<unsigned int> matches = findMatches(val, data);
	std::vector<Arcset_periodic> matchMembers;

	if(matches.size() == 0){
		astrohelion::printErr("Could not locate any matches. The family "
			"either has too few members to facilitate an accurate search or "
			"the desired trajectory does not exist.\n");
		return matchMembers;	// empty set
	}else{
		astrohelion::printColor(GREEN, "Located %zu matches; applying "
			"corrections\n", matches.size());
	}
	
	for(unsigned int n = 0; n < matches.size(); n++){
		unsigned int idx = matches[n];

		// Check to see if they are "close enough"
		if(std::abs(data->at(idx) - val) < matchTol){
			matchMembers.push_back(members[idx]);
			printf("  Candidate %d is close enough; no corrections required\n", n);
		}else{	// If not, employ corrections
			printf("  Correcting candidate %d...\n", idx);

			MultShootEngine corrector;
			corrector.setTol(1e-11);
			Arcset_periodic copyOrbit = members[idx];

			// Delete any constraints that might directly conflict with matchCon
			// There may still be indirect conflicts... ammend as necessary if 
			// those cases emerge
			std::vector<Constraint> cons;
			switch(matchCon.getAppType()){
				case ConstraintApp_tp::APP_TO_NODE:
					// Get constraints on relevant node, then clear them
					cons = copyOrbit.getNodeRef(matchCon.getID()).getConstraints();
					copyOrbit.getNodeRef(matchCon.getID()).clearConstraints();
					break;
				case ConstraintApp_tp::APP_TO_SEG:
					// get constraints on relevant segment, then clear them
					cons = copyOrbit.getSegRef(matchCon.getID()).getConstraints();
					copyOrbit.getSegRef(matchCon.getID()).clearConstraints();
					break;
				case ConstraintApp_tp::APP_TO_ARC:
					// get constraints only on arc, then clear them
					cons = copyOrbit.getArcConstraints();
					copyOrbit.clearArcConstraints();
					break;
			}
			// Re-add any constraints that do not conflict with matchCon
			for(unsigned int i = 0; i < cons.size(); i++){
				if(!cons[i].conflicts(matchCon)){
					copyOrbit.addConstraint(cons[i]);
				}
			}

			// Add match con
			copyOrbit.addConstraint(matchCon);

			// Do the corrections
			Arcset_periodic newOrbit(pSysData);
			try{
				corrector.multShoot(&copyOrbit, &newOrbit);
				matchMembers.push_back(newOrbit);
			}catch(DivergeException &e){
				printErr("  Unable to converge on a periodic solution "
					"for candidate %d...\n", n);
			}
		}
	}

	return matchMembers;

}//====================================================

/**
 *  @brief Populate a vector with a single coordinate from the initial state
 *  of each family member
 * 
 *  @param coordIx index of the coordinate within the state vector
 *  @param data pointer to the storage vector
 */
void Family_PO::getCoord(unsigned int coordIx, std::vector<double> *data) const{
	if(data){
		data->clear();
		data->assign(members.size(), NAN);
		std::vector<double> &dataRef = *data;
		for(unsigned int i = 0; i < members.size(); i++){
			std::vector<double> ic = members[i].getStateByIx(0);
			if(coordIx < ic.size())
				dataRef[i] = ic[coordIx];
			else
				throw Exception("Family_PO::getCoord: coordIx is larger than "
					"the state size");
		}
	}
}//====================================================

/**
 *  @brief Reverse the order of the family members.
 *  @details Family member data objects, their eigenvalues, and their
 *  eigenvectors are swapped in-place to reverse the order.
 */
void Family_PO::reverseOrder(){
	if(members.size() == 0)
		return;

	std::reverse(std::begin(members), std::end(members));
	std::reverse(std::begin(memberEigVecs), std::end(memberEigVecs));

	if(memberEigVecs.size() != members.size())
		throw Exception("Family_PO::reverseOrder: sizes of members and "
			"memberEigVecs are not consistent");

	unsigned int nE = memberEigVecs[0].rows();
	unsigned int nM = memberEigVals.size()/nE;

	if(nM != members.size())
		throw Exception("Family_PO::reverseOrder: sizes of members and "
			"memberEigVals are not consistent");

	// Swap groups of eigenvalues without changing order within each group
	for(unsigned int m = 0; m < nM/2; m++){
		for(unsigned int i = 0; i < nE; i++){
			std::iter_swap(memberEigVals.begin() + nE*m + i, 
				memberEigVals.end() - nE*(m+1) + i);
		}
	}
}//====================================================

//-----------------------------------------------------------------------------
//      File I/O
//-----------------------------------------------------------------------------

/**
 *  @brief Read family data from a Matlab file
 * 
 *  @param filename path to the file
 *  @param refLaws Reference to a vector of ControlLaw pointers. As control 
 *  laws are read from the Matlab file, unique control laws are constructed and 
 *  allocated on the stack. The user must manually delete the ControlLaw objects 
 *  to avoid memory leaks.
 *  
 *  @param bReconstruct whether or not to reconstruct each arc as it is read 
 *  from memory. "Reconstruction" is the process of propagating each segment to 
 *  populate the full segment state history.
 */
void Family_PO::readFromMat(const char *filename,
	std::vector<ControlLaw*> &refLaws, bool bReconstruct){

	mat_t *pMatFile = Mat_Open(filename, MAT_ACC_RDONLY);
	if(pMatFile == nullptr){
		char msg[256];
		sprintf(msg, "Family_PO::readFromMat: Could not load family from %s",
			filename);
		throw Exception(msg);
	}


	loadMembers(pMatFile, refLaws, bReconstruct);
	loadEigVals(pMatFile);
	loadMiscData(pMatFile);
	try{
		loadEigVecs(pMatFile);
	}catch(Exception &e){
		printWarn("Family_PO::readFromMat: Could not load eigenvectors from %s",
			filename);
	}
}//====================================================

/**
 *  @brief Save the family to a Matlab file
 *  @param filename file name (or path)
 */
void Family_PO::saveToMat(const char *filename) const{
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, nullptr, MAT_FT_DEFAULT);
	if(nullptr == matfp){
		astrohelion::printErr("Fam_cr3bp::saveToMat: Error creating MAT file\n");
	}else{
		// save things
		saveMembers(matfp);
		saveMiscData(matfp);
		saveEigVals(matfp);
		saveEigVecs(matfp);
		pSysData->saveToMat(matfp);
		saveTimestampToFile(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Load eigenvalues from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the eigenvalues
 *	@param matFile a pointer to the data file in question
 *	@throws Exception if the variable cannot be loaded
 */
void Family_PO::loadEigVals(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, VARNAME_FAM_EIGVAL);
	if(matvar == nullptr){
		throw Exception("Could not read eigenvalues into family");
	}else{
		unsigned int numMembers = matvar->dims[0], numVals = matvar->dims[1];

		if(members.size() != numMembers){
			throw Exception("Family_PO::loadEigVals: # eigenvalues is not "
				"same as number of members");
		}

		if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = 
				static_cast<mat_complex_split_t *>(matvar->data);

			if(splitVals != nullptr){
				// splitVals holds two void pointers to the real and imaginary 
				// parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				// Read data from column-major order matrix, store in row-major 
				// order vector
				for(unsigned int i = 0; i < numMembers; i++){
					for(unsigned int j = 0; j < numVals; j++){
						cdouble temp(realParts[j*numMembers + i], 
							imagParts[j*numMembers + i]);
						memberEigVals.push_back(temp);
					}
				}
			}
		}else{
			throw Exception("Family_PO::loadEigVals: Incompatible data file: "
				"unsupported data type/class");
		}
	}
	Mat_VarFree(matvar);
}//=============================================

/**
 *	@brief Load eigenvectors from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the 
 *	eigenvectors
 *	@param pMatFile a pointer to the data file in question
 *	@throws Exception if the variable cannot be loaded
 */
void Family_PO::loadEigVecs(mat_t* pMatFile){
	matvar_t *pMatvar = Mat_VarRead(pMatFile, VARNAME_FAM_EIGVEC);
	if(pMatvar == nullptr){
		throw Exception("Family_PO::loadEigVecs: Could not read data vector");
	}else{
		unsigned int numSteps = pMatvar->dims[2];

		if(members.size() == 0){
			throw Exception("Family_PO::loadEigVecs: Member vector has not "
				"been initialized!");
		}

		if(numSteps != members.size() ){
			throw Exception("Family_PO::loadEigVecs: Eigenvector vector does "
				"not have the same number of elements as the member vector");
		}

		if(pMatvar->dims[0] != pMatvar->dims[1]){
			throw Exception("Family_PO::loadEigVecs: Incompatible data file: "
				"Eigenvector matrix is not square.");
		}

		const unsigned int nE = pMatvar->dims[0];	// number of eigenvalues per family member

		if(pMatvar->class_type == MAT_C_DOUBLE && pMatvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = 
				static_cast<mat_complex_split_t *>(pMatvar->data);

			if(splitVals != nullptr){
				// splitVals holds two void pointers to the real and imaginary 
				// parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				for(unsigned int i = 0; i < numSteps; i++){
					std::vector<cdouble> vecData(nE*nE, 0);
					for(unsigned int j = 0; j < nE*nE; j++){
						vecData[j] = cdouble(realParts[j*numSteps + i], 
							imagParts[j*numSteps + i]);
					}
					memberEigVecs.push_back(Eigen::Map<MatrixXRcd>(\
						&(vecData.front()), nE, nE));
				}
			}
		}else{
			throw Exception("Family_PO::loadEigVecs: Incompatible data file: "
				"unsupported data type/class");
		}
	}
	Mat_VarFree(pMatvar);
}//====================================================

/**
 *  @brief Load family members from a Matlab file
 * 
 *  @param pMatFile pointer to the open Matlab file
 *  @param refLaws Reference to a vector of ControlLaw pointers. As control laws 
 *  are read from the Matlab file, unique control laws are constructed and 
 *  allocated on the stack. The user must manually delete the ControlLaw objects 
 *  to avoid memory leaks.
 *  
 *  @param bReconstruct whether or not to reconstruct each arc as it is read 
 *  from memory. "Reconstruction" is the process of propagating each segment to 
 *  populate the full segment state history.
 */
void Family_PO::loadMembers(mat_t *pMatFile, std::vector<ControlLaw*> &refLaws, 
	bool bReconstruct){

	matvar_t *pStruct = Mat_VarRead(pMatFile, VARNAME_FAM_MEMBER);
	if(pStruct == nullptr){
		char msg[256];
		sprintf(msg, "Family_PO::readFromMat: "
			"Could not read variable %s from file", VARNAME_FAM_MEMBER);
		throw Exception(msg);
	}else{
		if(pStruct->class_type == MAT_C_STRUCT && 
			pStruct->data_type == MAT_T_STRUCT){

			unsigned int numStructs = (pStruct->dims[0])*(pStruct->dims[1]);
			members = std::vector<Arcset_periodic>(numStructs,
				Arcset_periodic(pSysData));

			for(unsigned int s = 0; s < numStructs; s++){
				Arcset_periodic arc(pSysData);
				arc.readFromStruct(pStruct, s, refLaws);

				if(bReconstruct){
					Arcset_periodic fullArc(pSysData);
					reconstructArc(&arc, &fullArc);
					members[s] = fullArc;
				}else{
					members[s] = arc;
				}
			}
		}else{
			throw Exception("Family_PO::loadMembers: "
				"Family member variable does not have structure type/class");
		}
	}
}//====================================================

/**
 *  @brief Load miscellaneous family data from a Matlab file
 * 
 *  @param pMatFile pointer to an open Matlab file
 */
void Family_PO::loadMiscData(mat_t *pMatFile){
	name = readStringFromMat(pMatFile, VARNAME_NAME, MAT_T_UINT8, MAT_C_CHAR);
	int type = static_cast<int>(readDoubleFromMat(pMatFile, VARNAME_SORTTYPE));
	sortType = static_cast<FamSort_tp>(type);
}//====================================================

/**
 *  @brief Save family members to a matlab file
 *  @details Each family member is saved as an entry in a structure array
 * 
 *  @param pMatFile pointer to an open Matlab file
 */
void Family_PO::saveMembers(mat_t *pMatFile) const{
	if(members.size() == 0)
		return;

	const char *fieldNames[11] = {VARNAME_LINKTABLE,
		VARNAME_NODESTATE, VARNAME_STATE_DERIV, VARNAME_NODETIME, VARNAME_NODECTRL,
		VARNAME_SEGSTATE, VARNAME_SEGTIME, VARNAME_TOF, VARNAME_STM, VARNAME_SEGCTRL,
		VARNAME_CONSTRAINTS};
	const unsigned int numMembers = members.size();
	size_t dims[2] = {numMembers, 1};

	matvar_t *pStruct = Mat_VarCreateStruct(VARNAME_FAM_MEMBER, 2, dims, 
		fieldNames, 11);
	if(pStruct == nullptr){
		throw Exception("Family_PO::saveMembers: Could not create family "
			"member structure array");
	}

	for(unsigned int m = 0; m < numMembers; m++){
		members[m].saveToStruct(pStruct, m, Save_tp::SAVE_FRAME);
	}

	saveVar(pMatFile, pStruct, VARNAME_FAM_MEMBER, MAT_COMPRESSION_NONE);
}//====================================================

/**
 *	@brief Save eigenvalue data to a mat file
 *	@param pMatFile a pointer to the mat file in question
 *	@throws Exception if a family member does not have six eigenvalues
 */
void Family_PO::saveEigVals(mat_t *pMatFile) const{
	if(members.size() > 0){
		const unsigned int nE = memberEigVecs[0].rows();

		if(members.size() * nE != memberEigVals.size()){
			char msg[128];
			sprintf(msg, "Family_PO::saveEigVals: "
				"eigenvalue storage vector has %zu elements; "
				"expects %u*%zu", memberEigVals.size(), nE, members.size());
			throw Exception(msg);
		}

		// Separate all eigenvalues into real and complex parts
		std::vector<double> realParts(memberEigVals.size());
		std::vector<double> imagParts(memberEigVals.size());
		for(unsigned int i = 0; i < members.size(); i++){

			for(unsigned int j = 0; j < nE; j++){
				realParts[j*members.size() + i] = std::real(memberEigVals[i*nE + j]);
				imagParts[j*members.size() + i] = std::imag(memberEigVals[i*nE + j]);
			}
		}

		// create a special variable for them
		mat_complex_split_t splitVals = {&(realParts[0]), &(imagParts[0])};

		size_t dims[2] = {members.size(), nE};
		matvar_t *matvar = Mat_VarCreate(VARNAME_FAM_EIGVAL, MAT_C_DOUBLE, 
			MAT_T_DOUBLE, 2, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, matvar, VARNAME_FAM_EIGVAL, 
			MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	@brief Save eigenvector data to a mat file
 *	@param pMatFile a pointer to the mat file in question
 */
void Family_PO::saveEigVecs(mat_t *pMatFile) const{
	if(members.size() > 0){
		
		if(members.size() != memberEigVecs.size())
			throw Exception("Family_PO::saveEigVecs: eigenvector storage "
				"vector does not contain the appropriate number of elements!");

		const unsigned int nE = memberEigVecs[0].rows();
		std::vector<double> allVec_real(members.size()*nE*nE);
		std::vector<double> allVec_imag(members.size()*nE*nE);

		for(unsigned int i = 0; i < members.size(); i++){
			// get the transpose of the STM matrix; we need to store it in column-major order
			// and it's currently in row-major order
			MatrixXRcd P = memberEigVecs[i].transpose();
			
			// Retrieve the data from the matrix
			cdouble *matData = P.data();
			
			// Store that data in our huge vector
			for(unsigned int j = 0; j < nE*nE; j++){
				allVec_real[i*nE*nE + j] = std::real(matData[j]);
				allVec_imag[i*nE*nE + j] = std::imag(matData[j]);
			}
		}
		mat_complex_split_t splitVals = {&(allVec_real[0]), &(allVec_imag[0])};

		size_t dims[3] = {nE, nE, members.size()};
		matvar_t *pMatVar = Mat_VarCreate(VARNAME_FAM_EIGVEC, MAT_C_DOUBLE, 
			MAT_T_DOUBLE, 3, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, pMatVar, VARNAME_FAM_EIGVEC, 
			MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	@brief Save other useful information to a matlab file
 *	@param pMatFile the destination matlab file
 */
void Family_PO::saveMiscData(mat_t *pMatFile) const{
	// sortType
	int type = static_cast<int>(sortType);
	size_t dims[2] = {1,1};
	matvar_t *typeVar = Mat_VarCreate(VARNAME_SORTTYPE, MAT_C_INT32, 
		MAT_T_INT32, 2, dims, &type, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(pMatFile, typeVar, VARNAME_SORTTYPE, 
		MAT_COMPRESSION_NONE);

	// name
	char name_str[128];
	std::strcpy(name_str, name.c_str());
	dims[1] = name.length();
	matvar_t *nameVar = Mat_VarCreate(VARNAME_NAME, MAT_C_CHAR, MAT_T_UINT8, 2, 
		dims, &(name_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(pMatFile, nameVar, VARNAME_NAME, MAT_COMPRESSION_NONE);
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Make a copy of the family
 *  @param f reference to the source family
 */
void Family_PO::copyMe(const Family_PO &f){
	Family::copyMe(f);
	members = f.members;
	memberEigVals = f.memberEigVals;
	memberEigVecs = f.memberEigVecs;
}//====================================================

}// End of astrhelion namespace