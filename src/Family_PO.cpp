/**
 *  \file Family_PO.cpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version October 5, 2017
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

#include "Family_PO.hpp"

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "DynamicsModel.hpp"
#include "Exceptions.hpp"
#include "MultShootEngine.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors and Desctructor
//-----------------------------------------------------------------------------

Family_PO::Family_PO(const SysData *pSys) : Family(pSys) {}

Family_PO::Family_PO(const Family_PO &f): Family(f) { copyMe(f); }

Family_PO::~Family_PO() {}

//-----------------------------------------------------------------------------
//      Operators
//-----------------------------------------------------------------------------

Family_PO& Family_PO::operator= (const Family_PO &f){
	copyMe(f);
	return *this;
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Add a periodic orbit to the family
 *  \param arc reference to a periodic orbit; a copy is made and
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
 *  \brief Retrieve a copy of a family member
 * 
 *  \param ix index of the family member within the storage vector
 *  \return a copy of a family member
 */
Arcset_periodic Family_PO::getMember(int ix) const{
	while(ix < 0)
		ix += members.size();

	return members[ix];
}//====================================================

/**
 *  \brief Retrieve a reference to a family member
 * 
 *  \param ix index of the family member within the storage vector
 *  \return a reference to the family member
 */
Arcset_periodic& Family_PO::getMemberRef(int ix){
	while(ix < 0)
		ix += members.size();

	return members[ix];
}//====================================================

/**
 *  \brief Find a family member with a specific value of one of the initial states
 * 
 *  \param val the value of the initial state at the specified index, nondimensional units
 *  \param stateIx index of the state variable (e.g., x = 0, y = 1, etc.)
 * 
 *  \return A vector of family members that match the criterion
 */
std::vector<Arcset_periodic> Family_PO::getMemberByState(double val, unsigned int stateIx) const{
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
 *  \brief Find a family member with the specified time-of-flight
 * 
 *  \param tof time-of-flight, nondimensional units
 *  \return a vector of family members with matching times-of-flight
 */
std::vector<Arcset_periodic> Family_PO::getMemberByTOF(double tof) const{
	std::vector<double> allTOF;
	for(unsigned int n = 0; n < members.size(); n++){
		allTOF.push_back(members[n].getTotalTOF());
	}

	Constraint tofCon(Constraint_tp::TOF_TOTAL, 0, &tof, 1);

	return getMatchingMember(tof, &allTOF, tofCon);
}//====================================================

/**
 *  \brief Retrieve the number of periodic orbits in the family
 *  \return the number of periodic orbits in the family
 */
unsigned int Family_PO::getNumMembers() const { return members.size(); }

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *	\brief Locate all bifurcations in the family by analyzing the 
 *	eigenvalues
 *	\details Eigenvalues MUST be sorted, or this will yield completely bogus
 *	results
 */
std::vector<unsigned int> Family_PO::findBifurcations(){
	// Find the eigenvalues that are exactly (in theory) equal to one
	double onesErrTotal[3] = {0};
	unsigned int nE = memberEigVecs[0].rows();

	for(unsigned int m = 0; m < members.size(); m++){
		std::vector<cdouble> eigs(memberEigVals.begin() + m*nE, memberEigVals.begin() + (m+1)*nE);
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
		std::vector<cdouble> eigs(memberEigVals.begin() + m*nE, memberEigVals.begin() + (m+1)*nE);
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
				if(i != onesIx && stab[i]*prevStab[i] < 0 && std::abs(stab[i] - prevStab[i]) > okErr){
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
 *	\brief Sort all members' eigenvalues so they are in the same order.
 *	\details This is necessary before bifurcations can be accurately located.
 */
void Family_PO::sortEigs(){
	if(memberEigVecs.size() == 0)
		return;

	// Sort eigenvalues
	const unsigned int nE = memberEigVecs[0].rows();	// number of eigenvalues
	std::vector<unsigned int> sortedIxs = sortEig(memberEigVals, memberEigVecs);

	// Update all eigendata in the storage vectors
	std::vector<cdouble> sortedVals;
	MatrixXRcd sortedVecs = MatrixXRcd::Zero(nE, nE);
	unsigned int ix = 0, m = 0, c = 0;
	for(m = 0; m < members.size(); m++){
		// Construct new eigenvalue and eigenvector objects with the new index order
		for(c = 0; c < nE; c++){
			ix = sortedIxs[m*nE+c];
			sortedVals.push_back(memberEigVals[m*nE+ix]);

			sortedVecs.col(c) = memberEigVecs[m].col(ix);
		}
		
		// Update member with sorted vectors and values
		memberEigVecs[m] = sortedVecs;
		for(c = 0; c < nE; c++){
			memberEigVals[c] = sortedVals[c];
		}

		// Reset storage variables
		sortedVals.clear();
		sortedVecs = MatrixXRcd::Zero(nE,nE);
	}
}//====================================================

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
			throw Exception("Family_PO::sortMembers: JC Sorting not implemented here!");
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
		sortedEigVals.insert(sortedEigVals.end(), memberEigVals.begin() + nE*ix, memberEigVals.begin() + nE*(ix+1));
	}

	// Reassign the members vector
	members = sortedMembers;
	memberEigVecs = sortedEigVecs;
	memberEigVals = sortedEigVals;
}//====================================================

/**
 *	\brief Locate a family member with a specific attribute
 *	
 *	This function locates a family member or set of members that have a specific value
 *	for one of the variables of interest (e.g. coordinates, Jacobi, TOF). Exact matches
 *	and interpolated matches are returned; interpolated matches are computed using a
 *	differential corrections algorithm.
 *
 *	\param val the value the family member should have
 *	\param data a pointer to a vector containing the set of values to search for matches
 *	in. For example, if the <tt>val</tt> I pass in contains a specific TOF, then 
 *	<tt>data</tt> points to a vector containing the TOFs for the entire family, sorted
 *	according to this family's <tt>sortType</tt>.
 *	\param matchCon a constraint that can be applied in a corrections scheme that will
 *	ensure the corrected trajectory has the desired value for the variable of interest.
 *
 *	\return a vector of matches. If no matches are returned, the vector will be empty.
 */
std::vector<Arcset_periodic> Family_PO::getMatchingMember(double val,
	std::vector<double> *data, Constraint matchCon) const{

	// Locate possible candidates
	std::vector<unsigned int> matches = findMatches(val, data);
	std::vector<Arcset_periodic> matchMembers;

	if(matches.size() == 0){
		astrohelion::printErr("Could not locate any matches. The family either has too few members to facilitate an accurate search or the desired trajectory does not exist.\n");
		return matchMembers;	// empty set
	}else{
		astrohelion::printColor(GREEN, "Located %zu matches; applying corrections\n", matches.size());
	}
	
	for(unsigned int n = 0; n < matches.size(); n++){
		unsigned int idx = matches[n];

		// Check to see if they are "close enough"
		if(std::abs(data->at(idx) - val) < matchTol){
			matchMembers.push_back(members[idx]);
			printf("  Candidate %d is close enough; no corrections required\n", n);
		}else{	// If not, employ corrections
			printf("  Correcting candidate %d...\n", n);

			MultShootEngine corrector;
			corrector.setTol(1e-11);
			Arcset_periodic copyOrbit = members[idx];
			copyOrbit.addConstraint(matchCon);

			Arcset_periodic newOrbit(pSysData);

			try{
				corrector.multShoot(&copyOrbit, &newOrbit);
				matchMembers.push_back(newOrbit);
			}catch(DivergeException &e){
				printErr("  Unable to converge on a periodic solution for candidate %d...\n", n);
			}
		}
	}

	return matchMembers;

}//====================================================

/**
 *  \brief Populate a vector with a single coordinate from the initial state
 *  of each family member
 * 
 *  \param coordIx index of the coordinate within the state vector
 *  \param data pointer to the storage vector
 */
void Family_PO::getCoord(unsigned int coordIx, std::vector<double> *data) const{
	if(data){
		data->clear();
		for(unsigned int i = 0; i < members.size(); i++){
			std::vector<double> ic = members[i].getStateByIx(0);
			if(coordIx < ic.size())
				data->push_back(ic[coordIx]);
			else
				throw Exception("Family_PO::getCoord: coordIx is larger than the state size");
		}
	}
}//====================================================

//-----------------------------------------------------------------------------
//      File I/O
//-----------------------------------------------------------------------------

/**
 *  \brief Save the family to a file
 *  \param filename file name (or path)
 */
void Family_PO::saveToMat(const char *filename){
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
		// saveMembers(matfp);
		saveMiscData(matfp);
		saveEigVals(matfp);
		saveEigVecs(matfp);
		pSysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	\brief Load eigenvalues from the data file
 *
 *	NOTE: the vector of family members MUST be populated before loading the eigenvalues
 *	\param matFile a pointer to the data file in question
 *	\throws Exception if the variable cannot be loaded
 */
void Family_PO::loadEigVals(mat_t *matFile){
	matvar_t *matvar = Mat_VarRead(matFile, EIG_VAR_NAME);
	if(matvar == NULL){
		throw Exception("Could not read eigenvalues into family");
	}else{
		unsigned int numMembers = matvar->dims[0];
		if(matvar->dims[1] != 6){
			throw Exception("Incompatible data file: Data widths are different.");
		}

		if(members.size() != numMembers){
			throw Exception("Family_PO::loadEigVals: # eigenvalues is not same as number of members");
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
					for(int j = 0; j < 6; j++){
						cdouble temp(realParts[j*numMembers + i], imagParts[j*numMembers + i]);
						memberEigVals.push_back(temp);
					}
				}
			}
		}else{
			throw Exception("Family_PO::loadEigVals: Incompatible data file: unsupported data type/class");
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
void Family_PO::loadEigVecs(mat_t* pMatFile){
	matvar_t *pMatvar = Mat_VarRead(pMatFile, EIGVEC_VAR_NAME);
	if(pMatvar == NULL){
		throw Exception("Family_PO::loadEigVecs: Could not read data vector");
	}else{
		unsigned int numSteps = pMatvar->dims[2];

		if(members.size() == 0){
			throw Exception("Family_PO::loadEigVecs: Member vector has not been initialized!");
		}

		if(numSteps != members.size() ){
			throw Exception("Family_PO::loadEigVecs: Eigenvector vector does not have the same number of elements as the member vector");
		}

		if(pMatvar->dims[0] != pMatvar->dims[1]){
			throw Exception("Family_PO::loadEigVecs: Incompatible data file: Eigenvector matrix is not square.");
		}

		const unsigned int nE = pMatvar->dims[0];	// number of eigenvalues per family member

		if(pMatvar->class_type == MAT_C_DOUBLE && pMatvar->data_type == MAT_T_DOUBLE){
			// First cast the data to a special variable matio uses to store complex values
			mat_complex_split_t *splitVals = static_cast<mat_complex_split_t *>(pMatvar->data);

			if(splitVals != NULL){
				// splitVals holds two void pointers to the real and imaginary parts; cast them to doubles
				double *realParts = static_cast<double *>(splitVals->Re);
				double *imagParts = static_cast<double *>(splitVals->Im);

				for(unsigned int i = 0; i < numSteps; i++){
					std::vector<cdouble> vecData(nE*nE, 0);
					for(unsigned int j = 0; j < 36; j++){
						vecData[j] = cdouble(realParts[j*numSteps + i], imagParts[j*numSteps + i]);
					}
					memberEigVecs.push_back(Eigen::Map<MatrixXRcd>(&(vecData.front()), 6, 6));
				}
			}
		}else{
			throw Exception("Family_PO::loadEigVecs: Incompatible data file: unsupported data type/class");
		}
	}
	Mat_VarFree(pMatvar);
}//=============================================

/**
 *	\brief Save eigenvalue data to a mat file
 *	\param pMatFile a pointer to the mat file in question
 *	\throws Exception if a family member does not have six eigenvalues
 */
void Family_PO::saveEigVals(mat_t *pMatFile){
	if(members.size() > 0){
		const unsigned int nE = memberEigVecs[0].rows();

		if(members.size() * nE != memberEigVals.size())
			throw Exception("Family_PO::saveEigVals: eigenvalue storage vector does not contain the appropriate number of elements!");

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

		size_t dims[2] = {members.size(), 6};
		matvar_t *matvar = Mat_VarCreate(EIG_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, matvar, EIG_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save eigenvector data to a mat file
 *	\param pMatFile a pointer to the mat file in question
 */
void Family_PO::saveEigVecs(mat_t *pMatFile){
	if(members.size() > 0){
		
		if(members.size() != memberEigVecs.size())
			throw Exception("Family_PO::saveEigVecs: eigenvector storage vector does not contain the appropriate number of elements!");

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
		matvar_t *pMatVar = Mat_VarCreate(EIGVEC_VAR_NAME, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, &splitVals, MAT_F_COMPLEX);
		astrohelion::saveVar(pMatFile, pMatVar, EIGVEC_VAR_NAME, MAT_COMPRESSION_NONE);
	}
}//====================================================

/**
 *	\brief Save other useful information to a matlab file
 *	\param pMatFile the destination matlab file
 */
void Family_PO::saveMiscData(mat_t *pMatFile){
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

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

void Family_PO::copyMe(const Family_PO &f){
	Family::copyMe(f);
	members = f.members;
	memberEigVals = f.memberEigVals;
	memberEigVecs = f.memberEigVecs;
}//====================================================

}// End of astrhelion namespace