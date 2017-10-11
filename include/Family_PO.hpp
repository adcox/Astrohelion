/**
 * \file Family_PO.hpp
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrohelion.
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

#include <vector>

#include "Family.hpp"

#include "Arcset_periodic.hpp"

namespace astrohelion{

/**
 *  \ingroup traj
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param  [description]
 */
class Family_PO : public Family{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		Family_PO(const SysData*);
		Family_PO(const Family_PO&);
		virtual ~Family_PO();
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		Family_PO& operator= (const Family_PO&);
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		void addMember(const Arcset_periodic&);
		std::vector<unsigned int> findBifurcations();
		void saveToMat(const char*);
		void sortEigs();
		void sortMembers();
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		Arcset_periodic getMember(int) const;
		Arcset_periodic& getMemberRef(int);
		std::vector<Arcset_periodic> getMemberByState(double, unsigned int) const;
		std::vector<Arcset_periodic> getMemberByTOF(double) const;
		unsigned int getNumMembers() const;
		//\}
	protected:
		std::vector<Arcset_periodic> members {};		//!< Vector of periodic orbits
		std::vector<cdouble> memberEigVals {};			//!< Vector of eigenvalues (row-major order)
		std::vector<MatrixXRcd> memberEigVecs {};		//!< Vector of eigenvectors matrices

		const char* VARNAME_FAM_MEMBER = "Members";			//!< Variable name for the family members
		const char* VARNAME_FAM_EIGVAL = "Eigenvalues";		//!< Variable name for the eigenvalues
		const char* VARNAME_FAM_EIGVEC = "Eigenvectors";	//!< Variable name for the eigenvectors

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		std::vector<Arcset_periodic> getMatchingMember(double, std::vector<double>*, Constraint) const;
		void getCoord(unsigned int, std::vector<double>*) const;
		//\}

		/**
		 *  \name File I/O
		 *  \{
		 */
		// void loadMemberData(mat_t*);
		void loadEigVals(mat_t*);
		void loadEigVecs(mat_t*);
		
		void saveMembers(mat_t*);
		void saveMiscData(mat_t*);
		void saveEigVals(mat_t*);
		void saveEigVecs(mat_t*);
		//\}

		/**
		 *  \name File I/O
		 *  \{
		 */
		void copyMe(const Family_PO&);
		//\}
};

}// END of astrohelion namespace