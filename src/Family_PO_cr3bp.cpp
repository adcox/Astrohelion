/**
 *  @file Family_PO_cr3bp.cpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version October 18, 2017
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

#include "Family_PO_cr3bp.hpp"

#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "Constraint.hpp"
#include "SysData_cr3bp.hpp"

namespace astrohelion{

/**
 * @brief Construct a CR3BP family of periodic orbits
 * 
 * @param pSys pointer to a CR3BP system data object
 */
Family_PO_cr3bp::Family_PO_cr3bp(const SysData_cr3bp *pSys) : Family_PO(pSys) {}

/**
 * @brief Retrieve all family members that have Jacobi constant values matching
 * the specified value
 * 
 * @param JC Desired Jacobi constant value
 * @return A vector of all matches located. If no matches are located, an empty
 * vector is returned.
 */
std::vector<Arcset_periodic> Family_PO_cr3bp::getMemberByJacobi(double JC) const {
	// Get an array of all the jacobi values
	std::vector<double> allJC(members.size(), NAN);
	for(unsigned int n = 0; n < members.size(); n++){
		Arcset_cr3bp temp(members[n]);
		allJC[n] = temp.getJacobiByIx(0);
	}

	Constraint jacobiCon(Constraint_tp::JC, 0, &JC, 1);

	return getMatchingMember(JC, &allJC, jacobiCon);
}//====================================================

}// End of astrohelion namespace