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

#include "Family_PO_cr3bp.hpp"

#include "Arcset_cr3bp.hpp"
#include "Arcset_periodic.hpp"
#include "Constraint.hpp"
#include "SysData_cr3bp.hpp"

namespace astrohelion{

Family_PO_cr3bp::Family_PO_cr3bp(const SysData_cr3bp *pSys) : Family_PO(pSys) {}

std::vector<Arcset_periodic> Family_PO_cr3bp::getMemberByJacobi(double JC) const {
	// Get an array of all the jacobi values
	std::vector<double> allJC;
	for(unsigned int n = 0; n < members.size(); n++){
		Arcset_cr3bp temp(members[n]);
		allJC.push_back(temp.getJacobiByIx(0));
	}

	Constraint jacobiCon(Constraint_tp::JC, 0, &JC, 1);

	return getMatchingMember(JC, &allJC, jacobiCon);
}//====================================================

}// End of astrohelion namespace