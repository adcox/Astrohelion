/**
 *  @file Family_PO_cr3bp_lt.cpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version March 8, 2018
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

#include "Family_PO_cr3bp_lt.hpp"

#include "SysData_cr3bp_lt.hpp"

namespace astrohelion{

/**
 * @brief Construct a CR3BP-LT family of periodic orbits
 * 
 * @param pSys pointer to a CR3BP-LT system data object
 */
Family_PO_cr3bp_lt::Family_PO_cr3bp_lt(const SysData_cr3bp_lt *pSys) : Family_PO_cr3bp(pSys) {}

std::vector<Arcset_periodic> Family_PO_cr3bp_lt::getMemberByH_lt(double H) const{
	std::vector<double> allH(members.size(), NAN);
	for(unsigned int n = 0; n < members.size(); n++){
		Arcset_cr3bp_lt temp(members[n]);
		allH[n] = temp.getHltByIx(0);
	}

	Constraint con(Constraint_tp::HLT, 0, &H, 1);

	return getMatchingMember(H, &allH, con);
}//====================================================

}// End of astrohelion namespace