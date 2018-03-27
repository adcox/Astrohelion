/**
 * @file Family_PO_cr3bp.hpp
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Family_PO_cr3bp.hpp"

namespace astrohelion{

// Forward Declarations
class SysData_cr3bp_lt;

/**
 *  \ingroup fam cr3bp_lt
 *  @brief [brief description]
 *  @details [long description]
 * 
 *  @param  [description]
 *  @return [description]
 */
class Family_PO_cr3bp_lt : public Family_PO_cr3bp{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		Family_PO_cr3bp_lt(const SysData_cr3bp_lt*);
		//\}
		
		/**
		 *  \name Analysis Functions
		 *  \
		 */
		std::vector<Arcset_periodic> getMemberByH_lt(double) const;		
		std::vector<Arcset_periodic> getMemberByThrustMag(double) const;
		std::vector<Arcset_periodic> getMemberByThrustAngle(double, double beta = 0) const;
		//\}
};

}// End of astrohelion namespace