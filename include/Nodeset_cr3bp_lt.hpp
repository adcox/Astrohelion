/**
 * @file Nodeset_cr3bp_lt.hpp
 *
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

#pragma once

#include "Nodeset.hpp"

namespace astrohelion{

// Forward declarations
class SysData_cr3bp_lt;
class Traj_cr3bp_lt;
/**
 *  \ingroup traj cr3bp_lt
 *  \brief This derivative of the Nodeset contains information for the 
 *  low-thrust CR3BP
 *  
 *  \details [long description]
 * 
 * 	\author Andrew Cox
 * 	\version March 8, 2017
 * 	\copyright GNU GPL v3.0
 */
class Nodeset_cr3bp_lt : public Nodeset{

public:
	/**
	 *  \name *structors
	 *  \{
	 */
	Nodeset_cr3bp_lt(const SysData_cr3bp_lt*);
	Nodeset_cr3bp_lt(const SysData_cr3bp_lt*, const double[7], double, int, NodeDistro_tp type = NodeDistro_tp::TIME);
	Nodeset_cr3bp_lt(const SysData_cr3bp_lt*, std::vector<double>, double, int, NodeDistro_tp type = NodeDistro_tp::TIME);
	Nodeset_cr3bp_lt(Traj_cr3bp_lt, int, NodeDistro_tp);

	Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt&, int, int);
	Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt&);
	Nodeset_cr3bp_lt(const BaseArcset&);
	baseArcsetPtr create(const SysData*) const override;
	baseArcsetPtr clone() const override;

	//\}

};

}// End of astrohelion namespace