/**
 * \file Nodeset_cr3bp_lt.cpp
 * \brief Derivative of Nodeset object, specific to low-thrust CR3BP
 */
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "Nodeset_cr3bp_lt.hpp"

#include "SysData_cr3bp_lt.hpp"
#include "Traj_cr3bp_lt.hpp"

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys) : Nodeset(pSys){}

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys, const double ic[7], 
	double tof, int numNodes, NodeDistro_tp type, unsigned int ctrlLaw) : Nodeset(pSys){

	std::vector<double> ic_vec(ic, ic+7);
	initFromICs(ic_vec, 0, tof, numNodes, type, ctrlLaw);
}//====================================================

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys, std::vector<double> ic, 
	double tof, int numNodes, NodeDistro_tp type, unsigned int ctrlLaw) : Nodeset(pSys){

	initFromICs(ic, 0, tof, numNodes, type, ctrlLaw);
}//====================================================

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(Traj_cr3bp_lt traj, int numNodes,
	NodeDistro_tp type) : Nodeset(traj.getSysData()){

	initFromTraj(traj, numNodes, type);
}//====================================================

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt &orig, int first, int last) : Nodeset(orig, first, last) {}

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt &n) : Nodeset(n) {}

Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const BaseArcset &a) : Nodeset(a) {}

baseArcsetPtr Nodeset_cr3bp_lt::create( const SysData *pSys) const{
	const SysData_cr3bp_lt *crSys = static_cast<const SysData_cr3bp_lt*>(pSys);
	return baseArcsetPtr(new Nodeset_cr3bp_lt(crSys));
}//====================================================

baseArcsetPtr Nodeset_cr3bp_lt::clone() const{
	return baseArcsetPtr(new Nodeset_cr3bp_lt(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

}