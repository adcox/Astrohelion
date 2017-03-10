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

/**
 *  \brief Construct a default nodeset
 *  \param pSys pointer to system data object
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys) : Nodeset(pSys){}

/**
 *  \brief Construct a nodeset via numerical integration
 * 
 *  \param pSys pointer to system data object
 *  \param ic initial condition on the nodeset (7 elements)
 *  \param tof total time-of-flight on the nodeset
 *  \param numNodes number of nodes along the specified TOF
 *  \param type Nodeset type; default is TIME
 *  \param ctrlLaw control law applied; default is NO_CTRL
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys, const double ic[7], 
	double tof, int numNodes, NodeDistro_tp type, unsigned int ctrlLaw) : Nodeset(pSys){

	std::vector<double> ic_vec(ic, ic+7);
	initFromICs(ic_vec, 0, tof, numNodes, type, ctrlLaw);
}//====================================================

/**
 *  \brief Construct a nodeset via numerical integration
 *  \details [long description]
 * 
 *  \param pSys pointer to system data object
 *  \param ic initial condition on the nodeset
 *  \param tof total time-of-flight on the nodeset
 *  \param numNodes number of nodes along the specified TOF
 *  \param type Nodeset type; default is TIME
 *  \param ctrlLaw control law applied; default is NO_CTRL
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const SysData_cr3bp_lt *pSys, std::vector<double> ic, 
	double tof, int numNodes, NodeDistro_tp type, unsigned int ctrlLaw) : Nodeset(pSys){

	initFromICs(ic, 0, tof, numNodes, type, ctrlLaw);
}//====================================================

/**
 *  \brief Construct a nodeset from a trajectory
 *  \details [long description]
 * 
 *  \param traj Trajectory to break into nodes and segments
 *  \param numNodes Number of nodes to use
 *  \param type how the nodes are placed
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(Traj_cr3bp_lt traj, int numNodes,
	NodeDistro_tp type) : Nodeset(traj.getSysData()){

	initFromTraj(traj, numNodes, type);
}//====================================================

/**
 *  \brief Create a nodeset by clipping another between two specified points
 *  \details NOT YET IMPLEMENTED
 * 
 *  \param orig [description]
 *  \param first [description]
 *  \param last [description]
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt &orig, int first, int last) : Nodeset(orig, first, last) {}

/**
 *  \brief Copy Constructor
 *  \param n Another nodeset
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const Nodeset_cr3bp_lt &n) : Nodeset(n) {}

/**
 *  \brief Copy constructor
 *  \param a the base class
 */
Nodeset_cr3bp_lt::Nodeset_cr3bp_lt(const BaseArcset &a) : Nodeset(a) {}

/**
 *  \brief Create a new nodeset object on the stack
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param pSys pointer to a system data object
 *  \return a pointer to the newly created nodeset
 */
baseArcsetPtr Nodeset_cr3bp_lt::create( const SysData *pSys) const{
	const SysData_cr3bp_lt *crSys = static_cast<const SysData_cr3bp_lt*>(pSys);
	return baseArcsetPtr(new Nodeset_cr3bp_lt(crSys));
}//====================================================

/**
 *  \brief Create a new nodeset object on the stack that is a 
 *  duplicate of this object
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \return a pointer to the newly cloned nodeset
 */
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