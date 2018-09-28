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

#include "Arcset_cr3bp_lt.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "SysData_cr3bp_lt.hpp"

namespace astrohelion{

/**
 * @brief Construct a CR3BP-LT family of periodic orbits
 * 
 * @param pSys pointer to a CR3BP-LT system data object
 */
Family_PO_cr3bp_lt::Family_PO_cr3bp_lt(const SysData_cr3bp_lt *pSys) : Family_PO_cr3bp(pSys) {}

/**
 * @brief Retrieve a family member by the low-thrust Hamiltonian value at the 
 * first node (index = 0)
 * 
 * @param H Desired low-thrust Hamiltonian value (nondimensional)
 * @return A vector of periodic orbits that have the desired low-thrust 
 * Hamiltonian value at the first node
 */
std::vector<Arcset_periodic> Family_PO_cr3bp_lt::getMemberByH_lt(double H) const{
	std::vector<double> allH(members.size(), NAN);
	for(unsigned int n = 0; n < members.size(); n++){
		Arcset_cr3bp_lt temp(members[n]);
		allH[n] = temp.getHltByIx(0);
	}

	Constraint con(Constraint_tp::HLT, 0, &H, 1);

	return getMatchingMember(H, &allH, con);
}//====================================================

/**
 * @brief Retrieve all family members with a 2D thrust angle that matches the
 * specified angle
 * @details [long description]
 * 
 * @param alpha desired 2D thrust angle, in radians and between -PI and +PI.
 * Only the 2D components of the thrust vector are considered, so 3D thrust will
 * still return results.
 * 
 * @return a vector of periodic orbits with the specified 2D thrust angle
 */
std::vector<Arcset_periodic> Family_PO_cr3bp_lt::getMemberBy2DThrustAngle(
	double alpha) const{

	std::vector<double> allAlpha(members.size(), NAN);
	char msg[238];
	unsigned int ctrl_dim = 0, ctrl_tp = ControlLaw::NO_CTRL;
	for(unsigned int n = 0; n < members.size(); n++){
		ControlLaw_cr3bp_lt *pLaw = static_cast<ControlLaw_cr3bp_lt *>(
			members[n].getSegRefByIx(0).getCtrlLaw());

		if(pLaw){
			if(n == 0){
				ctrl_tp = pLaw->getType();
				ctrl_dim = pLaw->getNumStates();
			}else if(pLaw->getType() != ctrl_tp){
				snprintf(msg, 238, "Family_PO_cr3bp_lt::getMemberBy2DThrustAngle: "
					"Member %u, Segment 0 has thrust type %s, but previous"
					" members had type %s; cannot get consistent thrust"
					" angle.", n, pLaw->getTypeString().c_str(),
					ControlLaw::typeToString(ctrl_tp).c_str());
				throw Exception(msg);
			}

			const SysData_cr3bp_lt *pSys = 
				static_cast<const SysData_cr3bp_lt *>(pSysData);
			std::vector<double> q = 
				members[n].getSegRefByIx(0).getStateByRow(0);
			double t = members[n].getSegRefByIx(0).getTimeByIx(0);
			double a[3] = {0};
			pLaw->getOutput(t, &(q[0]), pSys, a, 3);

			// Get the thrust angle from the planar components of a
			if(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] > 1e-12){
				allAlpha[n] = atan2(a[1], a[0]);
			}else{
				snprintf(msg, 238, "Family_PO_cr3bp_lt::getMemberBy2DThrustAngle: "
					"Member %u, Segment 0  has thrust with zero magnitude;"
					" cannot get thrust angle.", n);
				throw Exception(msg);
			}
		}else{
			snprintf(msg, 238, "Family_PO_cr3bp_lt::getMemberBy2DThrustAngle: "
				"Member %u, Segment 0 has a null control law; cannot"
				" get thrust angle.", n);
			throw Exception(msg);
		}
	}

	std::vector<double> ctrl0(ctrl_dim, 0);
	ctrl0[0] = alpha;
	Constraint con(Constraint_tp::CTRL, 0, ctrl0);

	return getMatchingMember(alpha, &allAlpha, con);
}//====================================================

}// End of astrohelion namespace