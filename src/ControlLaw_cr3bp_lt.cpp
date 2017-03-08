/**
 * @file ControlLaw_cr3bp_lt.hpp
 * @brief Control Law for CR3BP-LT system header file 
 * 
 * @author Andrew Cox
 * @version March 3, 2017
 * @copyright GNU GPL v3.0
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
 
#include "ControlLaw_cr3bp_lt.hpp"

#include <cmath>

#include "Exceptions.hpp"
#include "SysData_cr3bp_lt.hpp"

namespace astrohelion{

ControlLaw_cr3bp_lt::ControlLaw_cr3bp_lt(){}

//------------------------------------------------------------------------------------------------------
//      Switchboard Functions
//------------------------------------------------------------------------------------------------------

void ControlLaw_cr3bp_lt::getLaw(double t, const double *s, const SysData *pSysData, unsigned int lawID,
	double *law, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSysData);

	switch(lawID){
		case Law_tp::CONST_C_2D_LEFT:
			getLaw_ConstC_2D_Left(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getLaw_ConstC_2D_Right(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::PRO_VEL:
			getLaw_Pro_Vel(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::ANTI_VEL:
			getLaw_Anti_Vel(t, s, pSysData_lt, law, len);
			break;
		default:
			ControlLaw::getLaw(t, s, pSysData, lawID, law, len);
	}
}//====================================================

/**
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param t [description]
 *  \param s [description]
 *  \param pSys [description]
 *  \param int [description]
 *  \param partials [description]
 *  \param int [description]
 */
void ControlLaw_cr3bp_lt::getPartials_State(double t, const double *s, const SysData *pSys, unsigned int lawID, double *partials, unsigned int len) const{
	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSys);
	switch(lawID){
		case Law_tp::CONST_C_2D_LEFT:
			getPartials_State_ConstC_2D(t, s, pSysData_lt, partials, len, -1);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getPartials_State_ConstC_2D(t, s, pSysData_lt, partials, len, 1);
			break;
		case Law_tp::PRO_VEL:
			// getLaw_Pro_Vel(t, s, pSysData, law, len);
			// break;
		case Law_tp::ANTI_VEL:
			// getLaw_Anti_Vel(t, s, pSysData, law, len);
			// break;
		default:
			ControlLaw::getPartials_State(t, s, pSys, lawID, partials, len);
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Control Laws
//------------------------------------------------------------------------------------------------------

void ControlLaw_cr3bp_lt::getLaw_ConstC_2D_Right(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_ConstC_2D_Right: law data length must be at least 3!");

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = s[4]/v;
	law[1] = -s[3]/v;
	law[2] = 0;

	(void) pSys;
	(void) t;
}//====================================================

void ControlLaw_cr3bp_lt::getLaw_ConstC_2D_Left(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_ConstC_2D_Left: law data length must be at least 3!");

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = -s[4]/v;
	law[1] = s[3]/v;
	law[2] = 0;

	(void) pSys;
	(void) t;
}//====================================================

void ControlLaw_cr3bp_lt::getLaw_Pro_Vel(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_Pro_Vel: law data length must be at least 3!");

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = s[3]/v;
	law[1] = s[4]/v;
	law[2] = s[5]/v;

	(void) pSys;
	(void) t;
}//====================================================

void ControlLaw_cr3bp_lt::getLaw_Anti_Vel(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_Anti_Vel: law data length must be at least 3!");

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = -s[3]/v;
	law[1] = -s[4]/v;
	law[2] = -s[5]/v;

	(void) pSys;
	(void) t;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Partial Derivatives of Control Laws
//------------------------------------------------------------------------------------------------------

void ControlLaw_cr3bp_lt::getPartials_State_ConstC_2D(double t, const double *s, const SysData_cr3bp_lt *pSys, double *partials, unsigned int len, int sign) const{
	if(std::abs(sign) != 1)
		sign = 1;	// +1 for RIGHT, -1 for LEFT

	if(len != 21)
		throw Exception("ControlLaw_cr3bp_lt::getPartials_State_ConstC_2D: Expects len = 21");

	// state s : [x, y, z, vx, vy, vz, m, ... stm_elements ...]
	// partials: row 1 = partials of a_x w.r.t. states, row 2 = partials of a_y w.r.t. states,
	// row 3 = partials of a_z w.r.t. states

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	double f = pSys->getThrust();

	// Initialize all partials to zero
	for(unsigned int i = 0; i < len; i++){ partials[i] = 0; }

	// dax/dvx
	partials[7*0 + 3] = -sign*f*s[3]*s[4]/(s[6]*pow(v,3));							// dax/dvx
	partials[7*0 + 4] = sign*(f/(s[6]*v) - f*s[4]*s[4]/(s[6]*pow(v,3)));			// dax/dvy
	partials[7*0 + 6] = -sign*f*s[4]/(s[6]*s[6]*v);									// dax/dm

	partials[7*1 + 3] = -sign*(f/(s[6]*v) - f*s[3]*s[3]/(s[6]*pow(v,3)));			// day/dvx
	partials[7*1 + 4] = sign*f*s[3]*s[4]/(s[6]*pow(v,3));							// day/dvy
	partials[7*1 + 6] = sign*f*s[3]/(s[6]*s[6]*v);									// day/dm

	(void) t;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

std::string ControlLaw_cr3bp_lt::lawIDToString(unsigned int id) const{
	switch(id){
		case Law_tp::CONST_C_2D_LEFT: return "Jacobi-Preserving, 2D, Left";
		case Law_tp::CONST_C_2D_RIGHT: return "Jacobi-Preserving, 2D, Right";
		case Law_tp::PRO_VEL: return "Prograde Velocity";
		case Law_tp::ANTI_VEL: return "Anti-Velocity";
		default:
			return ControlLaw::lawIDToString(id);
	}
}//====================================================

}
