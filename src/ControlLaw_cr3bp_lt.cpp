/**
 * \file ControlLaw_cr3bp_lt.hpp
 * \brief Control Law for CR3BP-LT system header file 
 * 
 * \author Andrew Cox
 * \version March 3, 2017
 * \copyright GNU GPL v3.0
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
#include "Utilities.hpp"

namespace astrohelion{

//------------------------------------------------------------------------------------------------------
//      Constructors
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Construct a default CR3BP low-thrust control law object
 *  
 *  \param id Control Law ID
 *  \param params a vector of parameters used by the control law. These parameters 
 *  must be thrust (in Newtons) and Specific Impulse (in seconds).
 */
ControlLaw_cr3bp_lt::ControlLaw_cr3bp_lt(unsigned int id, std::vector<double> params) : ControlLaw(id, params){}

/**
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param id Control law ID
 *  \param T Thrust value, Newtons
 *	\param I Specific Impulse (Isp), seconds
 */
ControlLaw_cr3bp_lt::ControlLaw_cr3bp_lt(unsigned int id, double T, double Isp){
	lawID = id;
	params.assign(2,0);
	params[0] = T;
	params[1] = Isp;

	init();
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Set and Get Functions
//------------------------------------------------------------------------------------------------------

/**
 *	\brief Get the spacecraft thrust in Newtons
 *	\return the thrust in Newtons
 */
double ControlLaw_cr3bp_lt::getThrust() const { return params[0]; }

/**
 *	\brief Get the specific impulse for the spacecraft
 *	\return the specific impulse, seconds
 */
double ControlLaw_cr3bp_lt::getIsp() const { return params[1]; }

/**
 *	\brief Set the spacecraft thrust
 *	\param f the thrust, in Newtons
 */
void ControlLaw_cr3bp_lt::setThrust(double f){ params[0] = f;}

/**
 *	\brief Set the specific impulse for the spacecraft
 *	\param Isp the specific impulse, seconds
 */
void ControlLaw_cr3bp_lt::setIsp(double Isp){ params[1] = Isp; }

//------------------------------------------------------------------------------------------------------
//      Dynamics Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the output of a control law
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSysData system data object
 *  \param lawID identifies the control law type
 *  \param law empty, initialized array to store the control law output in
 *  \param len the number of elements in the <tt>law</tt> array
 *  
 *  \throws Exception if the control law ID, <tt>lawID</tt>, is not recognized
 */
void ControlLaw_cr3bp_lt::getLaw_Accel(double t, const double *s, const SysData *pSysData,
	double *law, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSysData);

	switch(lawID){
		case Law_tp::CONST_C_2D_LEFT:
			getAccel_ConstC_2D(t, s, pSysData_lt, law, len, -1);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getAccel_ConstC_2D(t, s, pSysData_lt, law, len, 1);
			break;
		case Law_tp::PRO_VEL:
			getAccel_Along_Vel(t, s, pSysData_lt, law, len, 1);
			break;
		case Law_tp::ANTI_VEL:
			getAccel_Along_Vel(t, s, pSysData_lt, law, len, -1);
			break;
		case Law_tp::GENERAL_CONST_F:
			getAccel_GeneralDir(t, s, pSysData_lt, law, len);
			break;
		default:
			ControlLaw::getLaw_Accel(t, s, pSysData, law, len);
	}
}//====================================================

/**
 *  \brief Retrieve the partial derivatives of the control law with respect to state variables
 *  \details A set of partial derivatives of the control law outputs are computed with respect to the 
 *  states at the given time, state, in the specified system
 * 
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param lawID identifies the control law type
 *  \param partials empty, initialized array to store the control law derivatives in
 *  \param len number of elements in the <tt>law</tt> array
 *  
 *  \throws Exception if the control law ID, <tt>lawID</tt>, is not recognized
 */
void ControlLaw_cr3bp_lt::getLaw_AccelPartials(double t, const double *s, const SysData *pSys, 
	double *partials, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSys);
	switch(lawID){
		case Law_tp::CONST_C_2D_LEFT:
			getAccelPartials_ConstC_2D(t, s, pSysData_lt, partials, len, -1);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getAccelPartials_ConstC_2D(t, s, pSysData_lt, partials, len, 1);
			break;
		case Law_tp::GENERAL_CONST_F:
			getAccelPartials_GeneralDir(t, s, pSysData_lt, partials, len);
			break;
		case Law_tp::PRO_VEL:
			// getLaw_Pro_Vel(t, s, pSysData, law, len);
			// break;
		case Law_tp::ANTI_VEL:
			// getLaw_Anti_Vel(t, s, pSysData, law, len);
			// break;
		
		default:
			ControlLaw::getLaw_AccelPartials(t, s, pSys, partials, len);
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Control Laws
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the output of the Jacobi-preserving 2D control laws
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param law empty, initialized array to store the control law output in
 *  \param len number of elements in the <tt>law</tt> array
 *  \param sign specifies which of the ConstC_2D control laws to evaluate: +1 for RIGHT, -1 for LEFT
 */
void ControlLaw_cr3bp_lt::getAccel_ConstC_2D(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len, int sign) const{

	if(std::abs(sign) != 1)
		sign = sign/std::abs(sign);	// +1 for RIGHT, -1 for LEFT

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_ConstC_2D: law data length must be at least 3!");

	double charT = pSys->getCharT();
	double charL = pSys->getCharL();
	double f = (getThrust()/1000)*charT*charT/charL/pSys->getRefMass();    // nondimensional thrust

	double v = sqrt(s[3]*s[3] + s[4]*s[4]);
	law[0] = sign*(f/s[6])*s[4]/v;
	law[1] = -sign*(f/s[6])*s[3]/v;
	law[2] = 0;

	(void) pSys;
	(void) t;
}//====================================================

/**
 *  \brief Retrieve the output of the Jacobi-changing, parallel velocity control laws
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param law empty, initialized array to store the control law output in
 *  \param len number of elements in the <tt>law</tt> array
 *  \param sign specifies which of the control laws to evaluate: +1 for with-velocity, -1 for anti-velocity
 */
void ControlLaw_cr3bp_lt::getAccel_Along_Vel(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len, int sign) const{

	if(std::abs(sign) != 1)
		sign = sign/std::abs(sign);	// +1 for With-Velocity, -1 for Anti-Velocity

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_Pro_Vel: law data length must be at least 3!");

	double charT = pSys->getCharT();
	double charL = pSys->getCharL();
	double f = (getThrust()/1000)*charT*charT/charL/pSys->getRefMass();    // nondimensional thrust

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = sign*(f/s[6])*s[3]/v;
	law[1] = sign*(f/s[6])*s[4]/v;
	law[2] = sign*(f/s[6])*s[5]/v;

	(void) pSys;
	(void) t;
}//====================================================

void ControlLaw_cr3bp_lt::getAccel_GeneralDir(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	(void) t;

	if(len < 3)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_GeneralDir: law data length must be at least 3!");

	double charT = pSys->getCharT();
	double charL = pSys->getCharL();
	double f = (getThrust()/1000)*charT*charT/charL/pSys->getRefMass();    // nondimensional thrust

	// Direction is stored in the state variables after the core states
	unsigned int core_dim = pSys->getDynamicsModel()->getCoreStateSize();
	law[0] = (f/s[6])*s[core_dim+0];
	law[1] = (f/s[6])*s[core_dim+1];
	law[2] = (f/s[6])*s[core_dim+2];
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Partial Derivatives of Control Laws
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the partial derivatives of the control law with respect to state variables
 *  \details A set of partial derivatives of the control law outputs are computed with respect to the 
 *  states at the given time, state, in the specified system
 * 
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param partials empty, initialized array to store the control law derivatives in
 *  \param len number of elements in the <tt>law</tt> array
 *  \param sign specifies which of the ConstC_2D control laws to evaluate: +1 for RIGHT, -1 for LEFT
 */
void ControlLaw_cr3bp_lt::getAccelPartials_ConstC_2D(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len, int sign) const{

	if(std::abs(sign) != 1)
		sign = sign/std::abs(sign);	// +1 for RIGHT, -1 for LEFT

	if(len != 21)
		throw Exception("ControlLaw_cr3bp_lt::getPartials_State_ConstC_2D: Expects len = 21");

	// state s : [x, y, z, vx, vy, vz, m, ... stm_elements ...]
	// partials: row 1 = partials of a_x w.r.t. states, row 2 = partials of a_y w.r.t. states,
	// row 3 = partials of a_z w.r.t. states

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	double charT = pSys->getCharT();
	double charL = pSys->getCharL();
	double f = (getThrust()/1000)*charT*charT/charL/pSys->getRefMass();

	// Initialize all partials to zero
	for(unsigned int i = 0; i < len; i++){ partials[i] = 0; }

	partials[7*0 + 3] = -sign*f*s[3]*s[4]/(s[6]*pow(v,3));							// dax/dvx
	partials[7*0 + 4] = sign*(f/(s[6]*v) - f*s[4]*s[4]/(s[6]*pow(v,3)));			// dax/dvy
	partials[7*0 + 6] = -sign*f*s[4]/(s[6]*s[6]*v);									// dax/dm

	partials[7*1 + 3] = -sign*(f/(s[6]*v) - f*s[3]*s[3]/(s[6]*pow(v,3)));			// day/dvx
	partials[7*1 + 4] = sign*f*s[3]*s[4]/(s[6]*pow(v,3));							// day/dvy
	partials[7*1 + 6] = sign*f*s[3]/(s[6]*s[6]*v);									// day/dm

	(void) t;
}//====================================================

void ControlLaw_cr3bp_lt::getAccelPartials_GeneralDir(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len) const{

	if(len != 21)
		throw Exception("ControlLaw_cr3bp_lt::getPartials_State_GeneralDir: Expects len = 21");

	// State s : [x, y, z, vx, vy, vz, m, ux, uy, uz, stm_elements]
	// partials: row 1 = partials of ax w.r.t. states, row 2 = partials of ay w.r.t. states, etc.

	double charT = pSys->getCharT();
	double charL = pSys->getCharL();
	double f = (getThrust()/1000)*charT*charT/charL/pSys->getRefMass();

	// Initialize all partials to zero
	for(unsigned int i = 0; i < len; i++){ partials[i] = 0; }

	partials[7*0 + 6] = -f*s[7]/(s[6]*s[6]);	// dax/dm
	partials[7*1 + 6] = -f*s[8]/(s[6]*s[6]);	// day/dm
	partials[7*2 + 6] = -f*s[9]/(s[6]*s[6]);	// daz/dm

	(void) t;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

void ControlLaw_cr3bp_lt::init(){
	switch(lawID){
		case Law_tp::CONST_C_2D_LEFT:
		case Law_tp::CONST_C_2D_RIGHT:
		case Law_tp::PRO_VEL:
		case Law_tp::ANTI_VEL:
			numStates = 0;	// all directions are functions of other state variables; no need for new ones
			break;
		case Law_tp::GENERAL_CONST_F:
			numStates = 3;	// Three directions are free
			break;
		default:
			ControlLaw::init();
	}
}//====================================================

/**
 *  \brief Retrieve a string that represents the law ID
 * 
 *  \param id control law ID
 *  \return a string that represents the law ID
 */
std::string ControlLaw_cr3bp_lt::lawIDToString(unsigned int id) const{
	switch(id){
		case Law_tp::CONST_C_2D_LEFT: return "Jacobi-Preserving, 2D, Left";
		case Law_tp::CONST_C_2D_RIGHT: return "Jacobi-Preserving, 2D, Right";
		case Law_tp::PRO_VEL: return "Prograde Velocity";
		case Law_tp::ANTI_VEL: return "Anti-Velocity";
		case Law_tp::GENERAL_CONST_F: return "General Direction, Const. Thrust";
		default:
			return ControlLaw::lawIDToString(id);
	}
}//====================================================

}
