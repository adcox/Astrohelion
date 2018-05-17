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

#include "ControlLaw_cr3bp_lt.hpp"

#include <cmath>

#include "Arcset_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors
//-----------------------------------------------------------------------------

/**
 *  @brief Construct a default CR3BP low-thrust control law object
 *  
 *  @param id Control Law ID
 *  @param params a vector of parameters used by the control law. These 
 *  parameters  are generally sqrt(thrust mag.) (nondimensional) and Specific 
 *  Impulse (in seconds) for most control law implementations
 */
ControlLaw_cr3bp_lt::ControlLaw_cr3bp_lt(unsigned int id,
	std::vector<double> params) : ControlLaw(id, params){

	init();
}//====================================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve the name of the law type as a string
 *  @return the name of the law type as a string
 */
std::string ControlLaw_cr3bp_lt::getTypeString() const{
	return ControlLaw_cr3bp_lt::typeToString(lawType);
}

//-----------------------------------------------------------------------------
//      Analysis Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Compute the time derivative of spacecraft mass
 *  @details Computes the quantity \f$ \dot{m} \f$ for the CR3BP-LT equations of motion.
 * 
 *  @param t nondimensional time
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys pointer to the system data object
 *  @return the time derivative of spacecraft mass
 */
double ControlLaw_cr3bp_lt::get_dmdt(double t, const double *s,
	const SysData *pSys) const{

	switch(lawType & M_MASK){
		case CSI_VAR_M:
		{
			double f;
			unsigned int ix_param_shift = 0;
			switch(lawType & BASE_MASK){
				case GEN_INERT:
					// These laws store a reference angle as the first parameter
					// params: [theta0, ...]
					ix_param_shift++;
					break;
			}

			switch(lawType & F_MASK){
				case CONST_F:
					// These laws all store thrust magnitude and Isp as 
					// constant parameters:
					// params : [..., f, Isp]
					// nondimensional mass flow rate, simplified by cancelling
					// some of the constants
					f = params[ix_param_shift + 0];
					break;
				case VAR_F_BND:
				case VAR_F_UBND:
					// params : {..., fmax, Isp}
					// s = {x, y, z, vx, vy, vz, m, g, ...}
					f = getThrustMag(t, s, pSys);
					break;
				default:
					throw Exception("ControlLaw_cr3bp_lt::get_dmdt: "
						"unrecognized thrust parameterization");
			}

			return -f*pSys->getCharL()/(pSys->getCharT()*params[ix_param_shift + 1]*G_GRAV_0);
		}
		case CONST_M:
		default:
			return 0;
	}
}//====================================================

/**
 * @brief Compute the thrust magnitude in nondimensional units
 * @details This is NOT the same as the magnitude of the acceleration vector
 * 
 * @param t nondimensional time
 * @param s full state vector; contains core state, control states, STM 
 * elements, and extra states
 * @param pSys pointer to the system data object
 * @return the nondimensional thrust magnitude of this control law
 */
double ControlLaw_cr3bp_lt::getThrustMag(double t, const double *s, 
	const SysData *pSys) const{

	(void) t;
	(void) pSys;

	unsigned int ix_param_shift = 0;
	switch(lawType & BASE_MASK){
		case GEN_INERT:
			// This law stores an angle to orient rotating frame as first param
			// params = [theta0, ...]
			ix_param_shift++;
			break;
	}

	double g = 0;
	if((lawType & F_MASK) == CONST_F){
		return params[ix_param_shift];		// Thrust magnitude stored as first param
	}else{
		// Variable-Thrust Laws; First, get the thrust magnitude variable
		switch(lawType & BASE_MASK){
			case GENERAL:
			case GEN_INERT:
				// s = {x, y, z, vx, vy, vz, m, alpha, beta, g, ...}
				g = s[9];
				break;
			case VEL_PT:
			case CONST_C_2D:
				// s = {x, y, z, vx, vy, vz, m, g, ...}
				g = s[7];
				break;
		}

		switch(lawType & F_MASK){	
			case VAR_F_BND: 
				// params = {..., fmax, ...}
				return 0.5*params[ix_param_shift]*(sin(g) + 1);
			case VAR_F_UBND: 
				// params = {..., fmax, ...}
				return params[ix_param_shift]*g*g;
			default:
				throw Exception("ControlLaw_cr3bp_lt::getThrustMag: "
					"Control has unrecognized combo of BASE and F types");
		}
	}
}//====================================================

/**
 *  @brief Retrieve the output of the control law
 * 	@details A set of outputs are computed according to the specified control 
 * 	law, given the input time, state, and system data. These outputs are:
 * 	
 * 		[ a_x, a_y, a_z ]
 * 		
 * 	where a = [a_x, a_y, a_z] is the nondimensional acceleration vector that 
 * 	results from the low-thrust propulsion.
 * 	
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSysData system data object
 *  @param output empty, initialized array to store the control law output in. 
 *  Must be at least three elements to avoid memory errors.
 *  @param len the number of elements in the `output` array
 */
void ControlLaw_cr3bp_lt::getOutput(double t, const double *s, 
	const SysData *pSysData, double *output, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = 
		static_cast<const SysData_cr3bp_lt *>(pSysData);
	assert(pSysData_lt->getDynamicsModel()->getCoreStateSize() == 7);

	switch(lawType & BASE_MASK){
		case CONST_C_2D:
			getAccel_ConstC_2D(t, s, pSysData_lt, output, len);
			break;
		case VEL_PT:
			getAccel_AlongVel(t, s, pSysData_lt, output, len);
			break;
		case GENERAL:
			getAccel_GeneralRot(t, s, pSysData_lt, output, len);
			break;
		case GEN_INERT:
			getAccel_GeneralInert(t, s, pSysData_lt, output, len);
			break;
		default:
			ControlLaw::getOutput(t, s, pSysData, output, len);
	}
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control law output with 
 *  respect to the core state variables
 *  @details A set of partial derivatives of the control law outputs are 
 *  computed with respect to the core states at the given time, state, in the 
 *  specified system
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param partials empty, initialized array of zeros to store the control law 
 *  derivatives in
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_OutputWRTCoreState(double t,
	const double *s, const SysData *pSys, double *partials,
	unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt =
		static_cast<const SysData_cr3bp_lt *>(pSys);
	assert(pSysData_lt->getDynamicsModel()->getCoreStateSize() == 7);

	switch(lawType & BASE_MASK){
		case CONST_C_2D:
			getPartials_AccelWRTCore_ConstC_2D(t, s, pSysData_lt, partials, len);
			break;
		case VEL_PT:
			getPartials_AccelWRTCore_AlongVel(t, s, pSysData_lt, partials, len);
			break;
		case GENERAL:
			getPartials_AccelWRTCore_GeneralRot(t, s, pSysData_lt, partials, len);
			break;
		case GEN_INERT:
			getPartials_AccelWRTCore_GeneralInert(t, s, pSysData_lt, partials, len);
			break;
		default:
			ControlLaw::getPartials_OutputWRTCoreState(t, s, pSysData_lt,
				partials, len);
	}

}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the core state EOMs with
 *  respect to the control states
 *  @details If a nontrivial set of control states exists, these partial 
 *  derivatives form the right-hand block-column of the A matrix for the rows 
 *  associated with the core spacecraft state EOMs. I.e., these partial 
 *  derivatives do not include the partials of the control state derivatives 
 *  w.r.t. the control states; those partial derivatives are obtained from 
 *  getPartials_TimeDerivWRTAllState()
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial 
 *  derivatives
 *  of the control state derivatives
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrlState(double t, 
	const double *s, const SysData *pSys, double *partials,
	unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = 
		static_cast<const SysData_cr3bp_lt *>(pSys);
	assert(pSysData_lt->getDynamicsModel()->getCoreStateSize() == 7);

	switch(lawType & BASE_MASK){
		case GENERAL:
		case GEN_INERT:
			getPartials_EOMsWRTCtrl_GeneralDir(t, s, pSysData_lt, partials, len);
			break;
		default:
			switch(lawType & F_MASK){
				case CONST_F:
					// Other control laws default to the base behavior (all 
					// partials are zero)
					ControlLaw::getPartials_EOMsWRTCtrlState(t, s, pSys, 
						partials, len);
					break;
				default:
					getPartials_EOMsWRTCtrl_VarF(t, s, pSysData_lt, partials, 
						len);
			}
	}
	// if((lawType & BASE_MASK) == GENERAL){
	// 	getPartials_EOMsWRTCtrl_GeneralDir(t, s, pSysData_lt, partials, len);
	// 	return;
	// }

	// if((lawType & F_MASK) != CONST_F){
	// 	getPartials_EOMsWRTCtrl_VarF(t, s, pSysData_lt, partials, len);
	// 	return;
	// }else{
	// 	// Other control laws default to the base behavior (all partials are zero)
	// 	ControlLaw::getPartials_EOMsWRTCtrlState(t, s, pSys, partials, len);
	// }
	
}//====================================================

//-----------------------------------------------------------------------------
//      Control Laws
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve the output of the Jacobi-preserving 2D control laws
 * 	@details A set of outputs are computed according to the specified control 
 * 	law, given the input time, state, and system data.
 * 	
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param law empty, initialized array to store the control law output in
 *  @param len number of elements in the `law` array
 */
void ControlLaw_cr3bp_lt::getAccel_ConstC_2D(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *law, unsigned int len) const{

	double f = getThrustMag(t, s, pSys);
	int sign = (lawType & OP1_MASK) ? 1 : -1;	// +1 for RIGHT, -1 for LEFT

	if(len < numOutputs){
		char msg[64];
		sprintf(msg, "ControlLaw_cr3bp_lt::getAccel_ConstC_2D: "
			"law data length is %u but must be at least %u", len, numOutputs);
		throw Exception (msg);
	}

	double v = sqrt(s[3]*s[3] + s[4]*s[4]);
	law[0] = sign*(f/s[6])*s[4]/v;
	law[1] = -sign*(f/s[6])*s[3]/v;
	law[2] = 0;
}//====================================================

/**
 *  @brief Retrieve the output of the Jacobi-changing, parallel velocity 
 *  control laws
 * 	@details A set of outputs are computed according to the specified control 
 * 	law, given the input time, state, and system data.
 * 	
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param law empty, initialized array to store the control law output in
 *  @param len number of elements in the `law` array
 */
void ControlLaw_cr3bp_lt::getAccel_AlongVel(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *law, unsigned int len) const{

	double f = getThrustMag(t, s, pSys);
	int sign = (lawType & OP1_MASK) ? -1 : 1;	// +1 for PRO, -1 for ANTI

	if(len < numOutputs){
		char msg[64];
		sprintf(msg, "ControlLaw_cr3bp_lt::getAccel_AlongVel: "
			"law data length is %u but must be at least %u", len, numOutputs);
		throw Exception (msg);
	}

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = sign*(f/s[6])*s[3]/v;
	law[1] = sign*(f/s[6])*s[4]/v;
	law[2] = sign*(f/s[6])*s[5]/v;
}//====================================================

/**
 *  @brief Retrieve the output of the general direction (rotating) control law
 * 	@details A set of outputs are computed according to the specified control 
 * 	law, given the input time, state, and system data.
 * 	
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param law empty, initialized array to store the control law output in
 *  @param len number of elements in the `law` array
 */
void ControlLaw_cr3bp_lt::getAccel_GeneralRot(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *law, unsigned int len) const{

	double f = getThrustMag(t, s, pSys);
	// s = {x, y, z, vx, vy, vz, m, alpha, beta, ...}

	if(len < numOutputs){
		char msg[64];
		sprintf(msg, "ControlLaw_cr3bp_lt::getAccel_GeneralRot: "
			"law data length is %u but must be at least %u", len, numOutputs);
		throw Exception (msg);
	}

	law[0] = (f/s[6])*cos(s[8])*cos(s[7]);	// a_x
	law[1] = (f/s[6])*cos(s[8])*sin(s[7]);	// a_y
	law[2] = (f/s[6])*sin(s[8]);			// a_z
}//====================================================

/**
 *  @brief Retrieve the output of the general direction (inertial) control law
 * 	@details A set of outputs are computed according to the specified control 
 * 	law, given the input time, state, and system data.
 * 	
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param law empty, initialized array to store the control law output in
 *  @param len number of elements in the `law` array
 */
void ControlLaw_cr3bp_lt::getAccel_GeneralInert(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *law, unsigned int len) const{

	double f = getThrustMag(t, s, pSys);
	// s = {x, y, z, vx, vy, vz, m, psi, beta, ...}
	// params = [theta0, ...]

	if(len < numOutputs){
		char msg[64];
		sprintf(msg, "ControlLaw_cr3bp_lt::getAccel_GeneralInert: "
			"law data length is %u but must be at least %u", len, numOutputs);
		throw Exception (msg);
	}

	law[0] = (f/s[6])*cos(s[8])*cos(s[7] - params[0] - t);
	law[1] = (f/s[6])*cos(s[8])*sin(s[7] - params[0] - t);
	law[2] = (f/s[6])*sin(s[8]);
}//====================================================

//-----------------------------------------------------------------------------
//      Partial Derivatives of Control Laws
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve the partial derivatives of the control law output with 
 *  respect to state variables
 *  @details A set of partial derivatives of the control law outputs are 
 *  computed with respect to the states at the given time, state, in the 
 *  specified system
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param partials empty, initialized array of zeros to store the control law 
 *  derivatives in. The partials are stored in mxn row-major order array
 *  where each row represents the partials of one of the m control outputs
 *  with respect to each of the n core states
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_AccelWRTCore_ConstC_2D(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *partials, unsigned int len) const{


	if(len != numOutputs*7){
		char msg[128];
		sprintf(msg, "ControlLaw_cr3bp_lt::"
			"getPartials_AccelWRTCore_ConstC_2D: len = %u, expects len = %u",
			len, numOutputs*7);
		throw Exception(msg);
	}
	double f = getThrustMag(t, s, pSys);
	int sign = (lawType & OP1_MASK) ? 1 : -1;	// +1 for RIGHT, -1 for LEFT

	/*	CONST_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {}
	 *		params: {f, Isp}
	 *	VAR_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {g}
	 *		params: {fmax, Isp}
	 *		
	 *	partials:
	 *		row 0 = partials of a_x w.r.t. s
	 *		row 1 = partials of a_y w.r.t. s
	 *		row 2 = partials of a_z w.r.t. s
	 */
	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	double v3 = v*v*v;

	partials[7*0 + 3] = -sign*f*s[3]*s[4]/(s[6]*v3);				// dax/dvx
	partials[7*0 + 4] = sign*(f/(s[6]*v) - f*s[4]*s[4]/(s[6]*v3));	// dax/dvy
	partials[7*0 + 6] = -sign*f*s[4]/(s[6]*s[6]*v);					// dax/dm

	partials[7*1 + 3] = -sign*(f/(s[6]*v) - f*s[3]*s[3]/(s[6]*v3));	// day/dvx
	partials[7*1 + 4] = sign*f*s[3]*s[4]/(s[6]*v3);					// day/dvy
	partials[7*1 + 6] = sign*f*s[3]/(s[6]*s[6]*v);					// day/dm
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control law output with respect to state variables
 *  @details A set of partial derivatives of the control law outputs are computed with respect to the 
 *  states at the given time, state, in the specified system
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 * 	@param partials empty, initialized array of zeros to store the control law 
 *  derivatives in. The partials are stored in mxn row-major order array
 *  where each row represents the partials of one of the m control outputs
 *  with respect to each of the n core states
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_AccelWRTCore_AlongVel(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *partials, unsigned int len) const{

	if(len != numOutputs*7){
		char msg[128];
		sprintf(msg, "ControlLaw_cr3bp_lt::"
			"getPartials_AccelWRTCore_AlongVel: len = %u, expects len = %u",
			len, numOutputs*7);
		throw Exception(msg);
	}

	double f = getThrustMag(t, s, pSys);
	int sign = (lawType & OP1_MASK) ? -1 : 1;	// +1 for PRO, -1 for ANTI
	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	double v3 = v*v*v;

	/*	CONST_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {}
	 *		params: {f, Isp}
	 *	VAR_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {g}
	 *		params: {fmax, Isp}
	 *		
	 *	partials:
	 *		row 0 = partials of a_x w.r.t. s
	 *		row 1 = partials of a_y w.r.t. s
	 *		row 2 = partials of a_z w.r.t. s
	 */
	partials[7*0 + 3] =  sign*(f/s[6]) * (1.0/v - s[3]*s[3]/v3);	// dax/dvx
	partials[7*0 + 4] = -sign*(f/s[6]) * s[3]*s[4]/v3;				// dax/dvy
	partials[7*0 + 5] = -sign*(f/s[6]) * s[3]*s[5]/v3;				// dax/dvz
	partials[7*0 + 6] = -sign*(f/s[6]) * s[3]/(v*s[6]);				// dax/dm

	partials[7*1 + 3] = partials[7*0 + 4];							// day/dvz
	partials[7*1 + 4] =  sign*(f/s[6]) * (1.0/v - s[4]*s[4]/v3);	// day/dvy
	partials[7*1 + 5] = -sign*(f/s[6]) * s[4]*s[5]/v3;				// day/dvz
	partials[7*1 + 6] = -sign*(f/s[6]) * s[4]/(v*s[6]);				// day/dm

	partials[7*2 + 3] = partials[7*1 + 5];							// daz/dvx
	partials[7*2 + 4] = partials[7*2 + 5];							// daz/dvy
	partials[7*2 + 5] =  sign*(f/s[6]) * (1.0/v - s[5]*s[5]/v3);	// daz/dvz
	partials[7*2 + 6] = -sign*(f/s[6]) * s[5]/(v*s[6]);				// daz/dm
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control law output with 
 *  respect to state variables
 *  @details A set of partial derivatives of the control law outputs are 
 *  computed with respect to the states at the given time, state, in the 
 *  specified system
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 * 	@param partials empty, initialized array of zeros to store the control law 
 *  derivatives in. The partials are stored in mxn row-major order array
 *  where each row represents the partials of one of the m control outputs
 *  with respect to each of the n core states
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_AccelWRTCore_GeneralRot(double t, 
	const double *s, const SysData_cr3bp_lt *pSys, double *partials,
	unsigned int len) const{

	if(len != numOutputs*7){
		char msg[128];
		sprintf(msg, "ControlLaw_cr3bp_lt::"
			"getPartials_AccelWRTCore_GeneralRot: len = %u, expects len = %u",
			len, numOutputs*7);
		throw Exception(msg);
	}

	double f = getThrustMag(t, s, pSys);

	/*	CONST_F and CONST_MF laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {alpha, beta}
	 *		params: {f, Isp}
	 *		
	 *	VAR_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {alpha, beta, g}
	 *		params: {fmax, Isp}
	 *		
	 *	partials:
	 *		row 0 = partials of a_x w.r.t. s
	 *		row 1 = partials of a_y w.r.t. s
	 *		row 2 = partials of a_z w.r.t. s
	 */
	partials[7*0 + 6] = -f*cos(s[8])*cos(s[7])/(s[6]*s[6]);	// dax/dm
	partials[7*1 + 6] = -f*cos(s[8])*sin(s[7])/(s[6]*s[6]);	// day/dm
	partials[7*2 + 6] = -f*sin(s[8])/(s[6]*s[6]);			// daz/dm
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the control law output with 
 *  respect to state variables
 *  @details A set of partial derivatives of the control law outputs are 
 *  computed with respect to the states at the given time, state, in the 
 *  specified system
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 * 	@param partials empty, initialized array of zeros to store the control law 
 *  derivatives in. The partials are stored in mxn row-major order array
 *  where each row represents the partials of one of the m control outputs
 *  with respect to each of the n core states
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_AccelWRTCore_GeneralInert(double t, 
	const double *s, const SysData_cr3bp_lt *pSys, double *partials,
	unsigned int len) const{

	if(len != numOutputs*7){
		char msg[128];
		sprintf(msg, "ControlLaw_cr3bp_lt::"
			"getPartials_AccelWRTCore_GeneralInert: len = %u, expects len = %u",
			len, numOutputs*7);
		throw Exception(msg);
	}

	double f = getThrustMag(t, s, pSys);

	/*	CONST_F and CONST_MF laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {psi, beta}
	 *		params: {theta0, f, Isp}
	 *		
	 *	VAR_F laws:
	 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
	 *		ctrl: {psi, beta, g}
	 *		params: {theta0, fmax, Isp}
	 *		
	 *	partials:
	 *		row 0 = partials of a_x w.r.t. s
	 *		row 1 = partials of a_y w.r.t. s
	 *		row 2 = partials of a_z w.r.t. s
	 */
	partials[7*0 + 6] = -f*cos(s[8])*cos(s[7] - params[0] - t)/(s[6]*s[6]);	// dax/dm
	partials[7*1 + 6] = -f*cos(s[8])*sin(s[7] - params[0] - t)/(s[6]*s[6]);	// day/dm
	partials[7*2 + 6] = -f*sin(s[8])/(s[6]*s[6]);							// daz/dm
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the core state EOMs with respect
 *  to the control states
 *  @details The resulting partials are valid for the following law types
 *  * Law_tp::VAR_F_CONST_C_2D_LEFT
 *  * Law_tp::VAR_F_CONST_C_2D_RIGHT
 *  * Law_tp::VAR_F_PRO_VEL
 *  * Law_tp::VAR_F_ANTI_VEL
 *  These control laws all share the same control vector, which includes only 
 *  the thrust magnitude, so the partials are very similar. The partial 
 *  derivatives when the VAR_F_GENERAL law is employed are computed in 
 *  `getPartials_EOMsWRTCtrl_GeneralDir()`.
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial 
 *  derivatives. The array is nxm in row-major order. Each row stores the 
 *  partials of one of the n core states with respect to the m ctrl states.
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrl_VarF(double t, const double *s,
	const SysData_cr3bp_lt *pSys, double *partials, unsigned int len) const{

	if( (lawType & BASE_MASK) == CONST_C_2D || (lawType & BASE_MASK) == VEL_PT){

		if(len != numStates*7){
			char msg[128];
			sprintf(msg, "ControlLaw_cr3bp_lt::"
				"getPartials_EOMsWRTCtrl_VarF: len = %u, expects len = %u",
				len, numStates*7);
			throw Exception(msg);
		}

		// Get the output and divide by thrust magnitude to get the partials
		double a_lt[3] = {0};
		getOutput(t, s, pSys, a_lt, 3);
		double f = sqrt(a_lt[0]*a_lt[0] + a_lt[1]*a_lt[1] + a_lt[2]*a_lt[2]);


		/*	VAR_F laws:
		 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
		 *		ctrl: {g}
		 *		params: {fmax, Isp}
		 *		
		 *	partials:
		 *		row 0 = partials of v_x w.r.t. ctrl
		 *		row 1 = partials of v_y w.r.t. ctrl
		 *		row 2 = partials of v_z w.r.t. ctrl
		 *		row 3 = partials of xddot w.r.t. ctrl
		 *		row 4 = partials of yddot w.r.t. ctrl
		 *		row 5 = partials of zddot w.r.t. ctrl
		 *		row 6 = partials of mdot w.r.t. ctrl
		 */

		if((lawType & F_MASK) == VAR_F_BND){
			// a_lt = (1/2)*fmax*(sin(g) + 1)*u
			// params = {..., fmax, Isp}
			double coeff = 0.5*params[0]*cos(s[7]);
			if(f > 0){
				// partials of xddot (3), yddot (4), zddot(5) and mdot (5)
				// w.r.t. g
				partials[numStates*3 + 0] = coeff*a_lt[0]/(f*s[6]);
				partials[numStates*4 + 0] = coeff*a_lt[1]/(f*s[6]);
				partials[numStates*5 + 0] = coeff*a_lt[2]/(f*s[6]);
				partials[numStates*6 + 0] = -1*coeff*pSys->getCharL()/
					(params[1]*G_GRAV_0*pSys->getCharT());
			}
		}else if( (lawType & F_MASK) == VAR_F_UBND){
			// a_lt = fmax*g^2
			// params = {..., fmax, Isp}
			double coeff = std::abs(s[7]) > 0 ? 2/s[7] : 0;

			// partials of xddot (3), yddot (4), zddot(5) and mdot (5) w.r.t. g
			partials[numStates*3 + 0] = coeff*a_lt[0];
			partials[numStates*4 + 0] = coeff*a_lt[1];
			partials[numStates*5 + 0] = coeff*a_lt[2];
			partials[numStates*6 + 0] = -coeff*f*pSys->getCharL()/
				(params[1]*G_GRAV_0*pSys->getCharT());
		}else{
			printErr("Control Law: \n");
			print();
			throw Exception("ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrl_VarF: "
				"thrust parameterization is not supported.");
		}
	}else{
		printErr("Control Law: \n");
		print();
		throw Exception("ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrl_VarF: "
			"Law type does not match this function");
	}
}//====================================================

/**
 *  @brief Retrieve the partial derivatives of the core state EOMs with respect 
 *  to the control states
 * 
 *  @param t time parameter
 *  @param s full state vector; contains core state, control states, STM 
 *  elements, and extra states
 *  @param pSys system data object
 *  @param partials initialized array of zeros in which to store the partial 
 *  derivatives. The array is nxm in row-major order. Each row stores the 
 *  partials of one of the n core states with respect to the m ctrl states.
 *  @param len number of elements in the `partials` array
 */
void ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrl_GeneralDir(double t, 
	const double *s, const SysData_cr3bp_lt *pSys, double *partials, 
	unsigned int len) const{

	if( (lawType & BASE_MASK) == GENERAL || (lawType & BASE_MASK) == GEN_INERT){

		if(len != numStates*7){
			char msg[128];
			sprintf(msg, "ControlLaw_cr3bp_lt::"
				"getPartials_EOMsWRTCtrl_GeneralDir: len = %u, expects len = %u",
				len, numStates*7);
			throw Exception(msg);
		}

		double f = getThrustMag(t, s, pSys);

		// Shift the index by the number of parameters used by the base type
		unsigned int ix_param_shift = 0;
		double alpha = s[7];
		if((lawType & BASE_MASK) == GEN_INERT){
			ix_param_shift++;
			alpha -= params[0] + t;
		}

		/*	CONST_F and CONST_MF laws:
		 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
		 *		ctrl: {alpha, beta}
		 *		params: {..., f, Isp}
		 *		
		 *	VAR_F laws:
		 *		s: {x, y, z, vx, vy, vz, m, ... ctrl ... , ... stm ...}
		 *		ctrl: {alpha, beta, g}
		 *		params: {..., fmax, Isp}
		 *		
		 *	partials:
		 *		row 0 = partials of v_x w.r.t. ctrl
		 *		row 1 = partials of v_y w.r.t. ctrl
		 *		row 2 = partials of v_z w.r.t. ctrl
		 *		row 3 = partials of a_x w.r.t. ctrl
		 *		row 4 = partials of a_y w.r.t. ctrl
		 *		row 5 = partials of a_z w.r.t. ctrl
		 *		row 6 = partialx of mdot w.r.t. ctrl
		 */
		
		// Partials of xddot (3), yddot (4), zddot (5), and mdot (6) w.r.t.
		// alpha or psi (0) and beta (1)
		partials[numStates*3 + 0] = -f/s[6] * cos(s[8])*sin(alpha);
		partials[numStates*3 + 1] = -f/s[6] * sin(s[8])*cos(alpha);
		partials[numStates*4 + 0] = f/s[6] * cos(s[8])*cos(alpha);
		partials[numStates*4 + 1] = -f/s[6] * sin(s[8])*sin(alpha);
		partials[numStates*5 + 1] = f/s[6] * cos(s[8]);

		// Additional partials for variable-thrust paramaterizations
		if( (lawType & F_MASK) != CONST_F){
			double dfdg = 0;
			if((lawType & F_MASK) == VAR_F_BND){
				// a_lt = (1/2)*fmax*(sin(g) + 1)*u
				// params = {..., fmax, Isp}
				dfdg = 0.5*params[ix_param_shift + 0]*cos(s[9]);
				
				// partial of mdot w.r.t. g
				partials[numStates*6 + 2] = -dfdg*pSys->getCharL()/
					(params[ix_param_shift + 1]*G_GRAV_0*pSys->getCharT());
			}else if( (lawType & F_MASK) == VAR_F_UBND){
				// a_lt = fmax*g^2
				// params = {..., fmask, Isp}
				dfdg = 2*params[ix_param_shift + 0]*s[9];

				// partial of mdot w.r.t. g
				partials[numStates*6 + 2] = -dfdg*pSys->getCharL()/
					(params[ix_param_shift + 1]*G_GRAV_0*pSys->getCharT());
			}else{
				printErr("Control Law:\n");
				print();
				throw Exception("ControlLaw_cr3bp_lt::getPartials_EOMsWRTCtrl_VarF: "
					"thrust parameterization is not supported.");
			}

			// partials of xddot (3), yddot (4), and zddot (5) w.r.t. g
			partials[numStates*3 + 2] = dfdg*cos(s[8])*cos(alpha)/s[6];
			partials[numStates*4 + 2] = dfdg*cos(s[8])*sin(alpha)/s[6];
			partials[numStates*5 + 2] = dfdg*sin(s[8])/s[6];
		}
	}else{
		printErr("Control Law:\n");
		print();
		throw Exception("ControlLaw_cr3bp_lt::getAccel_GeneralRot: "
			"Law type is not general direction");
	}
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Initialize the control law
 *  @details This function sets variables that specify the number of control 
 *  states and output states. The base class functionality is called if the 
 *  law type is not specific to the this derived class
 */
void ControlLaw_cr3bp_lt::init(){
	unsigned int numParams = 0;	// constant parameters like Isp, fmax, etc.
	numStates = 0;		// control states like angles, thrust magnitude, etc.
	numOutputs = 3;		// {ax, ay, az}

	switch(lawType & BASE_MASK){
		case GEN_INERT:
			numParams++;		// Parameter for theta0
		case GENERAL:
			numStates += 2;		// (alpha or psi), beta
			break;
		case VEL_PT:
		case CONST_C_2D:
			break;				// no states or parameters yet
		default:
			ControlLaw::init();
			return;
	}

	switch(lawType & F_MASK){
		case CONST_F:
			numParams++;	// Parameter for f
			break;
		case VAR_F_BND:
		case VAR_F_UBND:
			numParams++;	// Parameter for fmax
			numStates++;	// Variable for f
			break;
	}

	switch(lawType & M_MASK){
		case CSI_VAR_M:
			numParams++;	// Parameter for Isp
			break;
	}

	if(params.size() != numParams){
		char msg[128];
		sprintf(msg, "ControlLaw_cr3bp_lt::init: "
			"Expect %d input params, but received %zu", numParams,
			params.size());
		throw Exception(msg);
	}
}//====================================================

/**
 *  @brief Retrieve a string that represents the law ID
 * 
 *  @param id control law ID
 *  @return a string that represents the law ID
 */
std::string ControlLaw_cr3bp_lt::typeToString(unsigned int id){
	std::string out = "";
	switch(id & BASE_MASK){
		case GENERAL: out.append("General Pointing (Rot), "); break;
		case GEN_INERT: out.append("General Pointing (Inert), "); break;
		case VEL_PT:
			out.append("Velocity-Pointing ");
			out.append((id & OP1_MASK) ? "(Anti), " : "(Pro), ");
			break;

		case CONST_C_2D:
			out.append("2D Jacobi-Preserving ");
			out.append((id & OP1_MASK) ? "(Right), " : "(Left), ");
			break;
		default: ControlLaw::typeToString(id);
	}

	switch(id & F_MASK){
		case CONST_F: out.append("Const. Thrust, "); break;
		case VAR_F_BND: out.append("Var. Thrust (BND), "); break;
		case VAR_F_UBND: out.append("Var. Thrust, (UBND), "); break;
	}

	switch(id & M_MASK){
		case CONST_M: out.append("Const. Mass"); break;
		case CSI_VAR_M: out.append("Var. Mass (CSI)"); break;
	}

	return out;
}//====================================================

/**
 *  @brief Convert dimensional thrust to nondimensional thrust in the specified system
 * 
 *  @param F dimensional thrust in Newtons
 *  @param pSys pointer to the CR3BP-LT system object
 * 
 *  @return nondimensional thrust
 */
double ControlLaw_cr3bp_lt::thrust_dim2nondim(double F, SysData_cr3bp_lt *pSys){
	return F*pSys->getCharT()*pSys->getCharT() / 
		(1000*pSys->getCharL()*pSys->getRefMass());
}//====================================================

/**
 *  @brief Convert nondimensional thrust to dimensional thrust in the specified system
 * 
 *  @param f nondimensional thrust value in the specified system
 *  @param pSys pointer to the CR3BP-LT system object
 * 
 *  @return dimensional thrust, in Newtons
 */
double ControlLaw_cr3bp_lt::thrust_nondim2dim(double f, SysData_cr3bp_lt *pSys){
	return f*1000*pSys->getCharL()*pSys->getRefMass() / 
		(pSys->getCharT()*pSys->getCharT());
}//====================================================

void ControlLaw_cr3bp_lt::print() const{
	printf("CR3BP-Low-Thrust Control Law\n  Type = %s\n", 
		typeToString(lawType).c_str());
	printf("  NumStates = %u\n  NumOutputs = %u\n", numStates, numOutputs);
	printf("  params = {");
	for(unsigned int i = 0; i < params.size(); i++){
		printf("%f", params[i]);
		if(i < params.size() - 1)
			printf(",  ");
	}
	printf("}\n");
}//====================================================

/**
 *  @brief Convert an arcset with an arbitrary set of control laws to leverage
 *  the specified control law on all segments
 *  @details Not all conversions are valid or well-defined. For those that are, 
 *  the input arcset is modified to have the propper number of control state variables
 *  in the Segment state vector and the control law objects and parameters are updated 
 *  on all Nodes and Segments.
 * 
 *  @param pArcset pointer to the arcset to convert
 *  @param pLaw pointer to the control law that should be leveraged
 */
void ControlLaw_cr3bp_lt::convertLaws(Arcset_cr3bp_lt *pArcset, 
	ControlLaw_cr3bp_lt *pLaw){
	
	if(pLaw == nullptr)
		throw Exception("ControlLaw_cr3bp_lt::convertLaws: "
			"Input control law is nullptr");

	switch(pLaw->getType()){
		case Law_tp::CONST_F_GENERAL:
			// call function
			convertTo_GeneralConstF(pArcset, pLaw);
			break;
		default:
			throw Exception("ControlLaw_cr3bp_lt::convertLaws: "
				"Conversion to the specified law type is not supported.");
	}
}//====================================================

/**
 *  @brief Convert all control data from an arcset to control data for the CONST_F_GENERAL
 *  law
 *  @details This conversion is well-defined from all of the simplified control laws, i.e.,
 *  CONST_FC_2D_LEFT, CONST_FC_2D_RIGHT, CONST_F_PRO_VEL, and CONST_F_ANTI_VEL because their directions are
 *  easily computed at each instant in time.
 * 
 *  @param pArcset Pointer to the arcset to be converted
 *  @param pNewLaw Pointer to the CONST_F_GENERAL control law
 */
void ControlLaw_cr3bp_lt::convertTo_GeneralConstF(Arcset_cr3bp_lt *pArcset, 
	ControlLaw_cr3bp_lt *pNewLaw){
	// pNewLaw is guaranteed by the calling function, convertLaws(), to be non-nullptr and have the correct law type

	std::vector<double> angles(2,0);	// Store in-plane and out-of-plane angles that describe general law
	std::vector<int> convertedNodes;	// Store IDs of nodes that have been converted
	Eigen::Vector3d lawOutput;			// Vector for law outputs
	const unsigned int coreDim = pArcset->getSysData()->getDynamicsModel()->getCoreStateSize();
	// const unsigned int extraDim = pArcset->getSysData()->getDynamicsModel()->getExtraStateSize();
	const unsigned int newCtrlDim = pNewLaw->getNumStates();

	for(unsigned int s = 0; s < pArcset->getNumSegs(); s++){
		// Use references for more concise code
		Segment &refSeg = pArcset->getSegRefByIx(s);
		Node &refOrigin = pArcset->getNodeRef(refSeg.getOrigin());

		const ControlLaw *pOldLaw = refSeg.getCtrlLaw();
		unsigned int oldLawType = pOldLaw ? pOldLaw->getType() : NO_CTRL;

		// Make sure we know how to convert
		bool knownConversion = false;
		switch(oldLawType){
			case to_underlying(Law_tp::CONST_FC_2D_LEFT):
			case to_underlying(Law_tp::CONST_FC_2D_RIGHT):
			case to_underlying(Law_tp::CONST_F_PRO_VEL):
			case to_underlying(Law_tp::CONST_F_ANTI_VEL):
			case NO_CTRL:
				knownConversion = true;
				break;
			case to_underlying(Law_tp::CONST_F_GENERAL):
				continue;	// Nothing to do for this segment
		}

		if(!knownConversion)
			throw Exception("ControlLaw_cr3bp_lt::convertTo_GeneralConstF: Conversion between input law type and CONST_F_GENERAL is undefined.");

		//--------------------------------------
		// Convert node control law information
		//--------------------------------------
		angles[0] = 0;	// Reset
		angles[1] = 0;
		if(pOldLaw){
			std::vector<double> originState = refOrigin.getState();
			pOldLaw->getOutput(refOrigin.getEpoch(), &(originState.front()), pArcset->getSysData(), lawOutput.data(), 3);
			pointingVecToAngles(lawOutput, &(angles[0]), &(angles[1]));
		}
		
		refOrigin.setExtraParamVec(PARAMKEY_CTRL, angles);
		convertedNodes.push_back(refOrigin.getID());

		//--------------------------------------
		// Convert segment control law information
		//--------------------------------------
		std::vector<double> oldSegStates = refSeg.getStateVector();
		const unsigned int oldStateWidth = refSeg.getStateWidth();
		const unsigned int numStates = oldSegStates.size() / oldStateWidth;
		const unsigned int oldCtrlDim = pOldLaw ? pOldLaw->getNumStates() : 0;

		// state width is the sum of coreDim + ctrlDim + (coreDim + ctrlDim)^2 + extraDim
		// The squared term is the number of elements in the STM
		const unsigned int newStateWidth = oldStateWidth - oldCtrlDim - 
			pow(coreDim + oldCtrlDim, 2) + newCtrlDim + pow(coreDim + newCtrlDim, 2);
		std::vector<double> newSegStates;
		newSegStates.reserve(numStates*newStateWidth);

		for(unsigned int i = 0; i < numStates; i++){
			// Copy over the core state
			newSegStates.insert(newSegStates.end(), oldSegStates.begin() + i*oldStateWidth, oldSegStates.begin() + i*oldStateWidth + coreDim);
			// If an old law exists, update law output and compute new states
			// If an old law does not exist, use angles{} computed from node
			if(pOldLaw){
				pOldLaw->getOutput(refSeg.getTimeByIx(i), &(oldSegStates[i*oldStateWidth]), pArcset->getSysData(), lawOutput.data(), 3);
				pointingVecToAngles(lawOutput, &(angles[0]), &(angles[1]));
			}
			// Add the new control states
			newSegStates.insert(newSegStates.end(), angles.begin(), angles.end());
			
			// Update the STM; new entries related to the control are left as zeros
			MatrixXRd oldSTM = Eigen::Map<MatrixXRd>(&(oldSegStates[0]) + i*oldStateWidth + coreDim + oldCtrlDim, coreDim + oldCtrlDim, coreDim + oldCtrlDim);
			MatrixXRd newSTM = MatrixXRd::Zero(coreDim + newCtrlDim, coreDim + newCtrlDim);
			newSTM.block(0, 0, coreDim, coreDim) = oldSTM.block(0, 0, coreDim, coreDim);	// Copy Core STM

			// Add STM states
			newSegStates.insert(newSegStates.end(), newSTM.data(), newSTM.data()+newSTM.size());

			// Add any remaining extra states
			newSegStates.insert(newSegStates.end(), oldSegStates.begin()+ i*oldStateWidth + coreDim + oldCtrlDim + oldSTM.size(), oldSegStates.begin() + (i+1)*oldStateWidth);
		}

		refSeg.setStateVector(newSegStates);
		refSeg.setStateWidth(newStateWidth);
		refSeg.setCtrlLaw(pNewLaw);
	}

	// Update all STMs
	pArcset->setSTMs_sequence();

	// Check to make sure all nodes have been converted
	for(unsigned int n = 0; n < pArcset->getNumNodes(); n++){
		if(std::find(convertedNodes.begin(), convertedNodes.end(), pArcset->getNodeByIx(n).getID()) == convertedNodes.end()){
			// Node has not been converted
			Node &refNode = pArcset->getNodeRefByIx(n);
			int linkedSegID = refNode.getLink(0) == Linkable::INVALID_ID ? refNode.getLink(1) : refNode.getLink(0);
			ControlLaw *pOldLaw = pArcset->getSegRef(linkedSegID).getCtrlLaw();

			angles[0] = 0;
			angles[1] = 0;
			if(pOldLaw){
				std::vector<double> nodeState = refNode.getState();
				pOldLaw->getOutput(refNode.getEpoch(), &(nodeState.front()), pArcset->getSysData(), lawOutput.data(), 3);
				pointingVecToAngles(lawOutput, &(angles[0]), &(angles[1]));
			}

			refNode.setExtraParamVec(PARAMKEY_CTRL, angles);
		}
	}
}//====================================================

/**
 *  @brief Convert a vector to spherical coordinates
 *  @details If the input vector has zero length, both angles are set to zero
 * 
 *  @param vec A 3D vector
 *  @param inPlane The in-plane angle of the vector; measured from x-axis, rotating about
 *  +z-axis
 *  @param outOfPlane The out-of-plane angle of the vector; measured from xy-plane, sign
 *  matches sign of z
 */
void ControlLaw_cr3bp_lt::pointingVecToAngles(Eigen::Vector3d vec, 
	double *inPlane, double *outOfPlane){

	if(vec.norm() == 0){
		if(inPlane)
			*inPlane = 0;

		if(outOfPlane)
			*outOfPlane = 0;
	}else{
		vec /= vec.norm();
		
		if(inPlane)
			*inPlane = atan2(vec[1], vec[0]);

		if(outOfPlane)
			*outOfPlane = atan2(vec[2], sqrt(vec[0]*vec[0] + vec[1]*vec[1]));
	}
}//====================================================

}// End of astrohelion namespace
