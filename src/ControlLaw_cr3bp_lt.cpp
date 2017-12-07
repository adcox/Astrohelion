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

#include "Arcset_cr3bp_lt.hpp"
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
 *  are generally thrust (nondimensional) and Specific Impulse (in seconds) for
 *  most control law implementations
 */
ControlLaw_cr3bp_lt::ControlLaw_cr3bp_lt(unsigned int id, std::vector<double> params) : ControlLaw(id, params){
	init();
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Set and Get Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Retrieve the name of the law type as a string
 *  \return the name of the law type as a string
 */
std::string ControlLaw_cr3bp_lt::getLawTypeString() const{ return ControlLaw_cr3bp_lt::lawTypeToString(lawType); }

//------------------------------------------------------------------------------------------------------
//      Analysis Functions
//------------------------------------------------------------------------------------------------------

double ControlLaw_cr3bp_lt::get_dmdt(double t, const double *s, const SysData *pSys) const{
	switch(lawType){
		case Law_tp::CONST_C_2D_LEFT:
		case Law_tp::CONST_C_2D_RIGHT:
		case Law_tp::PRO_VEL:
		case Law_tp::ANTI_VEL:
		case Law_tp::GENERAL_CONST_F:
			// These laws all store thrust magnitude and Isp as cosntant parameters
			// params : [f, Isp]
			// nondimensional mass flow rate; simplified a bit by cancelling some of the constants
			return -params[0]*pSys->getCharL()/(pSys->getCharT()*params[1]*G_GRAV_0);
		default:
			return 0;
	}
}//====================================================

/**
 *  \brief Retrieve the output of the control law
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data. These outputs are:
 * 	
 * 		[ a_x, a_y, a_z ]
 * 		
 * 	where a = [a_x, a_y, a_z] is the nondimensional acceleration vector that results from 
 * 	the low-thrust propulsion.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSysData system data object
 *  \param law empty, initialized array to store the control law output in. Must be at least
 *  four elements to avoid memory errors.
 *  \param len the number of elements in the <tt>law</tt> array
 */
void ControlLaw_cr3bp_lt::getLaw_Output(double t, const double *s, const SysData *pSysData,
	double *law, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSysData);

	switch(lawType){
		case Law_tp::CONST_C_2D_LEFT:
			getAccel_ConstC_2D(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getAccel_ConstC_2D(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::PRO_VEL:
			getAccel_AlongVel(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::ANTI_VEL:
			getAccel_AlongVel(t, s, pSysData_lt, law, len);
			break;
		case Law_tp::GENERAL_CONST_F:
			getAccel_GeneralDir(t, s, pSysData_lt, law, len);
			break;
		default:
			ControlLaw::getLaw_Output(t, s, pSysData, law, len);
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
 *  \param partials empty, initialized array to store the control law derivatives in
 *  \param len number of elements in the <tt>law</tt> array
 */
void ControlLaw_cr3bp_lt::getLaw_OutputPartials(double t, const double *s, const SysData *pSys, 
	double *partials, unsigned int len) const{

	const SysData_cr3bp_lt *pSysData_lt = static_cast<const SysData_cr3bp_lt *>(pSys);
	switch(lawType){
		case Law_tp::CONST_C_2D_LEFT:
			getAccelPartials_ConstC_2D(t, s, pSysData_lt, partials, len);
			break;
		case Law_tp::CONST_C_2D_RIGHT:
			getAccelPartials_ConstC_2D(t, s, pSysData_lt, partials, len);
			break;
		case Law_tp::GENERAL_CONST_F:
			getAccelPartials_GeneralDir(t, s, pSysData_lt, partials, len);
			break;
		case Law_tp::PRO_VEL:
			getAccelPartials_AlongVel(t, s, pSysData_lt, partials, len);
			break;
		case Law_tp::ANTI_VEL:
			getAccelPartials_AlongVel(t, s, pSysData_lt, partials, len);
			break;
		
		default:
			ControlLaw::getLaw_OutputPartials(t, s, pSys, partials, len);
	}
}//====================================================

/**
 *  \brief Retrieve the partial derivatives of the core state EOMs with respect to the control states
 *  \details If a nontrivial set of control states exists, these partial derivatives form the right-hand
 *  block-column of the A matrix for the rows associated with the core spacecraft state EOMs. I.e., these
 *  partial derivatives do not include the partials of the control state derivatives w.r.t. the control states;
 *  those partial derivatives are obtained from getLaw_StateDerivPartials()
 * 
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param partials initialized array of zeros in which to store the partial derivatives
 *  of the control state derivatives
 *  \param len number of elements in the <tt>partials</tt> array
 */
void ControlLaw_cr3bp_lt::getLaw_EOMPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const{
	switch(lawType){
		case Law_tp::GENERAL_CONST_F:
			getEOMPartials_GeneralDir(t, s, static_cast<const SysData_cr3bp_lt *>(pSys), partials, len);
			break;
		default:
			// Other control laws default to the base behavior (all partials are zero)
			ControlLaw::getLaw_EOMPartials(t, s, pSys, partials, len);
	}
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
 */
void ControlLaw_cr3bp_lt::getAccel_ConstC_2D(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(lawType != Law_tp::CONST_C_2D_LEFT && lawType != Law_tp::CONST_C_2D_RIGHT)
		printWarn("ControlLaw_cr3bp_lt::getAccelPartials_ConstC_2D: Law type is not one of the Const-C types");

	// +1 for RIGHT, -1 for LEFT
	int sign = lawType == CONST_C_2D_RIGHT ? 1 : -1;

	if(len < numOutputs)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_ConstC_2D: law data length must be at least 3!");

	// params[0] is nondimensional thrust
	double v = sqrt(s[3]*s[3] + s[4]*s[4]);
	law[0] = sign*(params[0]/s[6])*s[4]/v;
	law[1] = -sign*(params[0]/s[6])*s[3]/v;
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
 */
void ControlLaw_cr3bp_lt::getAccel_AlongVel(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(lawType != Law_tp::PRO_VEL && lawType != Law_tp::ANTI_VEL)
		printWarn("ControlLaw_cr3bp_lt::getAccelPartials_AlongVel: Law type is not one of the parallel-to-velocity types");

	// +1 for PRO, -1 for ANTI
	int sign = lawType == Law_tp::PRO_VEL ? 1 : -1;

	if(len < numOutputs)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_Pro_Vel: law data length must be at least 3!");

	// params[0] is nondimensional thrust
	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
	law[0] = sign*(params[0]/s[6])*s[3]/v;
	law[1] = sign*(params[0]/s[6])*s[4]/v;
	law[2] = sign*(params[0]/s[6])*s[5]/v;

	(void) pSys;
	(void) t;
}//====================================================

/**
 *  \brief Retrieve the output of the general direction, constant thrust control law
 * 	\details A set of outputs are computed according to the specified control law, given
 * 	the input time, state, and system data.
 * 	
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param law empty, initialized array to store the control law output in
 *  \param len number of elements in the <tt>law</tt> array
 */
void ControlLaw_cr3bp_lt::getAccel_GeneralDir(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *law, unsigned int len) const{

	if(lawType != Law_tp::GENERAL_CONST_F)
		printWarn("ControlLaw_cr3bp_lt::getAccel_GeneralDir: Law type is not general direction");

	if(len < numOutputs)
		throw Exception("ControlLaw_cr3bp_lt::getLaw_GeneralDir: law data length must be at least 3!");

	unsigned int core_dim = pSys->getDynamicsModel()->getCoreStateSize();
	double alpha = s[core_dim+0];
	double beta = s[core_dim+1];

	// Direction is stored in the state variables after the core states
	// params[0] is nondimensional thrust
	law[0] = (params[0]/s[6])*cos(beta)*cos(alpha);
	law[1] = (params[0]/s[6])*cos(beta)*sin(alpha);
	law[2] = (params[0]/s[6])*sin(beta);

	(void) pSys;
	(void) t;
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
 */
void ControlLaw_cr3bp_lt::getAccelPartials_ConstC_2D(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len) const{

	if(lawType != Law_tp::CONST_C_2D_LEFT && lawType != Law_tp::CONST_C_2D_RIGHT)
		printWarn("ControlLaw_cr3bp_lt::getAccelPartials_ConstC_2D: Law type is not one of the Const-C types");

	// +1 for RIGHT, -1 for LEFT
	int sign = lawType == CONST_C_2D_RIGHT ? 1 : -1;

	if(len != numOutputs*7)
		throw Exception("ControlLaw_cr3bp_lt::getAccelPartials_ConstC_2D: Expects len = 21");

	// state s : [x, y, z, vx, vy, vz, m, ... stm_elements ...]
	// params : [f, Isp]
	// partials: row 1 = partials of a_x w.r.t. states, row 2 = partials of a_y w.r.t. states,
	// row 3 = partials of a_z w.r.t. states

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);

	partials[7*0 + 3] = -sign*params[0]*s[3]*s[4]/(s[6]*pow(v,3));							// dax/dvx
	partials[7*0 + 4] = sign*(params[0]/(s[6]*v) - params[0]*s[4]*s[4]/(s[6]*pow(v,3)));	// dax/dvy
	partials[7*0 + 6] = -sign*params[0]*s[4]/(s[6]*s[6]*v);									// dax/dm

	partials[7*1 + 3] = -sign*(params[0]/(s[6]*v) - params[0]*s[3]*s[3]/(s[6]*pow(v,3)));	// day/dvx
	partials[7*1 + 4] = sign*params[0]*s[3]*s[4]/(s[6]*pow(v,3));							// day/dvy
	partials[7*1 + 6] = sign*params[0]*s[3]/(s[6]*s[6]*v);									// day/dm

	(void) pSys;
	(void) t;
}//====================================================

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
 */
void ControlLaw_cr3bp_lt::getAccelPartials_AlongVel(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len) const{

	if(lawType != Law_tp::PRO_VEL && lawType != Law_tp::ANTI_VEL)
		printWarn("ControlLaw_cr3bp_lt::getAccelPartials_AlongVel: Law type is not one of the parallel-to-velocity types");

	// +1 for PRO, -1 for ANTI
	int sign = lawType == Law_tp::PRO_VEL ? 1 : -1;

	if(len != numOutputs*7)
		throw Exception("ControlLaw_cr3bp_lt::getAccelPartials_AlongVel: Expects len = 21");

	// state s : [x, y, z, vx, vy, vz, m, ... stm_elements ...]
	// params : [f, Isp]
	// partials: row 1 = partials of a_x w.r.t. states, row 2 = partials of a_y w.r.t. states,
	// row 3 = partials of a_z w.r.t. states

	double v = sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);

	partials[7*0 + 3] =  sign*(params[0]/s[6]) * (1.0/v - s[3]*s[3]/pow(v,3));		// dax/dvx
	partials[7*0 + 4] = -sign*(params[0]/s[6]) * s[3]*s[4]/pow(v,3);				// dax/dvy
	partials[7*0 + 5] = -sign*(params[0]/s[6]) * s[3]*s[5]/pow(v,3);				// dax/dvz
	partials[7*0 + 6] = -sign*(params[0]/s[6]) * s[3]/(v*s[6]);						// dax/dm

	partials[7*1 + 3] = partials[7*0 + 4];									// day/dvz
	partials[7*1 + 4] =  sign*(params[0]/s[6]) * (1.0/v - s[4]*s[4]/pow(v,3));		// day/dvy
	partials[7*1 + 5] = -sign*(params[0]/s[6]) * s[4]*s[5]/pow(v,3);				// day/dvz
	partials[7*1 + 6] = -sign*(params[0]/s[6]) * s[4]/(v*s[6]);						// day/dm

	partials[7*2 + 3] = partials[7*1 + 5];									// daz/dvx
	partials[7*2 + 4] = partials[7*2 + 5];									// daz/dvy
	partials[7*2 + 5] =  sign*(params[0]/s[6]) * (1.0/v - s[5]*s[5]/pow(v,3));		// daz/dvz
	partials[7*2 + 6] = -sign*(params[0]/s[6]) * s[5]/(v*s[6]);						// daz/dm

	(void) pSys;
	(void) t;
}//====================================================

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
 */
void ControlLaw_cr3bp_lt::getAccelPartials_GeneralDir(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len) const{
	
	if(lawType != Law_tp::GENERAL_CONST_F)
		printWarn("ControlLaw_cr3bp_lt::getAccel_GeneralDir: Law type is not general direction");

	if(len != numOutputs*7)	// 7 core states
		throw Exception("ControlLaw_cr3bp_lt::getAccelPartials_GeneralDir: unexpected array length");

	// State s : [x, y, z, vx, vy, vz, m, ux, uy, uz]
	// params : [f, Isp]
	// partials: row 1 = partials of ax w.r.t. states, row 2 = partials of ay w.r.t. states, etc.

	double alpha = s[7];
	double beta = s[8];

	partials[7*0 + 6] = -params[0]*cos(beta)*cos(alpha)/(s[6]*s[6]);	// dax/dm
	partials[7*1 + 6] = -params[0]*cos(beta)*sin(alpha)/(s[6]*s[6]);	// day/dm
	partials[7*2 + 6] = -params[0]*sin(beta)/(s[6]*s[6]);	// daz/dm

	(void) pSys;
	(void) t;
}//====================================================

/**
 *  \brief Retrieve the partial derivatives of the core state EOMs with respect to the control states
 * 
 *  \param t time parameter
 *  \param s state vector
 *  \param pSys system data object
 *  \param partials initialized array of zeros in which to store the partial derivatives
 *  of the control state derivatives
 *  \param len number of elements in the <tt>partials</tt> array
 */
void ControlLaw_cr3bp_lt::getEOMPartials_GeneralDir(double t, const double *s, const SysData_cr3bp_lt *pSys,
	double *partials, unsigned int len) const{

	if(lawType != Law_tp::GENERAL_CONST_F)
		printWarn("ControlLaw_cr3bp_lt::getAccel_GeneralDir: Law type is not general direction");

	if(len != numStates*7)	// 7 core states
		throw Exception("ControlLaw_cr3bp_lt::getEOMPartials_GeneralDir: unexpected array length");

	// State s : [x, y, z, vx, vy, vz, m, ux, uy, uz]
	// params : [f, Isp]
	// partials: row 1 = partials of vx w.r.t. ctrl states, row 2 = partials of vy w.r.t. ctrl states, etc.

	double alpha = s[7];
	double beta = s[8];

	partials[numStates*3 + 0] = -params[0]/s[6] * cos(beta)*sin(alpha);		// partial of xddot w.r.t. alpha
	partials[numStates*3 + 1] = -params[0]/s[6] * sin(beta)*cos(alpha);		// partial of xddot w.r.t. beta
	partials[numStates*4 + 0] = params[0]/s[6] * cos(beta) * cos(alpha);	// partial of yddot w.r.t. alpha
	partials[numStates*4 + 1] = -params[0]/s[6] * sin(beta)*sin(alpha);		// partial of yddot w.r.t. beta
	partials[numStates*5 + 1] = params[0]/s[6] * cos(beta);					// partial of zddot w.r.t. beta

	(void) pSys;
	(void) t;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief Initialize the control law
 *  \details This function sets variables that specify the number of control states and
 *  output states. The base class functionality is called if the law type is not specific
 *  to the this derived class
 */
void ControlLaw_cr3bp_lt::init(){
	switch(lawType){
		case Law_tp::CONST_C_2D_LEFT:
		case Law_tp::CONST_C_2D_RIGHT:
		case Law_tp::PRO_VEL:
		case Law_tp::ANTI_VEL:
			numStates = 0;	// all directions are functions of other state variables; no need for new ones
			numOutputs = 3;
			break;
		case Law_tp::GENERAL_CONST_F:
			numStates = 2;	// Two angles to represent 3D unit vector
			numOutputs = 3;
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
std::string ControlLaw_cr3bp_lt::lawTypeToString(unsigned int id){
	switch(id){
		case Law_tp::CONST_C_2D_LEFT: return "Jacobi-Preserving, 2D, Left";
		case Law_tp::CONST_C_2D_RIGHT: return "Jacobi-Preserving, 2D, Right";
		case Law_tp::PRO_VEL: return "Prograde Velocity";
		case Law_tp::ANTI_VEL: return "Anti-Velocity";
		case Law_tp::GENERAL_CONST_F: return "General Direction, Const. Thrust";
		default:
			return ControlLaw::lawTypeToString(id);
	}
}//====================================================

/**
 *  \brief Convert dimensional thrust to nondimensional thrust in the specified system
 * 
 *  \param F dimensional thrust in Newtons
 *  \param pSys pointer to the CR3BP-LT system object
 * 
 *  \return nondimensional thrust
 */
double ControlLaw_cr3bp_lt::thrust_dim2nondim(double F, SysData_cr3bp_lt *pSys){
	return F*pSys->getCharT()*pSys->getCharT() / (1000*pSys->getCharL()*pSys->getRefMass());
}//====================================================

/**
 *  \brief Convert nondimensional thrust to dimensional thrust in the specified system
 * 
 *  \param f nondimensional thrust value in the specified system
 *  \param pSys pointer to the CR3BP-LT system object
 * 
 *  \return dimensional thrust, in Newtons
 */
double ControlLaw_cr3bp_lt::thrust_nondim2dim(double f, SysData_cr3bp_lt *pSys){
	return f*1000*pSys->getCharL()*pSys->getRefMass() / (pSys->getCharT()*pSys->getCharT());
}//====================================================

void ControlLaw_cr3bp_lt::print() const{
	printf("CR3BP-Low-Thrust Control Law\n  Type = %s\n", lawTypeToString(lawType).c_str());
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
 *  \brief Convert an arcset with an arbitrary set of control laws to leverage
 *  the specified control law on all segments
 *  \details Not all conversions are valid or well-defined. For those that are, 
 *  the input arcset is modified to have the propper number of control state variables
 *  in the Segment state vector and the control law objects and parameters are updated 
 *  on all Nodes and Segments.
 * 
 *  \param pArcset pointer to the arcset to convert
 *  \param pLaw pointer to the control law that should be leveraged
 */
void ControlLaw_cr3bp_lt::convertLaws(Arcset_cr3bp_lt *pArcset, ControlLaw_cr3bp_lt *pLaw){
	if(pLaw == nullptr)
		throw Exception("ControlLaw_cr3bp_lt::convertLaws: Input control law is nullptr");

	switch(pLaw->getLawType()){
		case Law_tp::GENERAL_CONST_F:
			// call function
			convertTo_GeneralConstF(pArcset, pLaw);
			break;
		default:
			throw Exception("ControlLaw_cr3bp_lt::convertLaws: Conversion to the specified law type is not supported.");
	}
}//====================================================

/**
 *  \brief Convert all control data from an arcset to control data for the GENERAL_CONST_F
 *  law
 *  \details This conversion is well-defined from all of the simplified control laws, i.e.,
 *  CONST_C_2D_LEFT, CONST_C_2D_RIGHT, PRO_VEL, and ANTI_VEL because their directions are
 *  easily computed at each instant in time.
 * 
 *  \param pArcset Pointer to the arcset to be converted
 *  \param pNewLaw Pointer to the GENERAL_CONST_F control law
 */
void ControlLaw_cr3bp_lt::convertTo_GeneralConstF(Arcset_cr3bp_lt *pArcset, ControlLaw_cr3bp_lt *pNewLaw){
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
		unsigned int oldLawType = pOldLaw ? pOldLaw->getLawType() : NO_CTRL;

		// Make sure we know how to convert
		bool knownConversion = false;
		switch(oldLawType){
			case to_underlying(Law_tp::CONST_C_2D_LEFT):
			case to_underlying(Law_tp::CONST_C_2D_RIGHT):
			case to_underlying(Law_tp::PRO_VEL):
			case to_underlying(Law_tp::ANTI_VEL):
			case NO_CTRL:
				knownConversion = true;
				break;
			case to_underlying(Law_tp::GENERAL_CONST_F):
				continue;	// Nothing to do for this segment
		}

		if(!knownConversion)
			throw Exception("ControlLaw_cr3bp_lt::convertTo_GeneralConstF: Conversion between input law type and GENERAL_CONST_F is undefined.");

		//--------------------------------------
		// Convert node control law information
		//--------------------------------------
		angles[0] = 0;	// Reset
		angles[1] = 0;
		if(pOldLaw){
			std::vector<double> originState = refOrigin.getState();
			pOldLaw->getLaw_Output(refOrigin.getEpoch(), &(originState.front()), pArcset->getSysData(), lawOutput.data(), 3);
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
		const unsigned int newStateWidth = oldStateWidth - oldCtrlDim - pow(coreDim + oldCtrlDim, 2) + newCtrlDim + pow(coreDim + newCtrlDim, 2);
		std::vector<double> newSegStates;
		newSegStates.reserve(numStates*newStateWidth);

		for(unsigned int i = 0; i < numStates; i++){
			// Copy over the core state
			newSegStates.insert(newSegStates.end(), oldSegStates.begin() + i*oldStateWidth, oldSegStates.begin() + i*oldStateWidth + coreDim);
			// If an old law exists, update law output and compute new states
			// If an old law does not exist, use angles{} computed from node
			if(pOldLaw){
				pOldLaw->getLaw_Output(refSeg.getTimeByIx(i), &(oldSegStates[i*oldStateWidth]), pArcset->getSysData(), lawOutput.data(), 3);
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
				pOldLaw->getLaw_Output(refNode.getEpoch(), &(nodeState.front()), pArcset->getSysData(), lawOutput.data(), 3);
				pointingVecToAngles(lawOutput, &(angles[0]), &(angles[1]));
			}

			refNode.setExtraParamVec(PARAMKEY_CTRL, angles);
		}
	}
}//====================================================

/**
 *  \brief Convert a vector to spherical coordinates
 *  \details If the input vector has zero length, both angles are set to zero
 * 
 *  \param vec A 3D vector
 *  \param inPlane The in-plane angle of the vector; measured from x-axis, rotating about
 *  +z-axis
 *  \param outOfPlane The out-of-plane angle of the vector; measured from xy-plane, sign
 *  matches sign of z
 */
void ControlLaw_cr3bp_lt::pointingVecToAngles(Eigen::Vector3d vec, double *inPlane, double *outOfPlane){
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
