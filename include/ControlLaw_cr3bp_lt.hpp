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

#pragma once

#include "ControlLaw.hpp"

namespace astrohelion{

// Forward declarations
class SysData;
class SysData_cr3bp_lt;

/**
 *  \ingroup model cr3bp_lt
 *  \brief Includes control laws specific to the CR3BP-LT problem
 *  
 *  <code>params</code> layout for various control laws:
 *  * CONST_C_2D_LEFT, CONST_C_2D_RIGHT, PRO_VEL, ANTI_VEL - params holds:
 *  	[0] Thrust (Newtons)
 *  	[1] Isp (Seconds)
 *  * GENERAL_CONST_F - params holds:
 *  	[0] Thrust (Newtons)
 *  	[1] Isp (Seconds)
 *  
 */
class ControlLaw_cr3bp_lt : public ControlLaw{
public:
	/**
	 *  \name Constructors
	 *  \{
	 */
	ControlLaw_cr3bp_lt(unsigned int id = NO_CTRL, std::vector<double> params = {});
	ControlLaw_cr3bp_lt(unsigned int, double, double);

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	double getThrust() const;
	double getIsp() const;

	void setThrust(double);
	void setIsp(double);
	//\}

	/**
	 *  \name Dynamics Functions
	 *  \{
	 */
	void getLaw_EOMPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	void getLaw_Output(double t, const double *s, const SysData *sysData, double *law, unsigned int len) const override;
	void getLaw_OutputPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	std::string lawIDToString(unsigned int) const;
	//\}

	/**
	 *  \brief Identify the control law
	 */
	enum Law_tp : unsigned int{
		CONST_C_2D_LEFT = 1,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control,
									 * thrust left w.r.t. velocity direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 */
		CONST_C_2D_RIGHT = 2,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control, 
									 * thrust right w.r.t. velocity direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 */
		PRO_VEL = 3,				/*!< Thrust along velocity vector (maximum energy increase)
									 * Not yet implemented
									 */
		ANTI_VEL = 4,				/*!< Thrust along anti-velocity vector (maximum energy decrease)
									 * Not yet implemented
									 */
		GENERAL_CONST_F = 5			/*!< Thrust in an arbitrary direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 */
	};
protected:

	void init() override;
	void getAccel_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccel_Along_Vel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccel_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getAccelPartials_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccelPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getEOMPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
};

}// End of astrohelion namespace