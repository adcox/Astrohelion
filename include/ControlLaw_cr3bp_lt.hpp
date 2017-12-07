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
#include "EigenDefs.hpp"

namespace astrohelion{

// Forward declarations
class Arcset_cr3bp_lt;
class SysData;
class SysData_cr3bp_lt;

/**
 *  \ingroup model cr3bp_lt
 *  \brief Includes control laws specific to the CR3BP-LT problem
 *  
 */
class ControlLaw_cr3bp_lt : public ControlLaw{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	ControlLaw_cr3bp_lt(unsigned int id = NO_CTRL, std::vector<double> params = {});
	// ControlLaw_cr3bp_lt(unsigned int, double, double);

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	std::string getLawTypeString() const override;
	// double getThrust() const;
	// double getThrust_dim(const SysData_cr3bp_lt*) const;
	// double getIsp() const;

	// void setThrust(double);
	// void setThrust_dim(double, const SysData_cr3bp_lt*);
	// void setIsp(double);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void getLaw_EOMPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	void getLaw_Output(double t, const double *s, const SysData *sysData, double *law, unsigned int len) const override;
	void getLaw_OutputPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	double get_dmdt(double t, const double *s, const SysData *pSys) const;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	static std::string lawTypeToString(unsigned int);
	static void convertLaws(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	static double thrust_dim2nondim(double, SysData_cr3bp_lt*);
	static double thrust_nondim2dim(double, SysData_cr3bp_lt*);
	void print() const;
	//\}

	/**
	 *  \brief Identify the control law
	 */
	enum Law_tp : unsigned int{
		CONST_C_2D_LEFT = 1,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control,
									 * thrust left w.r.t. velocity direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 * - params contains: { thrust (nondim), Isp (seconds) }
									 * - No additional states required for integration
									 */
		CONST_C_2D_RIGHT = 2,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control, 
									 * thrust right w.r.t. velocity direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 * - params contains: { thrust (nondim), Isp (seconds) }
									 * - No additional states required for integration
									 */
		PRO_VEL = 3,				/*!< Thrust along velocity vector (maximum energy increase). Thrust
									 * 	magnitude is constant.
									 * - params contains: { thrust (nondim), Isp (seconds) }
									 * - No additional states required for integration
									 */
		ANTI_VEL = 4,				/*!< Thrust along anti-velocity vector (maximum energy decrease). Thrust
									 * 	magnitude is constant.
									 * - params contains: { thrust (nondim), Isp (seconds) }
									 * - No additional states required for integration
									 */
		GENERAL_CONST_F = 5			/*!< Thrust in an arbitrary direction. Thrust magnitude is constant.
									 * - getLaw() returns 3 thrust directions
									 * - getPartials_State() returns 21 derivatives
									 * - params contains: { thrust (nondim), Isp (seconds) }
									 * - Requires two additional states {alpha, beta} for numerical integration
									 */
	};
protected:

	/**
	 *  \name *structors
	 *  \{
	 */
	void init() override;
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void getAccel_AlongVel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccel_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccel_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getAccelPartials_AlongVel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccelPartials_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getAccelPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getEOMPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	static void pointingVecToAngles(Eigen::Vector3d, double*, double*);
	static void convertTo_GeneralConstF(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	//\}
};

}// End of astrohelion namespace