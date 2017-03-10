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
 */
class ControlLaw_cr3bp_lt : public ControlLaw{
public:
	ControlLaw_cr3bp_lt();

	void getLaw(double t, const double *s, const SysData *sysData, unsigned int lawID, double *law, unsigned int len) const;
	void getPartials_State(double t, const double *s, const SysData *pSys, unsigned int lawID, double *partials, unsigned int len) const;

	std::string lawIDToString(unsigned int) const;

	enum Law_tp : unsigned int{
		CONST_C_2D_LEFT = 1,
		CONST_C_2D_RIGHT = 2,
		PRO_VEL = 3,
		ANTI_VEL = 4
	};
protected:

	void getLaw_ConstC_2D_Right(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
	void getLaw_Along_Vel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;

	void getPartials_State_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int, int) const;
};

}// End of astrohelion namespace