/**
 * @file ControlLaw.hpp
 * @brief Control Law header file
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

/* TODO
 *	- Write control laws for LowThrust
 *	- Test it out!
 */

#pragma once

#include <string>

#include "Core.hpp"

namespace astrohelion{

// Forward declarations
class SysData;

/**
 *  \ingroup model
 *  \brief Stores control law equations and derivatives for flexible applications throughout
 *  the astrohelion framework
 */
class ControlLaw : public Core{
public:
	ControlLaw();

	virtual void getLaw(double t, const double *s, const SysData *pSys, unsigned int lawID, double *law, unsigned int len) const;
	virtual void getPartials_State(double t, const double *s, const SysData *pSys, unsigned int lawID, double *partials, unsigned int len) const;
	virtual std::string lawIDToString(unsigned int) const;

	const static unsigned int NO_CTRL = 0;
};

}// End of astrohelion namespace