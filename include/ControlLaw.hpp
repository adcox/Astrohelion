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

#pragma once

#include "matio.h"
#include <string>
#include <vector>

#include "Core.hpp"

namespace astrohelion{

// Forward declarations
class SysData;

/**
 *  \ingroup model
 *  @brief Stores control law equations and derivatives for flexible applications throughout
 *  the astrohelion framework
 */
class ControlLaw : public Core{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	ControlLaw(unsigned int id = NO_CTRL, std::vector<double> params = {});
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	ControlLaw& operator =(const ControlLaw&);
	friend bool operator ==(const ControlLaw&, const ControlLaw&);
	friend bool operator !=(const ControlLaw&, const ControlLaw&);
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	unsigned int getLawType() const;
	virtual std::string getLawTypeString() const;
	unsigned int getNumOutputs() const;
	unsigned int getNumStates() const;
	double getParam(int) const;
	std::vector<double> getParams() const;
	const std::vector<double>& getParamsRef_const() const;

	void setLawType(unsigned int);
	void setParam(int, double);
	void setParams(double*, unsigned int);
	void setParams(std::vector<double>);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	virtual void getLaw_EOMPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const;
	virtual void getLaw_Output(double t, const double *s, const SysData *pSys, double *output, unsigned int len) const;
	virtual void getLaw_OutputPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const;
	virtual void getLaw_StateDeriv(double t, const double *s, const SysData *pSys, double *deriv, unsigned int len) const;
	virtual void getLaw_StateDerivPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	static std::string lawTypeToString(unsigned int);
	virtual void print() const;
	//\}
	
	const static unsigned int NO_CTRL = 0;	//!< Value to use for the control law ID when no control law is implemented

protected:
	/**
	 *  \name Utility Functions
	 *  \{
	 */
	virtual void init();
	void copyMe(const ControlLaw&);
	//\}

	unsigned int lawType = NO_CTRL;		//!< Value identifying the specific control forumalation to apply
	unsigned int numStates = 0;			//!< Number of control states
	unsigned int numOutputs = 0;		//!< Number of control outputs

	std::vector<double> params {};		//!< Parameters associated with the control law
};

}// End of astrohelion namespace