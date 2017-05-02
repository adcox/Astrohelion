/**
 *  \file LambertArcEngine.hpp
 *	\brief Generate Lambert Arcs via a Lambert Solver
 *	
 *	\author Andrew Cox
 *	\version September 30, 2016
 *	\copyright GNU GPL v3.0
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

#include "Core.hpp"
#include "Engine.hpp"

#include "Common.hpp"

#include <vector>

namespace astrohelion{

// Forward Declarations
class Arcset_2bp;
class SysData_2bp;

/**
 *	\ingroup engine
 *	\brief Generate Lambert Arcs in the 2BP
 *
 *	\author Andrew Cox
 *	\version September 9, 2016
 *	\copyright GNU GPL v3.0
 */
class LambertArcEngine : public Core, public Engine{

public:

	LambertArcEngine();
	LambertArcEngine(const LambertArcEngine&);

	LambertArcEngine& operator =(const LambertArcEngine&);

	Arcset_2bp getLambertArc(SysData_2bp*, std::vector<double>, std::vector<double>, double, unsigned int);

	void setMaxErr_tof(double);
	void setMaxErr_ta(double);
	void setMaxIts(unsigned int);

private:

	void cleanEngine();
	void reset();
	void copyMe(const LambertArcEngine&);

	/** Maximum acceptable error in time-of-flight, seconds */
	double tof_maxErr = 1e-3;

	/** Maximum acceptable error in true anomaly, rad */
	double ta_maxErr = 0.001*PI/180;

	/** Maximum number of iterations to solve the Lambert problem */
	unsigned int maxIts = 100;



};

}// End of Astrohelion namespace