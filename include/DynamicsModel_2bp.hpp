/**
 *  \file DynamicsModel_2bp.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version August 24, 2016
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

#include "DynamicsModel.hpp"

namespace astrohelion{

// Forward declarations
class SysData_2bp;

/**
 *	\ingroup model 2bp
 *	\brief Derivative of DynamicsModel, specific to the 2BP
 *
 *	The base class's methods provide a good framework for this system,
 *	so only minimal adjustments are needed.
 */
class DynamicsModel_2bp : public DynamicsModel{
public:

	/**
	 *  \name *structors
	 *  \{
	 */
	DynamicsModel_2bp();
	DynamicsModel_2bp(const DynamicsModel_2bp&);
	~DynamicsModel_2bp() {}
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	DynamicsModel_2bp& operator=(const DynamicsModel_2bp&);
	//\}

	/**
	 *  \name Core Analysis Functions
	 *  \{
	 */
	DynamicsModel::eom_fcn getFullEOM_fcn() const;
	DynamicsModel::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const SysData*) const;
	void getPrimPos(double, const SysData*, int, double*) const;
	std::vector<double> getPrimVel(double, const SysData*) const;
	void getPrimVel(double, const SysData*, int, double*) const;
	std::vector<double> getStateDeriv(double, std::vector<double>, EOM_ParamStruct*) const;
	//\}

	/**
	 *  \name Multiple Shooting Analysis Functions
	 *  \{
	 */
	void multShoot_initIterData(MultShootData *it) const;
	//\}
	
	/**
	 * \name Simulation Analysis Functions
	 */
		// Use the default Base class functions
	//\}

	/**
	 *  \name Equations of Motion
	 *  \{
	 */
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getUDDots(double, double, double, double, double*);
	//\}

};

}// END of Astrohelion namespace