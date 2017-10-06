/**
 *  \file DynamicsModel_cr3bp_lt.hpp
 *	\brief Dynamical model for a combined low-thrust CR3BP environment
 *	
 *	\author Andrew Cox
 *	\version March 3, 2017
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

#include "DynamicsModel_cr3bp.hpp"

namespace astrohelion{

// Forward Declarations
class SysData_cr3bp_lt;

/**
 *	\ingroup model cr3bp_lt
 *	\brief Derivative of DynamicsModel, specific to the low-thrust CR3BP
 */
class DynamicsModel_cr3bp_lt : public DynamicsModel_cr3bp{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	DynamicsModel_cr3bp_lt();
	DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt&);
	~DynamicsModel_cr3bp_lt() {}
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	DynamicsModel_cr3bp_lt& operator=(const DynamicsModel_cr3bp_lt&);
	//\}

	/**
	 *  \name Core Analysis Functions
	 *  \{
	 */
	static void getEquilibPt(const SysData_cr3bp_lt*, int, double, double, std::vector<double>*, Verbosity_tp verb = Verbosity_tp::NO_MSG);
	DynamicsModel::eom_fcn getFullEOM_fcn() const;
	DynamicsModel::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getStateDeriv(double, std::vector<double>, EOM_ParamStruct*) const;
	//\}

	/**
	 *  \name Equations of Motion
	 *  \{
	 */
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	//\}
	
	/**
	 *  \name Simulation Analysis Functions
	 *  \{
	 */
	std::vector<Event> sim_makeDefaultEvents(const SysData *pSys) const;
	//\}

	/**
	 *  \name Multiple Shooting Analysis Functions
	 *  \{
	 */
	void multShoot_initIterData(MultShootData *it) const override;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	ControlLaw* createControlLaw() const override;
	bool supportsControl(const ControlLaw*) const override;
	//\}
};

}