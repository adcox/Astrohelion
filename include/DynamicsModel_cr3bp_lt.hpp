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

#include "DynamicsModel.hpp"

namespace astrohelion{
/**
 *	\ingroup model cr3bp_lt
 *	\brief Derivative of DynamicsModel, specific to the low-thrust CR3BP
 */
class DynamicsModel_cr3bp_lt : public DynamicsModel{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	DynamicsModel_cr3bp_lt();
	DynamicsModel_cr3bp_lt(const DynamicsModel_cr3bp_lt&);
	~DynamicsModel_cr3bp_lt() {}
	//\}

	DynamicsModel_cr3bp_lt& operator=(const DynamicsModel_cr3bp_lt&);

	/**
	 *  \name Core Functions
	 *  \{
	 */
	DynamicsModel::eom_fcn getFullEOM_fcn() const;
	DynamicsModel::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const SysData*) const;
	std::vector<double> getPrimVel(double, const SysData*) const;
	std::vector<double> getStateDeriv(double, std::vector<double>, EOM_ParamStruct*) const;
	//\}

	/**
	 *  \name Static Calculations
	 *  \{
	 */
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static double getJacobi(const double[], double);
	//\}
	
	/**
	 *  \name Simulation Support Functions
	 *  \{
	 */
	bool sim_locateEvent(Event, Arcset*, const double*, double, double, EOM_ParamStruct*, Verbosity_tp) const;
	std::vector<Event> sim_makeDefaultEvents(const SysData *pSys) const;
	//\}

	/**
	 *  \name Multiple Shooting Support Functions
	 *  \{
	 */
	void multShoot_createOutput(const MultShootData*, const Arcset*, bool, Arcset*) const;
	void multShoot_initIterData(MultShootData *it) const override;
	//\}
};

}