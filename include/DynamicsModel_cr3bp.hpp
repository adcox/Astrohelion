/**
 *  @file DynamicsModel_cr3bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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
class SysData_cr3bp;

/**
 *	@ingroup model cr3bp
 *	@brief Derivative of DynamicsModel, specific to the CR3BP
 *
 *	The base class's methods provide a good framework for this system,
 *	so only minimal adjustments are needed. This model adds support for 
 *	Jacobi-Constant targeting.
 */
class DynamicsModel_cr3bp : public DynamicsModel{
public:
	/**
	 *  @name *structors
	 *  @{
	 */
	DynamicsModel_cr3bp();
	DynamicsModel_cr3bp(const DynamicsModel_cr3bp&);
	~DynamicsModel_cr3bp() {}
	//@}

	DynamicsModel_cr3bp& operator=(const DynamicsModel_cr3bp&);

	/**
	 *  @name Core Functions
	 *  @{
	 */
	DynamicsModel::eom_fcn getFullEOM_fcn() const;
	DynamicsModel::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const SysData*) const;
	std::vector<double> getPrimVel(double, const SysData*) const;
	//@}

	/**
	 *  @name Static Calculations
	 *  @{
	 */
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getEquilibPt(const SysData_cr3bp*, int, double, double[3]);
	static double getJacobi(const double[], double);
	static void getUDDots(double, double, double, double, double* ddots);
	//@}

	/**
	 *  @name Simulation Support Functions
	 *  @{
	 */
	void sim_saveIntegratedData(const double*, double, Traj*) const;
	bool sim_locateEvent(Event, Traj*, const double*, double, double, Verbosity_tp) const;
	//@}

	/**
	 *  @name Multiple Shooting Support Functions
	 *  @{
	 */
	void multShoot_applyConstraint(MultShootData*, Constraint, int) const override;
	void multShoot_createOutput(const MultShootData*, const Nodeset*, bool, Nodeset*) const;
	void multShoot_initIterData(MultShootData *it) const override;
	//@}
protected:

	/**
	 *  @name Multiple Shooting Support Functions
	 *  @{
	 */
	void multShoot_targetJC(MultShootData*, Constraint, int) const;
	void multShoot_targetPseudoArc(MultShootData*, Constraint, int) const;
	//@}
};

}// END of Astrohelion namespace