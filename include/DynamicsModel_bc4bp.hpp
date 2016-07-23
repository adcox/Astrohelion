/**
 *  @file DynamicsModel_bc4bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
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
class SysData_bc4bp;

/**
 *	@brief A derivative of the DynamicsModel class, specific to the BCR4BPR
 *
 *	This dynamic model overrides many of the base class's functions to add
 *	support for epoch-dependencies present in this non-autonomous system.
 *	It also adds support for SP targeting, which is only currently supported 
 *	in the BCR4BP.
 */
class DynamicsModel_bc4bp : public DynamicsModel{
public:
	DynamicsModel_bc4bp();
	DynamicsModel_bc4bp(const DynamicsModel_bc4bp&);
	~DynamicsModel_bc4bp() {}

	DynamicsModel_bc4bp& operator=(const DynamicsModel_bc4bp&);

	// Core Functions
	DynamicsModel::eom_fcn getFullEOM_fcn() const;
	DynamicsModel::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const SysData*) const;
	std::vector<double> getPrimVel(double, const SysData*) const;
	void sim_saveIntegratedData(const double*, double, Traj*) const;
	bool sim_locateEvent(Event, Traj*, const double*, double, double, Verbosity_tp) const;

	// Static Calculation Functions
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getPrimaryPos(double, const SysData_bc4bp*, double*);
	static void getPrimaryVel(double, const SysData_bc4bp*, double*);
	static void getPrimaryAccel(double, const SysData_bc4bp*, double*);
	static void orientAtEpoch(double, SysData_bc4bp*);
	
	// Multiple Shooting functions
	void multShoot_initDesignVec(MultShootData*, const Nodeset*) const;
	void multShoot_initIterData(MultShootData *it) const;
	void multShoot_scaleDesignVec(MultShootData*, const Nodeset*) const;
	void multShoot_createContCons(MultShootData*, const Nodeset*) const;
	void multShoot_getSimICs(const MultShootData*, const Nodeset*, int, double*, double*, double*) const;
	double multShoot_getSlackVarVal(const MultShootData*, Constraint) const;
	void multShoot_applyConstraint(MultShootData*, Constraint, int) const;
	void multShoot_createOutput(const MultShootData*, const Nodeset*, bool, Nodeset*) const;

protected:
	void multShoot_targetCont_PosVel(MultShootData*, Constraint, int) const;
	void multShoot_targetCont_Ex(MultShootData*, Constraint, int) const;
	void multShoot_targetCont_Ex_Seg(MultShootData*, Constraint, int) const;
	void multShoot_targetState(MultShootData*, Constraint, int) const;
	void multShoot_targetDeltaV(MultShootData*, Constraint, int) const;
	void multShoot_targetDist(MultShootData*, Constraint, int) const;
	double multShoot_targetDist_compSlackVar(const MultShootData*, Constraint) const;
	void multShoot_targetApse(MultShootData*, Constraint, int) const;
	void multShoot_targetSP(MultShootData*, Constraint, int) const;
	void multShoot_targetSP_mag(MultShootData*, Constraint, int) const;
	void multShoot_targetSP_dist(MultShootData*, Constraint, int) const;
	double multShoot_targetSPMag_compSlackVar(const MultShootData*, Constraint) const;
	double multShoot_targetSP_maxDist_compSlackVar(const MultShootData*, Constraint) const;
};

}