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
class SysData_bc4bp;

/**
 *	\ingroup model bc4bp
 *	@brief A derivative of the DynamicsModel class, specific to the BCR4BPR
 *
 *	This dynamic model overrides many of the base class's functions to add
 *	support for epoch-dependencies present in this non-autonomous system.
 *	It also adds support for SP targeting, which is only currently supported 
 *	in the BCR4BP.
 */
class DynamicsModel_bc4bp : public DynamicsModel{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	DynamicsModel_bc4bp();
	DynamicsModel_bc4bp(const DynamicsModel_bc4bp&);
	~DynamicsModel_bc4bp() {}
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	DynamicsModel_bc4bp& operator=(const DynamicsModel_bc4bp&);
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
	void getPrimAccel(double, const SysData*, int, double*) const;
	std::vector<double> getStateDeriv(double, std::vector<double>, EOM_ParamStruct*) const;
	//\}

	/**
	 *  \name Equations of Motion
	 *  \{
	 */
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void orientAtEpoch(double, SysData_bc4bp*);
	//\}
	
	/**
	 *  \name Simulation Analysis Functions
	 *  \{
	 */
	int sim_addNode(Node&, const double*, double, Arcset*, EOM_ParamStruct*, Event_tp) const;
	//\}

	/**
	 *  \name Multiple Shooting Analysis Functions
	 *  \{
	 */
	void multShoot_initDesignVec(MultShootData&) const override;
	void multShoot_initIterData(MultShootData& it) const override;
	void multShoot_createContCons(MultShootData&) const override;
	void multShoot_getSimICs(const MultShootData&, int, double*, double*, double*, double*) const override;
	double multShoot_getSlackVarVal(const MultShootData&, const Constraint&) const override;
	void multShoot_applyConstraint(MultShootData&, const Constraint&, int) const override;
	void multShoot_createOutput(const MultShootData&) const override;
	//\}

protected:

	/**
	 *  \name Multiple Shooting Analysis Functions
	 *  \{
	 */
	void multShoot_targetApse(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetApse_endSeg(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetCont_State(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetCont_Ex(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetCont_Ex_Seg(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetDeltaV(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetDist(MultShootData&, const Constraint&, int) const override;
	void multShoot_targetDist_endSeg(MultShootData&, const Constraint&, int) const override;
	double multShoot_targetDist_compSlackVar(const MultShootData&, const Constraint&) const override;
	double multShoot_targetDist_endSeg_compSlackVar(const MultShootData&, const Constraint&) const override;
	void multShoot_targetEpoch(MultShootData&, const Constraint&, int) const;
	void multShoot_targetSP(MultShootData&, const Constraint&, int) const;
	void multShoot_targetSP_mag(MultShootData&, const Constraint&, int) const;
	void multShoot_targetSP_dist(MultShootData&, const Constraint&, int) const;
	double multShoot_targetSPMag_compSlackVar(const MultShootData&, const Constraint&) const;
	double multShoot_targetSP_maxDist_compSlackVar(const MultShootData&, const Constraint&) const;
	void multShoot_targetState_endSeg(MultShootData&, const Constraint&, int) const override;
	//\}
};

}