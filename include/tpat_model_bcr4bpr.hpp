/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef H_MODEL_BCR4BPR
#define H_MODEL_BCR4BPR

#include "tpat_model.hpp"

// Forward declarations
class tpat_sys_data_bcr4bpr;

/**
 *	@brief A derivative of the tpat_model class, specific to the BCR4BPR
 *
 *	This dynamic model overrides many of the base class's functions to add
 *	support for epoch-dependencies present in this non-autonomous system.
 *	It also adds support for SP targeting, which is only currently supported 
 *	in the BCR4BP.
 */
class tpat_model_bcr4bpr : public tpat_model{
public:
	tpat_model_bcr4bpr();
	tpat_model_bcr4bpr(const tpat_model_bcr4bpr&);
	~tpat_model_bcr4bpr() {}

	tpat_model_bcr4bpr& operator=(const tpat_model_bcr4bpr&);

	// Core Functions
	tpat_model::eom_fcn getFullEOM_fcn() const;
	tpat_model::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const tpat_sys_data*) const;
	std::vector<double> getPrimVel(double, const tpat_sys_data*) const;
	void sim_saveIntegratedData(const double*, double, tpat_traj*) const;
	bool sim_locateEvent(tpat_event, tpat_traj*, const double*, double, double, tpat_verbosity_tp) const;

	// Static Calculation Functions
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getPrimaryPos(double, const tpat_sys_data_bcr4bpr*, double*);
	static void getPrimaryVel(double, const tpat_sys_data_bcr4bpr*, double*);
	static void getPrimaryAccel(double, const tpat_sys_data_bcr4bpr*, double*);
	static void orientAtEpoch(double, tpat_sys_data_bcr4bpr*);
	
	// Multiple Shooting functions
	void multShoot_initDesignVec(iterationData*, const tpat_nodeset*) const;
	void multShoot_initIterData(iterationData *it) const;
	void multShoot_scaleDesignVec(iterationData*, const tpat_nodeset*) const;
	void multShoot_createContCons(iterationData*, const tpat_nodeset*) const;
	void multShoot_getSimICs(const iterationData*, const tpat_nodeset*, int, double*, double*, double*) const;
	double multShoot_getSlackVarVal(const iterationData*, tpat_constraint) const;
	void multShoot_applyConstraint(iterationData*, tpat_constraint, int) const;
	void multShoot_createOutput(const iterationData*, const tpat_nodeset*, bool, tpat_nodeset*) const;

protected:
	void multShoot_targetPosVelCons(iterationData*, tpat_constraint, int) const;
	void multShoot_targetExContCons(iterationData*, tpat_constraint, int) const;
	void multShoot_targetState(iterationData*, tpat_constraint, int) const;
	void multShoot_targetDeltaV(iterationData*, tpat_constraint, int) const;
	void multShoot_targetDist(iterationData*, tpat_constraint, int) const;
	double multShoot_targetDist_compSlackVar(const iterationData*, tpat_constraint) const;
	void multShoot_targetApse(iterationData*, tpat_constraint, int) const;
	void multShoot_targetSP(iterationData*, tpat_constraint, int) const;
	void multShoot_targetSP_mag(iterationData*, tpat_constraint, int) const;
	void multShoot_targetSP_dist(iterationData*, tpat_constraint, int) const;
	double multShoot_targetSPMag_compSlackVar(const iterationData*, tpat_constraint) const;
	double multShoot_targetSP_maxDist_compSlackVar(const iterationData*, tpat_constraint) const;
};

#endif