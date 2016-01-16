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
	tpat_model::eom_fcn getFullEOM_fcn();
	tpat_model::eom_fcn getSimpleEOM_fcn();
	std::vector<double> getPrimPos(double, tpat_sys_data*);
	std::vector<double> getPrimVel(double, tpat_sys_data*);
	void sim_saveIntegratedData(double*, double, tpat_traj*);
	bool sim_locateEvent(tpat_event, tpat_traj*, double*, double, double, verbosity_t);

	// Multiple Shooting functions
	void multShoot_initDesignVec(iterationData*, tpat_nodeset*);
	void multShoot_scaleDesignVec(iterationData*);
	void multShoot_createContCons(iterationData*, tpat_nodeset*);
	void multShoot_getSimICs(iterationData*, tpat_nodeset*, int, double*, double*, double*);
	double multShoot_getSlackVarVal(iterationData*, tpat_constraint);
	void multShoot_applyConstraint(iterationData*, tpat_constraint, int);
	tpat_nodeset* multShoot_createOutput(iterationData*, tpat_nodeset*, bool);

protected:
	void multShoot_targetPosVelCons(iterationData*, tpat_constraint, int);
	void multShoot_targetExContCons(iterationData*, tpat_constraint, int);
	void multShoot_targetState(iterationData*, tpat_constraint, int);
	void multShoot_targetDeltaV(iterationData*, tpat_constraint, int);
	void multShoot_targetDist(iterationData*, tpat_constraint, int);
	double multShoot_targetDist_compSlackVar(iterationData*, tpat_constraint);
	void multShoot_targetSP(iterationData*, tpat_constraint, int);
	void multShoot_targetSP_mag(iterationData*, tpat_constraint, int);
	void multShoot_targetSP_dist(iterationData*, tpat_constraint, int);
	double multShoot_targetSPMag_compSlackVar(iterationData*, tpat_constraint);
	double multShoot_targetSP_maxDist_compSlackVar(iterationData*, tpat_constraint);
};

#endif