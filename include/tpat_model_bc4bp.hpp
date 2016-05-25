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
class TPAT_Sys_Data_BC4BP;

/**
 *	@brief A derivative of the TPAT_Model class, specific to the BCR4BPR
 *
 *	This dynamic model overrides many of the base class's functions to add
 *	support for epoch-dependencies present in this non-autonomous system.
 *	It also adds support for SP targeting, which is only currently supported 
 *	in the BCR4BP.
 */
class TPAT_Model_BC4BP : public TPAT_Model{
public:
	TPAT_Model_BC4BP();
	TPAT_Model_BC4BP(const TPAT_Model_BC4BP&);
	~TPAT_Model_BC4BP() {}

	TPAT_Model_BC4BP& operator=(const TPAT_Model_BC4BP&);

	// Core Functions
	TPAT_Model::eom_fcn getFullEOM_fcn() const;
	TPAT_Model::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const TPAT_Sys_Data*) const;
	std::vector<double> getPrimVel(double, const TPAT_Sys_Data*) const;
	void sim_saveIntegratedData(const double*, double, TPAT_Traj*) const;
	bool sim_locateEvent(TPAT_Event, TPAT_Traj*, const double*, double, double, TPAT_Verbosity_Tp) const;

	// Static Calculation Functions
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getPrimaryPos(double, const TPAT_Sys_Data_BC4BP*, double*);
	static void getPrimaryVel(double, const TPAT_Sys_Data_BC4BP*, double*);
	static void getPrimaryAccel(double, const TPAT_Sys_Data_BC4BP*, double*);
	static void orientAtEpoch(double, TPAT_Sys_Data_BC4BP*);
	
	// Multiple Shooting functions
	void multShoot_initDesignVec(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	void multShoot_initIterData(TPAT_MultShoot_Data *it) const;
	void multShoot_scaleDesignVec(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	void multShoot_createContCons(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	void multShoot_getSimICs(const TPAT_MultShoot_Data*, const TPAT_Nodeset*, int, double*, double*, double*) const;
	double multShoot_getSlackVarVal(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
	void multShoot_applyConstraint(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_createOutput(const TPAT_MultShoot_Data*, const TPAT_Nodeset*, bool, TPAT_Nodeset*) const;

protected:
	void multShoot_targetCont_PosVel(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetCont_Ex(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetCont_Ex_Seg(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetState(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetDeltaV(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetDist(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	double multShoot_targetDist_compSlackVar(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
	void multShoot_targetApse(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetSP(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetSP_mag(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetSP_dist(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	double multShoot_targetSPMag_compSlackVar(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
	double multShoot_targetSP_maxDist_compSlackVar(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
};

#endif