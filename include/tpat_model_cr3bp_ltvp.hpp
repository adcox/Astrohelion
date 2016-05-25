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
#ifndef H_MODEL_CR3BP_LTVP
#define H_MODEL_CR3BP_LTVP

#include "tpat_model.hpp"

/**
 *	@brief Derivative of TPAT_Model, specific to the CR3BP-LTVP
 *
 *	Under construction. Simulation is fully supported in this model,
 *	but the corrections process will fall back to the base model 
 *	behavior, which may produce unexpected results.
 */
class TPAT_Model_CR3BP_LTVP : public TPAT_Model{
public:
	TPAT_Model_CR3BP_LTVP();
	TPAT_Model_CR3BP_LTVP(const TPAT_Model_CR3BP_LTVP&);
	~TPAT_Model_CR3BP_LTVP() {}
	
	TPAT_Model_CR3BP_LTVP& operator=(const TPAT_Model_CR3BP_LTVP&);

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

	// Multiple Shooting Functions
	void multShoot_createOutput(const TPAT_MultShoot_Data*, const TPAT_Nodeset*, bool, TPAT_Nodeset*) const;
	void multShoot_initIterData(TPAT_MultShoot_Data *it) const;
};

#endif