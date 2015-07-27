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
 *	@brief Derivative of tpat_model, specific to the CR3BP-LTVP
 *
 *	Under construction. Simulation is fully supported in this model,
 *	but the corrections process will fall back to the base model 
 *	behavior, which may produce unexpected results.
 */
class tpat_model_cr3bp_ltvp : public tpat_model{
public:
	tpat_model_cr3bp_ltvp();
	tpat_model_cr3bp_ltvp(const tpat_model_cr3bp_ltvp&);
	~tpat_model_cr3bp_ltvp() {}
	
	tpat_model_cr3bp_ltvp& operator=(const tpat_model_cr3bp_ltvp&);

	// Core Functions
	tpat_model::eom_fcn getFullEOM_fcn();
	tpat_model::eom_fcn getSimpleEOM_fcn();
	std::vector<double> getPrimPos(double, tpat_sys_data*);
	std::vector<double> getPrimVel(double, tpat_sys_data*);
	void sim_saveIntegratedData(double*, double, tpat_traj*);
	bool sim_locateEvent(tpat_event, tpat_traj*, tpat_model*, double*, double, double, bool);

	// Corrector Functions
	tpat_nodeset* corrector_createOutput(iterationData*, tpat_nodeset*, bool);
};

#endif