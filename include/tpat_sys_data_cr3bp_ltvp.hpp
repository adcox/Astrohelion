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

#ifndef H_CR3BP_LTVP_SYS_DATA
#define H_CR3BP_LTVP_SYS_DATA

#include "tpat_sys_data_cr3bp.hpp"

#include "tpat_model_cr3bp_ltvp.hpp"
 
/**
 *	@brief System data object for the CR3BP Low Thrust, Velocity-Pointing system
 */
class tpat_sys_data_cr3bp_ltvp : public tpat_sys_data_cr3bp{

	public:
		tpat_sys_data_cr3bp_ltvp();
		tpat_sys_data_cr3bp_ltvp(std::string, std::string, double, double, double);
		tpat_sys_data_cr3bp_ltvp(const tpat_sys_data_cr3bp_ltvp&);
		tpat_sys_data_cr3bp_ltvp(const char*);
		
		tpat_sys_data_cr3bp_ltvp& operator=(const tpat_sys_data_cr3bp_ltvp&);
		
		const tpat_model* getModel() const;
		
		double getIsp() const;
		double getThrust() const;
		double getM0() const;

		void setIsp(double);
		void setIspDim(double);
		void setM0(double);
		void setM0Dim(double);
		void setThrust(double);
		void setThrustDim(double);

		void saveToMat(const char*) const;
		void saveToMat(mat_t*) const;
	private:
		/** The dynamic model that governs motion for this system*/
		tpat_model_cr3bp_ltvp model = tpat_model_cr3bp_ltvp();
		void readFromMat(mat_t*);
};

#endif