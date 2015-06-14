/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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
#ifndef __H_BCR4BPR_TRAJ_
#define __H_BCR4BPR_TRAJ_

#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_trajectory.hpp"

#include "matio.h"
 
/**
 *	@brief A derivative class of the tpat_trajectory object that
 *	contains trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_bcr4bpr_traj : public tpat_trajectory{

	public:
		tpat_bcr4bpr_traj();
		tpat_bcr4bpr_traj(int);
		tpat_bcr4bpr_traj(tpat_bcr4bpr_sys_data);
		tpat_bcr4bpr_traj(const tpat_bcr4bpr_traj&);

		// Operators
		tpat_bcr4bpr_traj& operator= (const tpat_bcr4bpr_traj&);
		friend tpat_bcr4bpr_traj operator +(const tpat_bcr4bpr_traj&, const tpat_bcr4bpr_traj&);

		// Set and Get Functions
		double getTheta0();
		double getPhi0();
		double getGamma();
		tpat_bcr4bpr_sys_data getSysData();
		tpat_sys_data::system_t getType() const;

		std::vector<double>* get_dqdT();
		std::vector<double> get_dqdT(int);
		
		void setLength();
		void setSysData(tpat_bcr4bpr_sys_data);
		void saveToMat(const char*);
	private:
		/** Derivatives of the state variables (pos, vel) with respect to Epoch Time; 
			used in corrections processes */
		std::vector<double> dqdT;

		/** A system data object specific to the BCR4BPR */
		tpat_bcr4bpr_sys_data sysData;

		void save_dqdT(mat_t*);
};

#endif
//END