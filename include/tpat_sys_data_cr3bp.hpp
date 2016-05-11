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
#ifndef H_CR3BP_SYS_DATA
#define H_CR3BP_SYS_DATA

#include "tpat_sys_data.hpp"

#include "tpat_model_cr3bp.hpp"

#include "matio.h"
#include <string>

/**
 *	@brief A derivative class of the tpat_sys_data object which
 *	contains information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_sys_data_cr3bp : public tpat_sys_data{
	public:
		tpat_sys_data_cr3bp();
		tpat_sys_data_cr3bp(std::string, std::string);
		tpat_sys_data_cr3bp(const tpat_sys_data_cr3bp&);
		tpat_sys_data_cr3bp(const char*);
		
		tpat_sys_data_cr3bp& operator=(const tpat_sys_data_cr3bp&);
		
		const tpat_model* getModel() const;
		double getMu() const;

		void saveToMat(const char*) const;
		void saveToMat(mat_t*) const;
		
	protected:
		void initFromPrimNames(std::string, std::string);
		void readFromMat(mat_t*);
		
	private:
		/** The dynamic model that governs motion for this system*/
		tpat_model_cr3bp model = tpat_model_cr3bp();
};

#endif
//END