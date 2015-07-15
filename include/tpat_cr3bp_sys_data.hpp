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
#ifndef __H_CR3BP_SYS_DATA_
#define __H_CR3BP_SYS_DATA_

#include "tpat_sys_data.hpp"
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
class tpat_cr3bp_sys_data : public tpat_sys_data{
	public:
		tpat_cr3bp_sys_data();
		tpat_cr3bp_sys_data(std::string, std::string);
		tpat_cr3bp_sys_data(const tpat_cr3bp_sys_data&);

		tpat_cr3bp_sys_data& operator=(const tpat_cr3bp_sys_data&);
		
		double getMu() const;

		void saveToMat(mat_t*);
		void readFromMat(mat_t*);
		
	private:
		void initFromPrimNames(std::string, std::string);
};

#endif
//END