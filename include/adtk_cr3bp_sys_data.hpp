/**
 *	A derivative class of the adtk_sys_data super-class. This class
 *	contains information specific to the CR3BP
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */

/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_CR3BP_SYS_DATA_
#define __H_CR3BP_SYS_DATA_

#include "adtk_sys_data.hpp"

class adtk_cr3bp_sys_data : public adtk_sys_data{
	private:
		/** Mass ratio between two primaries*/
		double mu;
		std::string P1;
		std::string P2;

	public:
		adtk_cr3bp_sys_data();
		adtk_cr3bp_sys_data(std::string P1, std::string P2);

		adtk_cr3bp_sys_data& operator=(const adtk_cr3bp_sys_data&);
		
		double getMu() const;
		std::string getPrimary(int n) const;	//We override this function, so re-declare it
};

#endif
//END