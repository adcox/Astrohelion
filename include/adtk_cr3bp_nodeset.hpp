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
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __H_CR3BP_NODESET__
#define __H_CR3BP_NODESET__

#include "adtk_nodeset.hpp"
#include "adtk_cr3bp_sys_data.hpp"

#include <vector>

// Forward declarations
class adtk_cr3bp_sys_data;

/**
 *	@brief This derivative of the adtk_nodeset contains additional information for 
 *	the BCR4BPR
 *
 *	enforces a node size of 6
 *
 *	@author Andrew Cox
 *	@version May 21, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_cr3bp_nodeset : public adtk_nodeset{
	public:
		// *structors
		adtk_cr3bp_nodeset();
		adtk_cr3bp_nodeset(adtk_cr3bp_sys_data data);
		adtk_cr3bp_nodeset(double[6], adtk_cr3bp_sys_data, double, int, node_distro_t);

		adtk_cr3bp_nodeset(const adtk_cr3bp_nodeset&);

		// Operators
		adtk_cr3bp_nodeset& operator =(const adtk_cr3bp_nodeset&);

		// Set and Gets
		adtk_sys_data* getSysData();
		void print() const;
	private:
		adtk_cr3bp_sys_data sysData; 	//!< Object holding information about the dynamical system
};

#endif