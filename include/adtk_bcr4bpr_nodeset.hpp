/**
 *	This derivative of the adtk_nodeset object contains additional information
 *	for the BCR4BP
 *
 *	Author: Andrew Cox
 *	Version: May 21, 2015
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

#ifndef __H_BCR4BPR_NODESET__
#define __H_BCR4BPR_NODESET__

#include "adtk_nodeset.hpp"
#include "adtk_bcr4bpr_constraint.hpp"

#include <vector>

// Forward declarations
class adtk_bcr4bpr_sys_data;

class adtk_bcr4bpr_nodeset : public adtk_nodeset{
	public:
		adtk_bcr4bpr_nodeset() : adtk_nodeset(6){}
		adtk_bcr4bpr_nodeset(double[6], adtk_bcr4bpr_sys_data, double, double, int,
			node_distro_t);

		adtk_bcr4bpr_constraint getConstraint(int);
		void addConstraint(adtk_bcr4bpr_constraint);

	private:
		std::vector<adtk_bcr4bpr_constraint> constraints;
		std::vector<double> epochs;
};

#endif