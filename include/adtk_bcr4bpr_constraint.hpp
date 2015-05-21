/**
 *	This derivative constraint class specifies the size of the node to be 7
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

#ifndef __H_BCR4BPR_CONSTRAINT__
#define __H_BCR4BPR_CONSTRAINT__

#include "adtk_constraint.hpp"

class adtk_bcr4bpr_constraint : public adtk_constraint{
	public:
		adtk_bcr4bpr_constraint() : adtk_constraint(7){}
		adtk_bcr4bpr_constraint(constraint_t t) : adtk_constraint(7, t){}
		adtk_bcr4bpr_constraint(constraint_t t, int i, std::vector<double> d) : adtk_constraint(7, t, i, d){}
		adtk_bcr4bpr_constraint(constraint_t t, int i, double d[7]) : adtk_constraint(7, t, i, d){}
};

#endif