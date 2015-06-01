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

#ifndef __H_CR3BP_CONSTRAINT__
#define __H_CR3BP_CONSTRAINT__

#include "adtk_constraint.hpp"

/**
 *	@brief This derivative of the adtk_constraint class applies specifically to 
 *	CR3BP corrections processes
 *
 *	specifies the size of the node to be 6
 *
 *	@author Andrew Cox
 *	@version May 21, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_cr3bp_constraint : public adtk_constraint{
	public:
		/**
		 *	Create a generic CR3BP constraint
		 */
		adtk_cr3bp_constraint() : adtk_constraint(6){}

		/**
		 *	Create a CR3BP constraint with specified type
		 *	@param t constraint type
		 */
		adtk_cr3bp_constraint(constraint_t t) : adtk_constraint(6, t){}

		/**
		 *	Create a CR3BP constraint with specified type, node, and data
		 *	@param t constraint type
		 *	@param i node index
		 *	@param d data (6 elements)
		 */
		adtk_cr3bp_constraint(constraint_t t, int i, std::vector<double> d) : adtk_constraint(6, t, i, d){}
		
		/**
		 *	Create a CR3BP constraint with specified type, node, and data
		 *	@param t constraint type
		 *	@param i node index
		 *	@param d data (6 elements)
		 */
		adtk_cr3bp_constraint(constraint_t t, int i, double d[6]) : adtk_constraint(6, t, i, d){}
};

#endif