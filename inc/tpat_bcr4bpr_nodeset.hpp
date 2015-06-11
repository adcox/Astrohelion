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

#ifndef __H_BCR4BPR_NODESET__
#define __H_BCR4BPR_NODESET__

#include "tpat_nodeset.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "matio.h"
 
#include <vector>

/**
 *	@brief This derivative of the tpat_nodeset object contains additional information
 *	for the BCR4BP
 *
 *	@author Andrew Cox
 *	@version May 21, 2015
 * 	@copyright GNU GPL v3.0
 */
class tpat_bcr4bpr_nodeset : public tpat_nodeset{
	public:
		tpat_bcr4bpr_nodeset() : tpat_nodeset(6){}	//!< Default, do-nothing constructor
		tpat_bcr4bpr_nodeset(tpat_bcr4bpr_sys_data);
		tpat_bcr4bpr_nodeset(double[6], tpat_bcr4bpr_sys_data, double, double, int);
		tpat_bcr4bpr_nodeset(double[6], tpat_bcr4bpr_sys_data, double, double, int,
			node_distro_t);
		tpat_bcr4bpr_nodeset(const tpat_bcr4bpr_nodeset&);
		~tpat_bcr4bpr_nodeset();
		

		tpat_bcr4bpr_nodeset& operator =(const tpat_bcr4bpr_nodeset&);
		friend tpat_bcr4bpr_nodeset operator +(const tpat_bcr4bpr_nodeset &, const tpat_bcr4bpr_nodeset &);
		
		std::vector<double>* getEpochs();
		
		double getEpoch(int) const;
		tpat_sys_data* getSysData();

		void appendEpoch(double);
		void print() const;
		void saveToMat(const char*);
	private:
		std::vector<double> epochs;			//!< Vector of epochs, one for each node
		tpat_bcr4bpr_sys_data sysData;		//!< System data object

		void initEpochs(int, double);
		void saveEpochs(mat_t*);
};

#endif