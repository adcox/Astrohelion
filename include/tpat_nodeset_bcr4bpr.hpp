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

#ifndef H_BCR4BPR_NODESET
#define H_BCR4BPR_NODESET

#include "tpat_nodeset.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "matio.h"
 
#include <vector>

/**
 *	@brief This derivative of the tpat_nodeset object contains additional information
 *	for the BCR4BP
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes and epoch times at each node are stored in separate vectors.
 *
 *	@author Andrew Cox
 *	@version May 21, 2015
 * 	@copyright GNU GPL v3.0
 */
class tpat_nodeset_bcr4bpr : public tpat_nodeset{
	public:
		tpat_nodeset_bcr4bpr() : tpat_nodeset(6){}	//!< Default, do-nothing constructor
		tpat_nodeset_bcr4bpr(tpat_sys_data_bcr4bpr);
		tpat_nodeset_bcr4bpr(double[6], tpat_sys_data_bcr4bpr, double, double, int);
		tpat_nodeset_bcr4bpr(std::vector<double>, tpat_sys_data_bcr4bpr, double, double, int);
		tpat_nodeset_bcr4bpr(double[6], tpat_sys_data_bcr4bpr, double, double, int,
			node_distro_t);
		tpat_nodeset_bcr4bpr(std::vector<double>, tpat_sys_data_bcr4bpr, double, double, int,
			node_distro_t);
		tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr&, int, int);
		tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr&);
		~tpat_nodeset_bcr4bpr();
		

		tpat_nodeset_bcr4bpr& operator =(const tpat_nodeset_bcr4bpr&);
		friend tpat_nodeset_bcr4bpr operator +(const tpat_nodeset_bcr4bpr &, const tpat_nodeset_bcr4bpr &);
		
		std::vector<double>* getEpochs();
		
		double getEpoch(int) const;
		tpat_sys_data* getSysData();

		void appendEpoch(double);
		void insertEpoch(int, double);
		void reverseOrder();
		void print() const;
		void saveToMat(const char*);

	private:
		std::vector<double> epochs;			//!< Vector of epochs, one for each node
		tpat_sys_data_bcr4bpr sysData;		//!< System data object

		void initEpochs(int, double);
		void saveEpochs(mat_t*);
};

#endif