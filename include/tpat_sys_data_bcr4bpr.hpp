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
#ifndef H_BCR4BPR_SYS_DATA
#define H_BCR4BPR_SYS_DATA

#include "tpat_sys_data.hpp"

#include "tpat_constants.hpp" 
#include "tpat_model_bcr4bpr.hpp"

#include "matio.h"

/**
 *	@brief This derivative class of the tpat_sys_data super-class
 *	contains information specific to the BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 18, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_sys_data_bcr4bpr : public tpat_sys_data{
	public:
		tpat_sys_data_bcr4bpr();
		tpat_sys_data_bcr4bpr(std::string, std::string, std::string);
		tpat_sys_data_bcr4bpr(const tpat_sys_data_bcr4bpr&);
		tpat_sys_data_bcr4bpr(const char*);
		
		tpat_sys_data_bcr4bpr& operator=(const tpat_sys_data_bcr4bpr&);
		
		const tpat_model* getModel() const;

		double getMu() const;
		double getNu() const;
		double getK() const;
		double getCharLRatio() const;

		double getEpoch0() const;
		double getTheta0() const;
		double getPhi0() const;
		double getGamma() const;

		void setEpoch0(double T);
		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);

		void saveToMat(const char*) const;
		void saveToMat(mat_t*) const;

		/** Time when geometry is at reference orientation (theta = phi = 0), seconds, J2000, UTC */
		static double REF_EPOCH;	// 2005/06/21 18:21:35

	private:
		/** The dynamic model that governs motion for this system*/
		tpat_model_bcr4bpr model = tpat_model_bcr4bpr();
		
		void initFromPrimNames(std::string, std::string, std::string);
		void readFromMat(mat_t*);
};

#endif
//END