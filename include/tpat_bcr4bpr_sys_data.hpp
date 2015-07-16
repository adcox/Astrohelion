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
#ifndef H_BCR4BPR_SYS_DATA
#define H_BCR4BPR_SYS_DATA

#include "tpat_sys_data.hpp"
#include "matio.h"
#include "tpat_constants.hpp"

/**
 *	@brief This derivative class of the tpat_sys_data super-class
 *	contains information specific to the BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 18, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_bcr4bpr_sys_data : public tpat_sys_data{
	public:
		tpat_bcr4bpr_sys_data();
		tpat_bcr4bpr_sys_data(std::string, std::string, std::string);
		tpat_bcr4bpr_sys_data(const tpat_bcr4bpr_sys_data&);

		tpat_bcr4bpr_sys_data& operator=(const tpat_bcr4bpr_sys_data&);
		
		double getMu() const;
		double getNu() const;
		double getK() const;
		double getCharLRatio() const;

		double getTheta0() const;
		double getPhi0() const;
		double getGamma() const;

		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);

		void saveToMat(mat_t*);
		void readFromMat(mat_t*);

		/** Time when geometry is at reference orientation (theta = phi = 0), seconds, J2000, UTC */
		static double REF_EPOCH;	// 2005/06/21 18:21:35

	private:
		void initFromPrimNames(std::string, std::string, std::string);
};

#endif
//END