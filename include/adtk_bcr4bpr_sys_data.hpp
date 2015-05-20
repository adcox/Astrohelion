/**
 *	A derivative class of the adtk_sys_data super-class. This class
 *	contains information specific to the BCR4BP, rotating coordinates
 *
 *	Author: Andrew Cox
 *
 *	Version: May 18, 2015
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
#ifndef __H_BCR4BPR_SYS_DATA_
#define __H_BCR4BPR_SYS_DATA_

#include "adtk_sys_data.hpp"

class adtk_bcr4bpr_sys_data : public adtk_sys_data{
	public:
		adtk_bcr4bpr_sys_data();
		adtk_bcr4bpr_sys_data(std::string, std::string, std::string);

		adtk_bcr4bpr_sys_data& operator=(const adtk_bcr4bpr_sys_data&);
		
		double getMu();
		double getNu();
		double getK();
		double getCharLRatio();
		std::string getPrimary(int n);	//We override this function, so re-declare it
	private:
		/** Mass ratio between P2 + P3 and total*/
		double mu;

		/** Mass ratio between P3 and total */
		double nu;

		/** Scaling constant, non-dim */
		double k;

		/** Ratio between P3's orbital radius and P2's orbital radius */
		double charLRatio;

		/** Name of P1 */
		std::string P1;

		/** Name of P2 */
		std::string P2;

		/** Name of P3 */
		std::string P3;
};

#endif
//END