/**
 *	The nodeset object is similar to a trajectory object, but a nodeset only contains a few
 *	distinct states, or "nodes" and is used in corrections processes to break a trajectory
 *	into smaller pieces, which can improve the corrector's performance.
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

#ifndef __H_NOTESET__
#define __H_NOTESET__

#include "adtk_constraint.hpp"

#include <vector>

class adtk_nodeset{
	public:
		adtk_nodeset();
		adtk_nodeset(const adtk_nodeset&);

		adtk_nodeset& operator =(const adtk_nodeset&);

		std::vector<double> getNode(int);
		double getTOF(int);
		double getNumNodes();
		int getNodeSize();

		void appendNode(std::vector<double>);
		void appendTOF(double);

	protected:
		int nodeSize = 6;
		std::vector<double> nodes;
		std::vector<double> tofs;

		adtk_constraint velCont;
		std::vector<adtk_constraint> constraints;

};

#endif