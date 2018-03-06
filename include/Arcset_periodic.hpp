/**
 * @file Arcset_periodic.hpp
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrohelion.
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "Arcset.hpp"

#include "EigenDefs.hpp"
#include "SysData.hpp"

namespace astrohelion{

/**
 *  \ingroup traj
 *  @brief [brief description]
 *  @details [long description]
 * 
 *  @param  [description]
 *  @return [description]
 */
class Arcset_periodic : public Arcset{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	Arcset_periodic(const SysData*);
	Arcset_periodic(const Arcset_periodic&);
	Arcset_periodic(const Arcset&);
	~Arcset_periodic();
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */

	void getEigData(std::vector<cdouble>*, MatrixXRcd*);
	MatrixXRd getMonodromy();

	double getXAmp() const;
	double getYAmp() const;
	double getZAmp() const;
	double getAmp(unsigned int) const;
	//\}
};

}// END of astrohelion namespace