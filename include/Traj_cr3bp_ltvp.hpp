/**
 *  @file Traj_cr3bp_ltvp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
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

#include "Traj.hpp"


namespace astrohelion{

// Forward Declarations
class SysData_cr3bp_ltvp;

/**
 *	@ingroup traj
 *	@brief A derivative class of the Traj object, which
 *	contains trajectory information specific to the Low Thrust, 
 *	Velocity-Pointing CR3BP
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class Traj_cr3bp_ltvp : public Traj{

public:
	/**
	 *  @name *structors
	 *  @{
	 */
	Traj_cr3bp_ltvp(const SysData_cr3bp_ltvp*);
	Traj_cr3bp_ltvp(const Traj_cr3bp_ltvp&);
	Traj_cr3bp_ltvp(const BaseArcset&);
	baseArcsetPtr create(const SysData*) const override;
	baseArcsetPtr clone() const;
	//@}

	/**
	 *  @name Set and Get Functions
	 *  @{
	 */
	double getJacobiByIx(int) const;
	double getMassByIx(int) const;

	void setJacobiByIx(int, double);
	void setMassByIx(int, double);
	//@}
	
	void saveToMat(const char*) const override;
private:
	void initExtraParam();
};

}// END of Astrohelion namespace