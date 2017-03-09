/**
 *  @file Traj.hpp
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

#include "BaseArcset.hpp"
#include "matio.h"

namespace astrohelion{

// Forward Declarations
class Nodeset;
class Traj;

/**
 * @brief Smart pointer to a Traj object
 */
typedef std::shared_ptr<Traj> trajPtr;

/**
 *	@ingroup traj
 *	@brief Contains information about a series of continuous states along a trajectory.
 *	
 *	The trajecotory object holds information about a trajectory in one package
 *	to make passing the data between engines and analysis tools simple. This class acts
 *	as a base class for model-specific trajectory objects. This general object stores the
 *	following information about a trajectory:
 *	* A vector of trajectory steps, each of which contains:
 *		* State  (x, y, z, vx, vy, vz [non-dimensional])
 *		* Accleration (ax, ay, az [non-dimensional])
 *		* State Transition Matrix (6x6 [non-dimensional])
 *		* Time since first step [non-dimensional]
 * 	* A system data object that describes the system the trajectory was integrated in
 * 	* A tolerance value that describes the minimum numerical accuracy used to create the trajectory
 * 	
 *	@author Andrew Cox
 *	@version August 29, 2015
 *	@copyright GNU GPL v3.0
 *	
 *	@see BaseArcset
 */
class Traj : public BaseArcset{

public:
	/**
	 *  @name *structors
	 *  @{
	 */
	Traj(const SysData*);
	Traj(const Traj&);
	Traj(const BaseArcset&);
	virtual ~Traj();
	virtual baseArcsetPtr create(const SysData*) const;
	virtual baseArcsetPtr clone() const;
	static Traj fromNodeset(Nodeset);
	//@}

	// Operators
	friend Traj operator +(const Traj&, const Traj&);
	virtual Traj& operator +=(const Traj&);

	/**
	 *  @name Set and Get Functions
	 *  @{
	 */
	virtual double getTotalTOF() const override;
	double getTimeByIx(int) const;
	void setTimeByIx(int, double);
	void shiftAllTimes(double);
	//@}

	// Utility Functions
	Nodeset discretize(int) const;
	void print() const;
	virtual void readFromMat(const char*);
	virtual void saveToMat(const char*) const;
};


}// END of Astrohelion namespace