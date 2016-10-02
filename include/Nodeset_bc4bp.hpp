/**
 *  @file Nodeset_bc4bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Nodeset.hpp"

#include "matio.h"


namespace astrohelion{

// Forward Declarations
class Event;
class SysData_bc4bp;

/**
 *	@ingroup traj bc4bp
 *	@brief his derivative of the Nodeset object contains additional information
 *	for the BCR4BP
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes and epoch times at each node are also stored
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class Nodeset_bc4bp : public Nodeset{

public:
	/**
	 *  @name *structors
	 *  @{
	 */
	Nodeset_bc4bp(const SysData_bc4bp*);
	Nodeset_bc4bp(const double[6], const SysData_bc4bp*, double, double, int);
	Nodeset_bc4bp(std::vector<double>, const SysData_bc4bp*, double, double, int);
	Nodeset_bc4bp(const double[6], const SysData_bc4bp*, double, double, int,
		NodeDistro_tp);
	Nodeset_bc4bp(std::vector<double>, const SysData_bc4bp*, double, double, int,
		NodeDistro_tp);
	// Nodeset_bc4bp(const Nodeset_bc4bp&, int, int);
	Nodeset_bc4bp(const Nodeset_bc4bp&);
	Nodeset_bc4bp(const BaseArcset&);
	baseArcsetPtr create(const SysData*) const override;
	baseArcsetPtr clone() const override;
	//@}
	
	// Operators

	// Set and Get Functions

	// Utility Functions
	
private:
	void initExtraParam() override;
};

}