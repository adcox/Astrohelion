/**
 *  @file Nodeset_cr3bp.hpp
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

#include "Nodeset.hpp"


namespace astrohelion{

// Forward Declarations
class SysData_cr3bp;
class Traj_cr3bp;

/**
 *	@ingroup traj cr3bp
 *	@brief This derivative of the Nodeset contains additional information for 
 *	the BCR4BPR
 *
 *	Nodes are 6-dimensional, with three position states and three velocity states. Times-
 *	of-flight between nodes are recorded in a separate vector.
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class Nodeset_cr3bp : public Nodeset{

public:
	/**
	 *  @name *structors
	 *  @{
	 */
	Nodeset_cr3bp(const SysData_cr3bp*);
	Nodeset_cr3bp(const SysData_cr3bp*, const double[6], double, int, NodeDistro_tp type = NodeDistro_tp::TIME);
	Nodeset_cr3bp(const SysData_cr3bp*, std::vector<double>, double, int, NodeDistro_tp type = NodeDistro_tp::TIME);
	Nodeset_cr3bp(Traj_cr3bp, int);
	Nodeset_cr3bp(Traj_cr3bp, int, NodeDistro_tp);
	Nodeset_cr3bp(const Nodeset_cr3bp&, int, int);
	Nodeset_cr3bp(const Nodeset_cr3bp&);
	Nodeset_cr3bp(const BaseArcset&);
	baseArcsetPtr create(const SysData*) const override;
	baseArcsetPtr clone() const override;
	//@}

	// Operators

	/**
	 *  @name Set and Get Functions
	 *  @{
	 */
	double getJacobi(int) const;
	double getJacobiByIx(int) const;
	void setJacobi(int, double);
	void setJacobiByIx(int, double);
	//@}
protected:
	// Utility
	void initExtraParam() override;
};

}