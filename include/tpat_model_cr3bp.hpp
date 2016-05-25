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
#ifndef H_MODEL_CR3BP
#define H_MODEL_CR3BP

#include "tpat_model.hpp"

// Forward declarations
class TPAT_Sys_Data_CR3BP;

/**
 *	@brief Derivative of TPAT_Model, specific to the CR3BP
 *
 *	The base class's methods provide a good framework for this system,
 *	so only minimal adjustments are needed. This model adds support for 
 *	Jacobi-Constant targeting.
 */
class TPAT_Model_CR3BP : public TPAT_Model{
public:
	TPAT_Model_CR3BP();
	TPAT_Model_CR3BP(const TPAT_Model_CR3BP&);
	~TPAT_Model_CR3BP() {}

	TPAT_Model_CR3BP& operator=(const TPAT_Model_CR3BP&);

	// Core Functions
	TPAT_Model::eom_fcn getFullEOM_fcn() const;
	TPAT_Model::eom_fcn getSimpleEOM_fcn() const;
	std::vector<double> getPrimPos(double, const TPAT_Sys_Data*) const;
	std::vector<double> getPrimVel(double, const TPAT_Sys_Data*) const;
	
	// Static Calculation Functions
	static int fullEOMs(double, const double[], double[], void*);
	static int simpleEOMs(double, const double[], double[], void*);
	static void getEquilibPt(const TPAT_Sys_Data_CR3BP*, int, double, double[3]);
	static double getJacobi(const double[], double);
	static void getUDDots(double, double, double, double, double* ddots);

	// Simulation Engine Functions
	void sim_saveIntegratedData(const double*, double, TPAT_Traj*) const;
	bool sim_locateEvent(TPAT_Event, TPAT_Traj*, const double*, double, double, TPAT_Verbosity_Tp) const;

	// Multiple Shooting Functions
	void multShoot_applyConstraint(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_createOutput(const TPAT_MultShoot_Data*, const TPAT_Nodeset*, bool, TPAT_Nodeset*) const;
	void multShoot_initIterData(TPAT_MultShoot_Data *it) const;
	
protected:
	void multShoot_targetJC(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	void multShoot_targetPseudoArc(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
};

#endif