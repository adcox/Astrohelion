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

#ifndef H_FAMILY_MEMBER
#define H_FAMILY_MEMBER

#include "tpat_constants.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

class tpat_traj_cr3bp;
class tpat_sys_data_cr3bp;

/**
 *	@brief A data object to store information about a family member
 *	
 *	The purpose of this object is to store just enough information to re-create the
 *	orbit using the simulation engine with additional variables to describe the orbit
 *	for categorization purposes (e.g. the getMemberBy___ functions in tpat_family_cr3bp)
 */
class tpat_family_member_cr3bp{
	public:
		tpat_family_member_cr3bp(){}
		tpat_family_member_cr3bp(double*, double, double, double, double, double);
		tpat_family_member_cr3bp(const tpat_traj_cr3bp);
		tpat_family_member_cr3bp(const tpat_family_member_cr3bp&);
		~tpat_family_member_cr3bp();

		tpat_family_member_cr3bp& operator= (const tpat_family_member_cr3bp&);

		std::vector<cdouble> getEigVals() const;
		std::vector<double> getIC() const;
		double getTOF() const;
		double getJacobi() const; 
		double getXAmplitude() const;
		double getYAmplitude() const;
		double getZAmplitude() const;

		void setEigVals(std::vector<cdouble>);
		void setIC(std::vector<double>);
		void setJacobi(double);
		void setTOF(double);
		void setXAmplitude(double);
		void setYAmplitude(double);
		void setZAmplitude(double);
		tpat_traj_cr3bp toTraj(const tpat_sys_data_cr3bp*);

	protected:
		/** Vector of 6 eigenvalues; initialized to NAN by default */
		std::vector<cdouble> eigVals {{NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}};
		std::vector<double> IC {0,0,0,0,0,0};	//!< Initial state for this trajectory, non-dim units
		double TOF = NAN;			//!< Time of flight for traj., non-dim units
		double JC = NAN;			//!< Jacobi constant for traj., non-dimensional
		double xAmplitude = NAN;		//!< Max width in x-direction, non-dim units
		double yAmplitude = NAN;		//!< Max width in y-direction, non-dim units
		double zAmplitude = NAN;		//!< Max width in z-direction, non-dim units

		void copyMe(const tpat_family_member_cr3bp&);
};

#endif