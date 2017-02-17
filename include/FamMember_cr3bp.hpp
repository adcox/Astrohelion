/**
 *  @file FamMember_cr3bp.hpp
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

#include "Core.hpp"
 
#include "Common.hpp"
#include "EigenDefs.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace astrohelion{

class Traj_cr3bp;
class SysData_cr3bp;

/**
 *	@ingroup traj cr3bp
 *	@brief A data object to store information about a family member
 *	
 *	The purpose of this object is to store just enough information to re-create the
 *	orbit using the simulation engine with additional variables to describe the orbit
 *	for categorization purposes (e.g. the getMemberBy___ functions in Fam_cr3bp)
 */
class FamMember_cr3bp : public Core{
	public:
		/**
		 *  @name *structors
		 *  @{
		 */
		FamMember_cr3bp(){}
		FamMember_cr3bp(double*, double, double, double, double, double);
		FamMember_cr3bp(const Traj_cr3bp);
		FamMember_cr3bp(const FamMember_cr3bp&);
		~FamMember_cr3bp();
		//@}

		FamMember_cr3bp& operator= (const FamMember_cr3bp&);

		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		std::vector<cdouble> getEigVals() const;
		MatrixXRcd getEigVecs() const;
		std::vector<double> getIC() const;
		double getTOF() const;
		MatrixXRd getSTM() const;
		double getJacobi() const; 
		double getXAmplitude() const;
		double getYAmplitude() const;
		double getZAmplitude() const;

		void setEigVals(std::vector<cdouble>);
		void setEigVecs(MatrixXRcd);
		void setIC(std::vector<double>);
		void setJacobi(double);
		void setSTM(MatrixXRd);
		void setTOF(double);
		void setXAmplitude(double);
		void setYAmplitude(double);
		void setZAmplitude(double);
		//@}

		Traj_cr3bp toTraj(const SysData_cr3bp*);

	protected:
		/** Vector of 6 eigenvalues; initialized to NAN by default */
		MatrixXRd stm = MatrixXRd::Identity(6,6);
		MatrixXRcd eigVecs = MatrixXRcd::Zero(6,6);
		std::vector<cdouble> eigVals {{NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}, {NAN,0}};
		std::vector<double> IC {0,0,0,0,0,0};		//!< Initial state for this trajectory, non-dim units
		double TOF = NAN;							//!< Time of flight for traj., non-dim units
		double JC = NAN;							//!< Jacobi constant for traj., non-dimensional
		double xAmplitude = NAN;					//!< Max width in x-direction, non-dim units
		double yAmplitude = NAN;					//!< Max width in y-direction, non-dim units
		double zAmplitude = NAN;					//!< Max width in z-direction, non-dim units

		void copyMe(const FamMember_cr3bp&);
};

}