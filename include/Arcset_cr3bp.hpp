/**
 *  @file Arcset_cr3bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 1, 2017
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Arcset.hpp"

namespace astrohelion{

// Forward Declarations
class SysData_cr3bp;

/**
 *	\ingroup traj cr3bp
 *	@brief Supplies functions specific to CR3BP arcsets
 *
 *	@author Andrew Cox
 *	@version June 9, 2017
 *	@copyright GNU GPL v3.0
 */
class Arcset_cr3bp : public Arcset{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	Arcset_cr3bp(const SysData_cr3bp*);
	Arcset_cr3bp(const Arcset_cr3bp&);
	Arcset_cr3bp(const BaseArcset&);
	baseArcsetPtr create(const SysData*) const override;
	baseArcsetPtr clone() const override;
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
		// Nothing Here!
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	double getJacobi(int);
	double getJacobi(int) const;
	double getJacobiByIx(int);
	double getJacobiByIx(int) const;

	void setJacobi(int, double);
	void setJacobiByIx(int, double);
	//\}

protected:
	/**
	 *  \name File I/O
	 *  \{
	 */
		void saveCmds_toFile(mat_t*, Save_tp saveTp = Save_tp::SAVE_ALL) const override;
		void saveCmds_toStruct(matvar_t*, unsigned int, Save_tp saveTp = Save_tp::SAVE_ALL) const override;
		void readCmds_fromFile(mat_t*, std::vector<ControlLaw*>&) override;
		void readCmds_fromStruct(matvar_t*, unsigned int, std::vector<ControlLaw*>&) override;
	//\}
};

}// End of astrohelion namespace