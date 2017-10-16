/**
 * \file PseudoArcEngine.hpp
 * \brief
 * 
 * \author Andrew Cox
 * \version October 13, 2017
 * \copyright GNU GPL v3.0
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

#include "ContinuationEngine.hpp"
#include "EigenDefs.hpp"

namespace astrohelion{

// Forward declarations
class Arcset_cr3bp;
class Arcset_periodic;
class Family_PO;
enum class Mirror_tp;
class MultShootData;

class PseudoArcEngine : public ContinuationEngine{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	PseudoArcEngine();
	PseudoArcEngine(const PseudoArcEngine&);
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	PseudoArcEngine& operator=(const PseudoArcEngine&);
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	double getStepSize();

	void setStepSize(double);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void generatePO_cr3bp(Family_PO*, const Arcset_cr3bp*, Mirror_tp, std::vector<int>);
	void getNextPACGuess_cr3bp(Arcset_cr3bp*, const Eigen::VectorXd&, const Eigen::VectorXd&, double, MultShootData*);
	bool checkPACSolution_cr3bp(const Arcset_periodic*, const MultShootData*, const Eigen::VectorXd&, double, bool*);
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	void reset() override;
	//\}

private:
	/**
	 *  \name Utility Functions
	 *  \{
	 */
	void copyMe(const PseudoArcEngine&);
	void cleanEngine() override;
	//\}

	double stepSize = 0.001;
};

}// End of astrohelion namespace