/**
 * \file NatParamEngine.hpp
 * \brief
 * 
 * \author Andrew Cox
 * \version June 27, 2017
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

#include <vector>

#include "ContinuationEngine.hpp"

namespace astrohelion{

// Forward Declarations
class Arcset;
class Arcset_cr3bp;
class Arcset_periodic;
class Family_PO;
enum class Mirror_tp;

/**
 *  \ingroup engine fam
 *  \brief Perform natural parameter continuation on various types of structures
 *  in different models
 */
class NatParamEngine : public ContinuationEngine{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	NatParamEngine();
	NatParamEngine(const NatParamEngine&);
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	NatParamEngine& operator=(const NatParamEngine&);
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	void setCurveFitMem(unsigned int);
	void setNumSimple(unsigned int);
	void setSlopeThresh(double);
	void setStep_simple(double);
	void setStep_fitted_1(double);
	void setStep_fitted_2(double);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void continueSymmetricPO_cr3bp(Family_PO*, const Arcset_cr3bp*, std::vector<Mirror_tp>, std::vector<unsigned int>, std::vector<unsigned int>);
	void continuePO_cr3bp(Family_PO*, const Arcset_cr3bp*, std::vector<double>, std::vector<unsigned int>, std::vector<unsigned int>);
	void continuePO(const Arcset*, Arcset*, Arcset*, std::vector<Arcset>&, std::vector<double>, std::vector<unsigned int>, std::vector<unsigned int>);
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
	void copyMe(const NatParamEngine&);
	void cleanEngine() override;
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	bool decreaseStepSize();
	std::vector<double> familyCont_LS(unsigned int, double, std::vector<unsigned int>, std::vector<double>);
	void increaseStepSize(unsigned int);
	bool updateIC(const Arcset&, std::vector<double>*, double*, const std::vector<unsigned int>&,
		const std::vector<unsigned int>&, const std::vector<Arcset>&);
	//\}


	unsigned int curveFitMem = 5;			//!< Number of points to use with Least-Squares algorithm
	double deltaVar1 = 1;					//!< Change in independent variable #1 during continuation
	double deltaVar2 = 1;					//!< Change in independent variable #2 during continuation
	double indVarSlope = 1;					//!< Slope of the indVar1 vs indVar2 line
	std::vector<unsigned int> fixStates {};	//!< A list of indices that represent the states that are constrained (i.e., fixed) during the current continuation step
	unsigned int numSimple = 3;				//!< Number of simply-continued family members
	unsigned int orbitCount = 0;			//!< Count of how many orbits have been computed
	double stepScaleFactor = 2;				//!< A multiplying factor to scale the step size up/down during adaptive step sizing
	double slopeThresh = 1;			//!< Minimum slope for stepping in indVar1; else step in indVar2
	double step_simple = 5e-4;		//!< Step size in the independent variable when using simple continuation
	double step_fitted_1 = 5e-3;	//!< Step size in first ind. var. when using advanced continuation
	double step_fitted_2 = 5e-3;	//!< Step size in second ind. var. when using advanced continuation
};

}// end of astrohelion namespace