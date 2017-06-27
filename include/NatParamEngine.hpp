/**
 * \file NatParamEngine.hpp
 * \brief
 * 
 * \author Andrew Cox
 * \version June 27, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

namespace astrohelion{


class NatParamEngine : public ContinuationEngine{
public:
	NatParamEngine();

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	void setSlopeThresh(double);
	void setStep_simple(double);
	void setStep_fitted_1(double);
	void setStep_fitted_2(double);
	//\}

	/**
	 *  \name Continuation Functions
	 *  \{
	 */
	void continuePeriodic_cr3bp(Fam_cr3bp*, const Arcset_cr3bp*, std::vector<Mirror_tp>, std::vector<int>, std::vector<int>, int order = 1);
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	void copyMe(const NatParamEngine&);
	void reset();
	void cleanEngine();
	//\}

private:
	int curveFitMem = 5;			//!< Number of points to use with Least-Squares algorithm
	int numSimple = 3;				//!< Number of simply-continued family members
	double slopeThresh = 1;			//!< Minimum slope for stepping in indVar1; else step in indVar2
	double step_simple = 5e-4;		//!< Step size in the independent variable when using simple continuation
	double step_fitted_1 = 5e-3;	//!< Step size in first ind. var. when using advanced continuation
	double step_fitted_2 = 5e-3;	//!< Step size in second ind. var. when using advanced continuation
};

}// end of astrohelion namespace