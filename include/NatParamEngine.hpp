/**
 * \file NatParamEngine.hpp
 * \brief
 * 
 * \author Andrew Cox
 * \version June 27, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

#include <vector>

#include "ContinuationEngine.hpp"

namespace astrohelion{

// Forward Declarations
class Family_PO;
class Arcset_cr3bp;
class Arcset_periodic;
enum class Mirror_tp;

/**
 *  \ingroup engine
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
	void generateSymmetricPO_cr3bp(Family_PO*, const Arcset_cr3bp*, std::vector<Mirror_tp>, std::vector<unsigned int>, std::vector<unsigned int>, int order = 1);
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
	unsigned int curveFitMem = 5;			//!< Number of points to use with Least-Squares algorithm
	unsigned int numSimple = 3;				//!< Number of simply-continued family members
	double slopeThresh = 1;			//!< Minimum slope for stepping in indVar1; else step in indVar2
	double step_simple = 5e-4;		//!< Step size in the independent variable when using simple continuation
	double step_fitted_1 = 5e-3;	//!< Step size in first ind. var. when using advanced continuation
	double step_fitted_2 = 5e-3;	//!< Step size in second ind. var. when using advanced continuation
};

}// end of astrohelion namespace