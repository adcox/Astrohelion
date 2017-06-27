/**
 * \file ContinuationEngine.hpp
 * \brief
 * 
 * \author Andrew Cox
 * \version June 27, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

#include "Core.hpp"

namespace astrohelion{

// Forward Declarations


/**
 *	\brief Type of continuation to use when generating a family
 */
enum class Continuation_tp{
	NAT_PARAM,	//!< Use natural parameter continuation
	PSEUDO_ARC	//!> Use pseudo arclength continuation
};

class ContinuationEngine : public Core, public Engine{
	public:
		/**
	 	* \name *structors
	 	* \{
	 	*/
		ContinuationEngine(){}
		
		ContinuationEngine(const ContinuationEngine &engine){
			copyMe(engine);
			return *this;
		}//============================================

		~ContinuationEngine(){}
		//\}

		void setMaxStepSize(double step){ maxStepSize = step; }
		void setMinStepSize(double step){ minStepSize = step; }
		void setNumNodes(int num){ numNodes = num; }
		void setNumOrbits(int num){ numOrbits = num; }
		void setTol(double tol){ this->tol = tol; }

		void reset(){
			tol = 1e-12;
			numOrbits = 500;
			numNodes = 3;
		}//============================================

	private:
		
		double minStepSize = 1e-6;		//!< Minimum allowable step size
		double maxStepSize = 0.05;		//!< Maximum allowable step size
		double numOrbits = 500;			//!< Maximum number of orbits in the continuation
		double numNodes = 3;			//!< Number of nodes to use on each trajectory
		double tol = 1e-12;				//!< Tolerance for corrections

		void copyMe(const ContinuationEngine &engine){
			tol = engine.tol;
			numOrbits = engine.numOrbits;
			numNodes = engine.numNodes;
		}//============================================

		void cleanEngine(){}
};

}// End of astrohelion namespace




