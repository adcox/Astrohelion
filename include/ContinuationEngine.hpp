/**
 * \file ContinuationEngine.hpp
 * \brief Define the basic behavior and common variables for continuation engines
 * \details This class is fully defined in the header; there is no accompanying *.cpp file.
 * 
 * \author Andrew Cox
 * \version October 3, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

#include "Core.hpp"
#include "Engine.hpp"

namespace astrohelion{

/**
 *	\brief Type of continuation to use when generating a family
 */
enum class Continuation_tp{
	NAT_PARAM,	//!< Use natural parameter continuation
	PSEUDO_ARC	//!> Use pseudo arclength continuation
};

/**
 *  \ingroup engine fam
 *  \brief An engine object that runs continuation processes
 *  \details This class contains the template to specify tolerances, 
 *  min and max step sizes, and common behaviors like reset() and 
 *  cleanEngine()
 */
class ContinuationEngine : public Core, public Engine{
	public:
		/**
	 	* \name *structors
	 	* \{
	 	*/

	 	/**
	 	 *  \brief Default constructor
	 	 */
		ContinuationEngine(){}
		
		/**
		 *  \brief Copy Constructor
		 *  \param engine reference to an engine object
		 */
		ContinuationEngine(const ContinuationEngine &engine){
			copyMe(engine);
		}//============================================

		/**
		 *  \brief Default destructor
		 */
		~ContinuationEngine(){}
		//\}

		/**
	 	* \name Operators
	 	* \{
	 	*/
		ContinuationEngine& operator=(const ContinuationEngine &c){
			copyMe(c);
			return *this;
		}//============================================
	 	//\}
	 	
		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		/**
		 *  \brief Set the maximum step size
		 * 
		 *  \param step the largest permissible step size in 
		 *  whatever parameter is leveraged in the continuation process.
		 *  Most methods leverage dynamic step-size adjustment, so
		 *  if the algorithm is converging quickly, it will increase
		 *  the step size up to this maximum amount. Thus, a smaller maximum
		 *  may yield a finer resolution in the family
		 */
		void setMaxStepSize(double step){ maxStepSize = step; }

		/**
		 *  \brief Set the minimum step size
		 * 
		 *  \param step the smallest permissible step size in 
		 *  whatever parameter is leveraged in the continuation process.
		 *  Most methods leverage dynamic step-size adjustment, so if 
		 *  convergence is not achieved with this minimum step size,
		 *  the continuation is typically abandoned.
		 */
		void setMinStepSize(double step){ minStepSize = step; }

		/**
		 *  \brief Set the number of orbits to generate
		 *  \details This is the maximum number of orbits that my be generated via
		 *  continuation. If the continuation algorithm does not locate the end of
		 *  the family (or fails to converge), the continuation process will stop
		 *  with the specified limit.
		 * 
		 *  \param num maximum number of orbits to generate via continuation
		 */
		void setNumOrbits(unsigned int num){ numOrbits = num; }

		/**
		 *  \brief Set the iteration threshold at which to increase the step size
		 *  \details If a solution is converged with fewer iterations than this quantity,
		 *  the step size is increased
		 * 
		 *  \param c the step count threshold at which to increase the step size
		 */
		void setStepCountIncrease(unsigned int c){ stepCount_increase = c; }

		/**
		 *  \brief Set the step scale factor
		 *  \details When the step size is increased, the current step size is multiplied by this value.
		 *  Conversely, decreasing the step size divides the step size by this value.
		 * 
		 *  \param s the step scale factor
		 */
		void setStepScaleFactor(double s){ stepScaleFactor = s; }

		/**
		 *  \brief Set the numerical tolerance
		 *  \details This tolerance is passed to multiple-shooting and
		 *  simulation engines.
		 * 
		 *  \param tol the numerical tolerance
		 */
		void setTol(double tol){ this->tol = tol; }
		//\}

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		/**
		 *  \brief Reset all continuation parameters to their default values
		 */
		virtual void reset(){
			tol = 1e-12;
			minStepSize = 1e-6;
			maxStepSize = 0.05;
			numOrbits = 500;
		}//============================================
		//\}

	protected:
		
		double minStepSize = 1e-6;		//!< Minimum allowable step size
		double maxStepSize = 0.05;		//!< Maximum allowable step size
		unsigned int numOrbits = 500;	//!< Maximum number of orbits in the continuation
		unsigned int stepCount_increase = 5;	//!< Increase the step size if convergence is reached in fewer than this many iterations
		double stepScaleFactor = 2;		//!< A multiplying factor to scale the step size up/down during adaptive step sizing
		double tol = 1e-12;				//!< Tolerance for corrections

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		/**
		 *  \brief Copy all member variables associated with the continuation engine
		 *  and its parent classes.
		 * 
		 *  \param engine reference to an existing engine object
		 */
		virtual void copyMe(const ContinuationEngine &engine){
			Engine::copyBaseEngine(engine);
			tol = engine.tol;
			minStepSize = engine.minStepSize;
			maxStepSize = engine.maxStepSize;
			numOrbits = engine.numOrbits;
			stepCount_increase = engine.stepCount_increase;
			stepScaleFactor = engine.stepScaleFactor;
		}//============================================

		/**
		 *  \brief Resets variables specific to a continuation process
		 *  \details This does <emph>not</emph> reset the continuation parameters
		 *  such as step size, tolerance, number of orbits, etc.
		 */
		virtual void cleanEngine(){}
		//\}
};

}// End of astrohelion namespace




