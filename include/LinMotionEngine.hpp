/**
 *  @file LinMotionEngine.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version September 28, 2017
 *	@copyright GNU GPL v3.0
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

#include "Core.hpp"
#include "Engine.hpp"

namespace astrohelion{

/**
 *	\ingroup engine
 *	@brief type of linear motion
 *
 *	This type tells the engine what kind of linearization to produce
 */
class LinMotion_tp{
public:
	static const unsigned int NONE = 0;			//!< No motion specified
	static const unsigned int HYP = 10;			//!< Stable/Unstable mode (real eigenvalues)
	static const unsigned int OSC = 20;			//!< Oscillatory mode (imaginary eigenvalues)
	static const unsigned int STAB_OSC = 30; 	//!< Stable Oscillatory (complex with negative real part)
	static const unsigned int UNSTAB_OSC = 40;	//!< Unstable Oscillatory (complex with positive real part)
};

/**
 *	@brief An engine that will generate a trajectory in the linearized EOMs
 *	@details This class is merely a shell for the guts that are model-specific
 *	
 *	@see LinMotionEngine_cr3bp
 */
class LinMotionEngine : public Core, public Engine{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */

		/**
		 *	@brief Default, do-nothing constructor
		 */
		LinMotionEngine() {}

		/**
		 *  @brief Default, do-nothing desctructor
		 */
		virtual ~LinMotionEngine() {}

		//\}

		// Operators
		//  - None!

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		double getNumRevs() const;
		double getTimeStep() const;
		double getTol() const;

		void setNumRevs(double);
		void setTimeStep(double);
		void setTol(double);
		virtual const char* getTypeStr(unsigned int) const;
		//\}
		
	protected:
		/** @brief step size between points on linear motion trajectory */
		double t_step = 0.001;

		/** @brief Number of revolutions to propagate */
		double revs = 1;

		/** @brief tolerance for numerical methods, like locating Lagrange points */
		double tol = 1e-14;

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		virtual void cleanEngine();
		virtual void reset();
		//\}
};

}// END of Astrohelion namespace