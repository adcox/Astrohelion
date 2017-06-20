/**
 *  \file Engine.hpp
 *  \brief Contains the base behavior for engine objects
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

#include "Common.hpp"

namespace astrohelion{ 
// Forward Declarations

/**
 *	\ingroup engine
 *	\brief Contains the base behavior for engine objects
 *
 *	\author Andrew Cox
 *	\version September 9, 2016
 *	\copyright GNU GPL v3.0
 */
class Engine{

public:
	
	virtual ~Engine(){}

	/**
	 *  \brief Retrieve the verbosity (i.e., how many messages will be printed)
	 *  of the Engine object
	 *  \return The verbosity of the engine
	 */
	Verbosity_tp getVerbosity() const{
		return verbosity;
	}//================================================

	/**
	 *  \brief Set the verbosity (i.e., how many messages will be printed)
	 *  of the Engine object
	 * 
	 *  \param v The verbosity of the engine
	 */
	void setVerbosity(Verbosity_tp v){
		verbosity = v;
	}//================================================
	
	/**
	 *  \brief Resets the engine, including any parameters the user
	 *  may have set
	 */
	virtual void reset() = 0;

protected:

	/**
	 *  \brief Make this Engine object the same as the input Engine
	 *  \details This function is used in copy constructors to avoid
	 *  code duplication
	 * 
	 *  \param e Reference to an Engine object
	 */
	void copyBaseEngine(const Engine &e){
		verbosity = e.verbosity;
		bIsClean = e.bIsClean;
	}//================================================

	/**
	 *  \brief Reset data storage parts of the engine, but not parameters the user
	 *  has set (i.e., tolerances, etc.)
	 */
	virtual void cleanEngine() = 0;

	/** Describes the number of messages the engine should output by default */
	Verbosity_tp verbosity = Verbosity_tp::SOME_MSG;

	/** Whether or not the engine is "clean"; if it isn't data is being stored that may be cleaned out
		by the cleanEngine() function */
	bool bIsClean = true;

private:

};

}// END of Astrohelion namespace