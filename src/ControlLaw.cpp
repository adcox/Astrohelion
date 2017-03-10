/**
 *  \file ControlLaw.cpp
 *  \author Andrew Cox
 *  \version March 3, 2017
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

#include "ControlLaw.hpp"
#include "Exceptions.hpp"
#include "SysData.hpp"

namespace astrohelion{

ControlLaw::ControlLaw(){}

//------------------------------------------------------------------------------------------------------
//      Switchboard Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param t [description]
 *  \param s [description]
 *  \param pSys [description]
 *  \param int [description]
 *  \param law [description]
 *  \param int [description]
 */
void ControlLaw::getLaw(double t, const double *s, const SysData *pSys, unsigned int lawID, double *law, unsigned int len) const{
	switch(lawID){
		case NO_CTRL:
			// Handle default case with no control
			for(unsigned int i = 0; i < len; i++){
				law[i] = 0;
			}
			break;
		default:
			throw Exception("ControlLaw::GetLaw: Unrecognized lawID");
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

/**
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param t [description]
 *  \param s [description]
 *  \param pSys [description]
 *  \param int [description]
 *  \param partials [description]
 *  \param int [description]
 */
void ControlLaw::getPartials_State(double t, const double *s, const SysData *pSys, unsigned int lawID, double *partials, unsigned int len) const{
	switch(lawID){
		case NO_CTRL:
			// Handle default case with no control
			for(unsigned int i = 0; i < len; i++){
				partials[i] = 0;
			}
			break;
		default:
			throw Exception("ControlLaw::GetPartials_State: Unrecognized lawID");
	}
	(void) t;
	(void) s;
	(void) pSys;
}//====================================================

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

/**
 *  \brief [brief description]
 *  \details [long description]
 * 
 *  \param int [description]
 *  \return [description]
 */
std::string ControlLaw::lawIDToString(unsigned int id) const{
	switch(id){
		case NO_CTRL: return "NONE";
		default: return "UNDEFINED";
	}
}//====================================================

}