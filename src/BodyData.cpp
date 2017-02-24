/**
 *	@file BodyData.cpp
 *
 *	@brief Data class that holds data about orbitRad body
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
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

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#include "BodyData.hpp"

#include "Common.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

namespace astrohelion{
/*
 *********** CONSTRUCTOR FUNCTIONS *******************************
 */

/**
 * @brief Default constructor; sets all fields to 0
 */
BodyData::BodyData(){
	bodyRad = 0;
	mass = 0;
	gravParam = 0;
	orbitRad = 0;
	minFlyByAlt = 0;
	name = "Unknown";
	ID = 999999;
	parent = "Unknown";
}//====================================================

/**
 *  \brief Construct a BodyData object from preloaded data.
 *  \details Data about a body is loaded from a map that pulls data
 *  from the body_data.xml file.
 * 
 *  \param n name of the body or dynamical object (e.g., a system barycenter).
 *  \throws Exception if the body with the specified name does not have a SPICE
 *  ID or if the body has not been loaded
 */
BodyData::BodyData(std::string n){
	// Get the ID associated with the name; throws an exception if it can't find one
	ID = getSpiceIDFromName(n.c_str());

	// Get a reference to the loaded body data list and make sure the body we're looking for is there
	std::map<int, Body_Data> &allBodyData = this->initializer().allBodyData;
	if(allBodyData.find(ID) == allBodyData.end()){
		throw Exception("BodyData::BodyData: Body has not been loaded\n");
	}

	// Load parameters from the list
	bodyRad = allBodyData[ID].bodyRad;
	gravParam = allBodyData[ID].gravParam;
	orbitRad = allBodyData[ID].orbitRad;
	name = allBodyData[ID].name;
	parent = allBodyData[ID].parent;
	mass = gravParam/G;
}//====================================================

/**
 *  @brief Constructor
 *	@param m mass, kg
 *	@param R mean orbital bodyRad, km
 *	@param r mean bodyRad, km
 *	@param mu gravitational parameter, km^3/s^2
 * 	@param n name of the body
 *	@param p name of the body"s parent
 */
BodyData::BodyData(double m, double R, double r, double mu, std::string n, std::string p){
	mass = m;
	orbitRad = R;
	bodyRad = r;
	gravParam = mu;
	name = n;
	parent = p;
}//-------------------------

/*
 *********** SET AND GET FUNCTIONS *******************************
 */

/**
 *	@brief Retrieve the mean bodyRad of the body, km
 *	@return the mean bodyRad of the body, km
 */
double BodyData::getBodyRad(){return bodyRad;}

/**
 *	@brief Retrieve the mass of the body, kg
 *	@return the mass of the body, kg
 */
double BodyData::getMass(){return mass;}

/**
 *	@brief Retrieve the gravitational parameter for the body, km^3/s^2
 *	@return the gravitational parameter for the body, km^3/s^2
 */
double BodyData::getGravParam(){return gravParam;}

/**
 *	@brief Retrieve
 *	@return the mean orbital bodyRad of this body, km
 */
double BodyData::getOrbitRad(){return orbitRad;}

/**
 *	@brief Retrieve the minimum fly-by altitude for this body, km
 * 	@return the minimum fly-by altitude for this body, km
 */
double BodyData::getMinFlyBy(){return minFlyByAlt;}

/**
 *	@brief Retrieve the name of the body
 *	@return the name of the body
 */
std::string BodyData::getName(){return name;}

/**
 *	@brief Retrieve the ID (SPICE or HORIZONS ID) associated with this body
 *	@return the ID associated with this body
 */
int BodyData::getID(){ return ID; }

/**
 *	@brief Retrieve the name of the parent body
 *	@return the name of the parent body. If there is no parent, "None" is returned
 */
std::string BodyData::getParent(){return parent;}

/**
 *	@brief Set the mean bodyRad of the body
 *	@param r the bodyRad of the body, km
 */
void BodyData::setBodyRad(double r){bodyRad = r;}

/**
 *	@brief Set the mass of the body
 *	@param m the mass of the body, kg
 */
void BodyData::setMass(double m){mass = m;}

/**
 *	@brief Set the mean orbital bodyRad of the body
 *	@param R the orbital bodyRad of the body, km
 */
void BodyData::setOrbitRad(double R){orbitRad = R;}

/**
 *	@brief Set the gravitational parameter of the body
 *	@param mu the gravitational parameter of the body, km^3/s^2
 */
void BodyData::setGravParam(double mu){gravParam = mu;}

/**
 *	@brief Set the name of the body
 *	@param s the name the body
 */
void BodyData::setName(std::string s){name = s;}

/**
 *	@brief Set the mean bodyRad of the body
 *	@param s the name of this body"s parent. For example, Earth"s parent body is the Sun,
 *	and the Moon"s parent body is Earth
 */
void BodyData::setParent(std::string s){parent = s;}

void BodyData::print(){
	printf("Body data for %s:\n", name.c_str());
	printf("  SPICE ID = %d\n", ID);
	printf("  Parent = %s\n", parent.c_str());
	printf("  GM = %.6e\n", gravParam);
	printf("  Circ. Orbit Rad = %.6e km\n", orbitRad);
	printf("  Body Rad = %.6e km\n", bodyRad);
	printf("  Mass = %.6e kg\n", mass);
}//====================================================

}// END of Astrohelion namespace
