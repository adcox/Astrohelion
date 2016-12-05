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
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "BodyData.hpp"

#include "Common.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <cmath>
#include <exception>

namespace astrohelion{
/*
 *********** CONSTRUCTOR FUNCTIONS *******************************
 */

/**
 * @brief Default constructor; sets all fields to 0
 */
BodyData::BodyData(){
	radius = 0;
	mass = 0;
	gravParam = 0;
	orbitRad = 0;
	minFlyByAlt = 0;
	name = "Unknown";
	ID = 999999;
	parent = "Unknown";
}//====================================================

/**
 *	@brief Construct a body by specifying its name; the other parameters are loaded from memory
 *	@param n the name of the body; NOTE: USE ALL LOWERCASE
 *	@throws Exception if the body cannot be located
 */
BodyData::BodyData(std::string n){
	name = n;
	minFlyByAlt = 0;

	transform(n.begin(), n.end(), n.begin(), ::tolower);
	
	// Begin the behomoth of if statements...
	if(n.compare("earth") == 0){
		name = "Earth";
		ID = 399;
        parent = "Sun";
        gravParam = 3.9860043543609593e5;
        radius = 6378.137;
        orbitRad = 1.4959784234299874e8;	// Average SMA from frequency analysis of SPICE data between 1700 and 2200
        mass = gravParam/G;
        // minFlyByAlt = 2000;
        return;
	}

	if(n.compare("earth_barycenter") == 0){
		name = "Earth Barycenter";
		ID = 3;
		parent = "Sun";
		gravParam = 4.0350323550225981e5;
		radius = 0;
		orbitRad = 1.4959784418519253e8;	// Average SMA from frequency analysis of SPICE data between 1700 and 2200
		mass = gravParam/G;
		return;
	}

	if(n.compare("sun") == 0){
		name = "Sun";
		ID = 10;
		parent = "N/A";
		gravParam = 1.3271244004193930e11;	// Average SMA from frequency analysis of SPICE data between 1700 and 2200
		radius = 695990;
		orbitRad = 0;
		mass = gravParam/G;
		return;
    }
    
    if(n.compare("moon") == 0){
		name = "Moon";
		ID = 301;
		parent = "Earth";
		gravParam = 4.9028000661637961e3;
		radius = 1737.40;
		orbitRad = 3.8474793244470679e5;	// Average SMA from frequency analysis of SPICE data between 1700 and 2200
		mass = gravParam/G;
		return;
    }

    if(n.compare("mercury") == 0){
		name = "Mercury";
		ID = 199;
		parent = "Sun";
		gravParam = 2.2031780000000021e4;
		radius = 2439.7;
		orbitRad = 5.7909137745690152e7;	// Average SMA from frequency analysis of SPICE data between 1700 and 2200
		mass = gravParam/G;
		return;
    }
    if(n.compare("venus") == 0){
		name = "Venus";
		ID = 299;
		parent = "Sun";
		gravParam = 3.248585920000000e5;
		radius = 6051.8;
		orbitRad = 1.0820890093249036e8;
		mass = gravParam/G;
		return;
    }
    
    if(n.compare("mars") == 0){
		name = "Mars";
		ID = 499;
		parent = "Sun";
		gravParam = 4.282837362069909e4;
		radius = 3396.19;
		orbitRad = 2.289352472157e8;
		mass = gravParam/G;
		return;
    }
    if(n.compare("jupiter") == 0){
		name = "Jupiter";
		parent = "Sun";
		gravParam = 1.266865349218008e8;
		radius = 71492.00;
		orbitRad = 7.793390706261e9;
		mass = 1.898E27;
		return;
    }
    if(n.compare("saturn") == 0){
		name = "Saturn";
		ID = 699;
		parent = "Sun";
		gravParam = 3.793120749865224e7;
		radius = 60268.00;
		orbitRad = 1.429737744187e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("uranus") == 0){
		name = "Uranus";
		ID = 799;
		parent = "Sun";
		gravParam = 5.793951322279009e6;
		radius = 25559;
		orbitRad = 2.873246924686e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("neptune") == 0){
		name = "Neptune";
		ID = 899;
		parent = "Sun";
		gravParam = 6.835099502439672e6;
		radius = 25559;
		orbitRad = 4.499806241072e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("pluto") == 0){
		name = "Pluto";
		ID = 999;
		parent = "Sun";
		gravParam = 8.696138177608748e2;
		radius = 1195.00;
		orbitRad = 6.183241717355e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("ganymede") == 0){
		name = "Ganymede";
		ID = 503;
		parent = "Sun";
		gravParam = 9.887834453334144e3;
		radius = 2634.1;
		orbitRad = 1070400;
		mass = gravParam/G;
		return;
    }
    if(n.compare("triton") == 0){
		name = "Triton";
		ID = 801;
		parent = "Neptune";
		gravParam = 1.427598140725034e3;
		radius = 1353.4;
		orbitRad = 354759;
		mass = gravParam/G;
		return;
    }
    if(n.compare("titan") == 0){
		name = "Titan";
		ID = 606;
		parent = "Saturn";
		gravParam = 8.978138845307376e3;
		radius = 2576;
		orbitRad = 1221870;
		mass = gravParam/G;
		return;
    }
    if(n.compare("charon") == 0){
		name = "Charon";
		ID = 901;
		parent = "Pluto";
		radius = 603.5;
		orbitRad = 17536;
		gravParam = 1.058799888601881e2;
		mass = gravParam/G;
		return;
    }
    if(n.compare("europa") == 0){
		name = "Europa";
		ID = 502;
		parent = "Jupiter";
		radius = 1560;
		orbitRad = 671101.9638;
		gravParam = 3.202738774922892e3;
		mass = gravParam/G;
		return;
    }
    if(n.compare("binary") == 0){
		name = "Binary";
		parent = "Binary";
		gravParam = 1;
		radius = 1;
		orbitRad = 40;
		mass = 1;
		return;
    }

    // Function should have returned by this point if the body was found
    throw Exception("Could not locate body");
}// END of constructor using body name -------------------------------------

/**
 *  @brief Constructor
 *	@param m mass, kg
 *	@param R mean orbital radius, km
 *	@param r mean radius, km
 *	@param mu gravitational parameter, km^3/s^2
 * 	@param n name of the body
 *	@param p name of the body"s parent
 */
BodyData::BodyData(double m, double R, double r, double mu, std::string n, std::string p){
	mass = m;
	orbitRad = R;
	radius = r;
	gravParam = mu;
	name = n;
	parent = p;
}//-------------------------

/*
 *********** SET AND GET FUNCTIONS *******************************
 */

/**
 *	@brief Retrieve the mean radius of the body, km
 *	@return the mean radius of the body, km
 */
double BodyData::getRadius(){return radius;}

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
 *	@return the mean orbital radius of this body, km
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
 *	@brief Set the mean radius of the body
 *	@param r the radius of the body, km
 */
void BodyData::setRadius(double r){radius = r;}

/**
 *	@brief Set the mass of the body
 *	@param m the mass of the body, kg
 */
void BodyData::setMass(double m){mass = m;}

/**
 *	@brief Set the mean orbital radius of the body
 *	@param R the orbital radius of the body, km
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
 *	@brief Set the mean radius of the body
 *	@param s the name of this body"s parent. For example, Earth"s parent body is the Sun,
 *	and the Moon"s parent body is Earth
 */
void BodyData::setParent(std::string s){parent = s;}

}// END of Astrohelion namespace
