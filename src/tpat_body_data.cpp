/**
 *	@file tpat_body_data.cpp
 *
 *	@brief Data class that holds data about orbitRad body
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "tpat.hpp"

#include "tpat_body_data.hpp"

#include "tpat_constants.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

#include <algorithm>
#include <cmath>
#include <exception>

/*
 *********** CONSTRUCTOR FUNCTIONS *******************************
 */

/**
 * @brief Default constructor; sets all fields to 0
 */
tpat_body_data::tpat_body_data(){
	radius = 0;
	mass = 0;
	gravParam = 0;
	orbitRad = 0;
	minFlyByAlt = 0;
	name = "Unknown";
	ID = 999999;
	parent = "Unknown";
}//------------------------

/**
 *	@brief Construct a body by specifying its name; the other parameters are loaded from memory
 *	@param n the name of the body; NOTE: USE ALL LOWERCASE
 */
tpat_body_data::tpat_body_data(std::string n){
	name = n;
	minFlyByAlt = 0;

	transform(n.begin(), n.end(), n.begin(), ::tolower);
	
	// Begin the behomoth of if statements...
	if(n.compare("earth") == 0){
		name = "Earth";
		ID = 399;
        parent = "Sun";
        gravParam = 3.9860043543609598e5;
        radius = 6378.137;
        orbitRad = 1.496188026821e8;
        mass = gravParam/G;
        minFlyByAlt = 2000;
        return;
	}

	if(n.compare("sun") == 0){
		name = "Sun";
		ID = 10;
		parent = "N/A";
		gravParam = 1.3271244004193938e11;
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
		orbitRad = 3.850003793780e5;
		mass = gravParam/G;
		return;
    }

    if(n.compare("mercury") == 0){
		name = "Mercury";
		ID = 199;
		parent = "Sun";
		gravParam = 2.2031780000000021e4;
		radius = 2439.7;
		orbitRad = 5.913352656662e7;
		mass = gravParam/G;
		return;
    }
    if(n.compare("venus") == 0){
		name = "Venus";
		ID = 299;
		parent = "Sun";
		gravParam = 3.2485859200000006e5;
		radius = 6051.8;
		orbitRad = 1.082112582876e8;
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
    throw tpat_exception("Could not locate body");
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
tpat_body_data::tpat_body_data(double m, double R, double r, double mu, std::string n, std::string p){
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
double tpat_body_data::getRadius(){return radius;}

/**
 *	@brief Retrieve the mass of the body, kg
 *	@return the mass of the body, kg
 */
double tpat_body_data::getMass(){return mass;}

/**
 *	@brief Retrieve the gravitational parameter for the body, km^3/s^2
 *	@return the gravitational parameter for the body, km^3/s^2
 */
double tpat_body_data::getGravParam(){return gravParam;}

/**
 *	@brief Retrieve
 *	@return the mean orbital radius of this body, km
 */
double tpat_body_data::getOrbitRad(){return orbitRad;}

/**
 *	@brief Retrieve the minimum fly-by altitude for this body, km
 * 	@return the minimum fly-by altitude for this body, km
 */
double tpat_body_data::getMinFlyBy(){return minFlyByAlt;}

/**
 *	@brief Retrieve the name of the body
 *	@return the name of the body
 */
std::string tpat_body_data::getName(){return name;}

/**
 *	@brief Retrieve the ID (SPICE or HORIZONS ID) associated with this body
 *	@return the ID associated with this body
 */
int tpat_body_data::getID(){ return ID; }

/**
 *	@brief Retrieve the name of the parent body
 *	@return the name of the parent body. If there is no parent, "None" is returned
 */
std::string tpat_body_data::getParent(){return parent;}

/**
 *	@brief Set the mean radius of the body
 *	@param r the radius of the body, km
 */
void tpat_body_data::setRadius(double r){radius = r;}

/**
 *	@brief Set the mass of the body
 *	@param m the mass of the body, kg
 */
void tpat_body_data::setMass(double m){mass = m;}

/**
 *	@brief Set the mean orbital radius of the body
 *	@param R the orbital radius of the body, km
 */
void tpat_body_data::setOrbitRad(double R){orbitRad = R;}

/**
 *	@brief Set the gravitational parameter of the body
 *	@param mu the gravitational parameter of the body, km^3/s^2
 */
void tpat_body_data::setGravParam(double mu){gravParam = mu;}

/**
 *	@brief Set the name of the body
 *	@param s the name the body
 */
void tpat_body_data::setName(std::string s){name = s;}

/**
 *	@brief Set the mean radius of the body
 *	@param s the name of this body"s parent. For example, Earth"s parent body is the Sun,
 *	and the Moon"s parent body is Earth
 */
void tpat_body_data::setParent(std::string s){parent = s;}
