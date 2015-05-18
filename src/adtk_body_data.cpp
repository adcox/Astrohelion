/**
 *	@file adtk_body_data.cpp
 *
 *	Data class that holds data about orbitRad body
 */
#include "adtk_body_data.hpp"

#include "adtk_constants.hpp"

#include <exception>
#include <iostream>
#include <math.h>


using namespace std;

/*
*********** CONSTRUCTOR FUNCTIONS *******************************
*/

/**
 * Default constructor; sets all fields to 0
 */
adtk_body_data::adtk_body_data(){
	radius = 0;
	mass = 0;
	gravParam = 0;
	orbitRad = 0;
	minFlyByAlt = 0;
	name = "Unknown";
	parent = "Unknown";
}//------------------------

/**
 *	Construct a body by specifying its name; the other parameters are loaded from memory
 *	@param n the name of the body; NOTE: USE ALL LOWERCASE
 */
adtk_body_data::adtk_body_data(std::string n){
	name = n;
	minFlyByAlt = 0;

	// Begin the behomoth of if statements...
	if(n.compare("earth") == 0){
		name = "Earth";
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
		parent = "N/A";
		gravParam = 1.3271244004193938e11;
		radius = 695990;
		orbitRad = 0;
		mass = gravParam/G;
		return;
    }
    
    if(n.compare("moon") == 0){
		name = "Moon";
		parent = "Earth";
		gravParam = 4.9028000661637961e3;
		radius = 1737.40;
		orbitRad = 3.850003793780e5;
		mass = gravParam/G;
		return;
    }

    if(n.compare("mercury") == 0){
		name = "Mercury";
		parent = "Sun";
		gravParam = 2.2031780000000021e4;
		radius = 2439.7;
		orbitRad = 5.913352656662e7;
		mass = gravParam/G;
		return;
    }
    if(n.compare("venus") == 0){
		name = "Venus";
		parent = "Sun";
		gravParam = 3.2485859200000006e5;
		radius = 6051.8;
		orbitRad = 1.082112582876e8;
		mass = gravParam/G;
		return;
    }
    
    if(n.compare("mars") == 0){
		name = "Mars";
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
		parent = "Sun";
		gravParam = 3.793120749865224e7;
		radius = 60268.00;
		orbitRad = 1.429737744187e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("uranus") == 0){
		name = "Uranus";
		parent = "Sun";
		gravParam = 5.793951322279009e6;
		radius = 25559;
		orbitRad = 2.873246924686e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("neptune") == 0){
		name = "Neptune";
		parent = "Sun";
		gravParam = 6.835099502439672e6;
		radius = 25559;
		orbitRad = 4.499806241072e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("pluto") == 0){
		name = "Pluto";
		parent = "Sun";
		gravParam = 8.696138177608748e2;
		radius = 1195.00;
		orbitRad = 6.183241717355e9;
		mass = gravParam/G;
		return;
    }
    if(n.compare("ganymede") == 0){
		name = "Ganymede";
		parent = "Sun";
		gravParam = 9.887834453334144e3;
		radius = 2634.1;
		orbitRad = 1070400;
		mass = gravParam/G;
		return;
    }
    if(n.compare("triton") == 0){
		name = "Triton";
		parent = "Neptune";
		gravParam = 1.427598140725034e3;
		radius = 1353.4;
		orbitRad = 354759;
		mass = gravParam/G;
		return;
    }
    if(n.compare("titan") == 0){
		name = "Titan";
		parent = "Saturn";
		gravParam = 8.978138845307376e3;
		radius = 2576;
		orbitRad = 1221870;
		mass = gravParam/G;
		return;
    }
    if(n.compare("charon") == 0){
		name = "Charon";
		parent = "Pluto";
		radius = 603.5;
		orbitRad = 17536;
		gravParam = 1.058799888601881e2;
		mass = gravParam/G;
		return;
    }
    if(n.compare("europa") == 0){
		name = "Europa";
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
    cout << "adtk_body_data constructor :: Could not locate body " << n << "\n";
    throw;
}// END of constructor using body name -------------------------------------

/**
* Constructor
*	@param m mass, kg
*	@param R mean orbital radius, km
*	@param r mean radius, km
*	@param mu gravitational parameter, km^3/s^2
* 	@param n name of the body
*	@param p name of the body"s parent
*/
adtk_body_data::adtk_body_data(double m, double R, double r, double mu, std::string n, std::string p){
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
*	@return the mean radius of the body, km
*/
double adtk_body_data::getRadius(){return radius;}

/**
*	@return the mass of the body, kg
*/
double adtk_body_data::getMass(){return mass;}

/**
*	@return the gravitational parameter for the body, km^3/s^2
*/
double adtk_body_data::getGravParam(){return gravParam;}

/**
*	@return the mean orbital radius of this body, km
*/
double adtk_body_data::getOrbitRad(){return orbitRad;}

/**
* 	@return the minimum fly-by altitude for this body, km
*/
double adtk_body_data::getMinFlyBy(){return minFlyByAlt;}

/**
*	@return the name of the body
*/
std::string adtk_body_data::getName(){return name;}

/**
*	@return the name of the parent body. If there is no parent, "None" is returned
*/
std::string adtk_body_data::getParent(){return parent;}

/**
*	Set the mean radius of the body
*	@param r the radius of the body, km
*/
void adtk_body_data::setRadius(double r){radius = r;}

/**
*	Set the mass of the body
*	@param m the mass of the body, kg
*/
void adtk_body_data::setMass(double m){mass = m;}

/**
*	Set the mean orbital radius of the body
*	@param R the orbital radius of the body, km
*/
void adtk_body_data::setOrbitRad(double R){orbitRad = R;}

/**
*	Set the gravitational parameter of the body
*	@param mu the gravitational parameter of the body, km^3/s^2
*/
void adtk_body_data::setGravParam(double mu){gravParam = mu;}

/**
*	Set the name of the body
*	@param s the name the body
*/
void adtk_body_data::setName(std::string s){name = s;}

/**
*	Set the mean radius of the body
*	@param s the name of this body"s parent. For example, Earth"s parent body is the Sun,
*	and the Moon"s parent body is Earth
*/
void adtk_body_data::setParent(std::string s){parent = s;}
