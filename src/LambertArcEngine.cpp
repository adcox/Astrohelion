/**
 *  \file LambertArcEngine.cpp
 *	\brief Generate Lambert Arcs in the 2BP
 *
 *	\author Andrew Cox
 *	\version 
 *	\copyright GNU GPL v3.0
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

#include "LambertArcEngine.hpp"

#include "Exceptions.hpp" 
#include "SysData_2bp.hpp"
#include "Traj_2bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

LambertArcEngine::LambertArcEngine(){}

LambertArcEngine::LambertArcEngine(const LambertArcEngine &e){
	copyMe(e);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

LambertArcEngine& LambertArcEngine::operator =(const LambertArcEngine &e){
	copyMe(e);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Set the maximum allowable TOF error, seconds
 *  @details Default value is 1 millisecond
 * 
 *  \param err maximum allowable TOF error, seconds
 */
void LambertArcEngine::setMaxErr_tof(double err){ tof_maxErr = err; }

/**
 *  \brief Set the maximum allowable true anomaly error, radians
 *  @details Default value is 0.001 deg equivalent
 * 
 *  \param err maximum allowable true anomaly error, radians
 */
void LambertArcEngine::setMaxErr_ta(double err){ ta_maxErr = err; }

/**
 *  \brief Set the maximum allowable iterations to find the semi-major axis
 *  @details Default value is 100
 * 
 *  \param its maximum allowable iterations to find the semi-major axis
 */
void LambertArcEngine::setMaxIts(unsigned int its){ maxIts = its; }

//-----------------------------------------------------
//      Lambert Arc Generation Functions
//-----------------------------------------------------

Traj_2bp LambertArcEngine::getLambertArc(SysData_2bp *pSys, std::vector<double> r1, std::vector<double> r2, double tof, unsigned int type){
	if(r1.size() < 3 || r2.size() < 3)
        throw Exception("r2bp_getLambertArc: r1 or r2 has fewer than three elements");

    Traj_2bp traj(pSys);
    std::vector<double> state1(3,0), state2(3,0);   // three zero velocity components

    // Copy the positions into two state vectors; velocity components are initially all set to zero
    state1.insert(state1.begin(), r1.begin(), r1.end());
    state2.insert(state2.begin(), r2.begin(), r2.end());
    
    // Create nodes to contain the states
    Node node1(state1, 0);
    Node node2(state2, tof);

    // Define some step sizes and tolerances for the iteration scheme
    double a_step = pSys->getMu()/(2e5);    // Initial step size for semi-major axis

    // Get magnitudes of vectors and angle between them via dot product
    double r1_mag = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    double r2_mag = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);

    // Compute r, h, and theta unit vectors for both positions on the trajectory
    double r1_unit[] = {r1[0]/r1_mag, r1[1]/r1_mag, r1[2]/r1_mag};          // unit vector, nondimensional
    double r2_unit[] = {r2[0]/r2_mag, r2[1]/r2_mag, r2[2]/r2_mag};          // unit vector, nondimensional
    
    // Do cross product, then divide by magnitude
    double h_unit[] = {	r1[1]*r2[2] - r1[2]*r2[1],      // h = r1 x r2 (direction only!), nondimensional
                        r1[2]*r2[0] - r1[0]*r2[2],
                        r1[0]*r2[1] - r1[1]*r2[0]};
    double mag = sqrt(h_unit[0]*h_unit[0] + h_unit[1]*h_unit[1] + h_unit[2]*h_unit[2]);
    h_unit[0] /= mag;	h_unit[1] /= mag;	h_unit[2] /= mag;

    double theta1_unit[] = {h_unit[1]*r1[2] - h_unit[2]*r1[1],     // unit vector, theta1 = h x r1, nondimensional
                            h_unit[2]*r1[0] - h_unit[0]*r1[2],
                            h_unit[0]*r1[1] - h_unit[1]*r1[0]};
    mag = sqrt(theta1_unit[0]*theta1_unit[0] + theta1_unit[1]*theta1_unit[1] + theta1_unit[2]*theta1_unit[2]);
    theta1_unit[0] /= mag; 	theta1_unit[1] /= mag; 	theta1_unit[2] /= mag;

    double theta2_unit[] = {h_unit[1]*r2[2] - h_unit[2]*r2[1],     // unit vector, theta2 = h x r2, nondimensional
                            h_unit[2]*r2[0] - h_unit[0]*r2[2],
                            h_unit[0]*r2[1] - h_unit[1]*r2[0]};
    mag = sqrt(theta2_unit[0]*theta2_unit[0] + theta2_unit[1]*theta2_unit[1] + theta2_unit[2]*theta2_unit[2]);
    theta2_unit[0] /= mag; 	theta2_unit[1] /= mag; 	theta2_unit[2] /= mag;

    // Determine the orientation of the orbit plane and the angular distance of both points from the ascending node
    double inc = acos(boundValue(h_unit[2], -1.0, 1.0));		// inclination, rad, choose positive
    double raan = resolveAngle(asin(boundValue(h_unit[0]/sin(inc), -1.0, 1.0)), 
    	acos(boundValue(-h_unit[1]/sin(inc), -1.0, 1.0)));
    double theta1 = resolveAngle(asin(boundValue(r1_unit[2]/sin(inc), -1.0, 1.0)),
    	acos(boundValue(theta1_unit[2]/sin(inc), -1.0, 1.0)));
    double theta2 = resolveAngle(asin(boundValue(r2_unit[2]/sin(inc), -1.0, 1.0)),
    	acos(boundValue(theta2_unit[2]/sin(inc), -1.0, 1.0)));

    node1.setExtraParam("inc", inc);
    node1.setExtraParam("raan", raan);

    node2.setExtraParam("inc", inc);
    node2.setExtraParam("raan", raan);

    double eta = acos((r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2])/(r1_mag*r2_mag)); // double angle possible; consider one value only for now
    double transAngle = eta;

    // Type 2 transfers: use the opposite side of the transfer arc
    if(type == 2)
        transAngle = 2*PI - eta;

    // "Space Triangle" Geometry
    double c = sqrt(r1_mag*r1_mag + r2_mag*r2_mag - 2*r1_mag*r2_mag*cos(eta));
    double s = 0.5*(r1_mag + r2_mag + c);

    // Compute the parabolic TOF (limiting case between elliptic and hyperbolic)
    double tof_par = type == 1 ? 
        sqrt(2.0/pSys->getMu())*(pow(s, 1.5) - pow(s - c, 1.5))/3.0 :     // Type 1
        sqrt(2.0/pSys->getMu())*(pow(s, 1.5) + pow(s - c, 1.5))/3.0;      // Type 2

    if(verbosity >= Verbosity_tp::ALL_MSG){
	    printf("Space triangle:\n  r1 = %.4f km\n  r2 = %.4f km\n  c = %.4f km\n  s = %.4f km\n", r1_mag, r2_mag, c, s);
	    printf("  Transfer Angle = %.2f deg\n", transAngle*180/PI);
	    printf("\nParabolic Transfer:\n  tof_par = %.2f s\n", tof_par);
	}

	double a = NAN, e = NAN, p = NAN, tof_guess = 0, alpha = 0, beta = 0;
	unsigned int subType = 0;
    if(tof > tof_par){  // Elliptical transfer
        double a_min = 0.5*s;
        double alpha0 = 2*asin(sqrt(s/(2*a_min)));
        double beta0 = 2*asin(sqrt((s - c)/(2*a_min)));

        double tof_min = type == 1 ?
            sqrt(pow(a_min,3)/pSys->getMu())*(alpha0 - sin(alpha0) - beta0 + sin(beta0)) :  // Type 1
            sqrt(pow(a_min,3)/pSys->getMu())*(alpha0 - sin(alpha0) + beta0 - sin(beta0));   // Type 2

        subType = tof < tof_min ? 1 : 2;    // 1 = A, 2 = B
        double a_guess = tof < tof_min ? a_min : a_min*pow(tof/tof_min, 0.25);

        if(verbosity >= Verbosity_tp::ALL_MSG){
        	printf("Transfer Type: %u%c\n", type, subType == 1 ? 'A' : 'B');
        	printf("Initial guess for a = %.4f km\n", a_guess);
        }
        // Iterate to find semi-major axis for the transfer arc
        unsigned int count = 0;
        double last_tof_guess = 0, last_a_guess = 0;
        while(std::abs(tof - tof_guess) > tof_maxErr && count < maxIts){
            last_a_guess = a_guess;
            last_tof_guess = tof_guess;

            a_guess += a_step;
            count++;

            // Compute the principle values of alpha and beta
            alpha = 2*asin(sqrt(s/(2*a_guess)));
            beta = 2*asin(sqrt((s-c)/(2*a_guess)));

            // Adjust alpha and beta to their correct quadrants
            if(subType == 2){
                alpha = 2*PI - alpha;
            }
            if(type == 2){
                beta *= -1;
            }

            tof_guess = sqrt(pow(a_guess,3)/pSys->getMu())*(alpha - beta - sin(alpha) + sin(beta));

            if(count == 1){
                // Adjust a_step if it is high or low and adjust step direction to go the right way
                if(subType == 2){
                    if( (tof_guess > tof && a_step > 0) || (tof_guess < tof && a_step < 0) )
                        a_step *= -1;
                }else{
                    if( (tof_guess < tof && a_step > 0) || (tof_guess > tof && a_step < 0) )
                        a_step *= -1;
                }
            }else{
                if(subType == 2 && a_step < 0){
                    // Make sure the guessing doesn't pass a_min
                    if(a_guess < a_min)
                        break;
                }

                // If the search passes the correct value, turn around and take smaller steps!
                if(sign(tof - last_tof_guess) != sign(tof - tof_guess))
                    a_step *= -0.1;
            }

            if(std::abs(tof - tof_guess) > tof_maxErr){
            	// try a big step
            	double m = (tof_guess - last_tof_guess)/(a_guess - last_a_guess);
            	a_guess = (tof - last_tof_guess)/m + last_a_guess;
            }

            printVerb(verbosity >= Verbosity_tp::DEBUG, "It %03u: a_step = %.4f km, a = %.4f km, tof_err = %.4f sec\n",
            	count, a_step, a_guess, tof - tof_guess);

        }

        printVerb(verbosity >= Verbosity_tp::ALL_MSG, "\nSteps Taken: %u\n", count);

        if(std::abs(tof - tof_guess) > tof_maxErr)
            throw DivergeException("No match for the specified TOF could be found!");

        a = a_guess;
        double p_option1 = (4*a*(s - r1_mag)*(s - r2_mag)/(c*c)) * sin(0.5*(alpha + beta))*sin(0.5*(alpha + beta));
        double p_option2 = (4*a*(s - r1_mag)*(s - r2_mag)/(c*c)) * sin(0.5*(alpha - beta))*sin(0.5*(alpha - beta));

        // 1A and 2B take the max value, 1B and 2A take the min value
        p = type == subType? std::max(p_option1, p_option2) : std::min(p_option1, p_option2);

        printVerb(verbosity >= Verbosity_tp::ALL_MSG, 
        	"\nMinimum Energy Transfer:\n  a_min = %.4f km\n  alpha0 = %.4f deg\n  beta0 = %.4f deg\n  tof_min = %.2f s\n",
        	a_min, alpha0*180/PI, beta0*180/PI, tof_min);
    }else{  // tof < tof_par: Hyperbolic transfer
    	subType = 3;	// hyperbola, "H"
    	a_step = -1*std::abs(a_step);	// Make sure step is negative for hyperbolic

    	// Iterate to find semi-major axis for the transfer arc
        unsigned int count = 0;
        double last_tof_guess = 0, a_guess = 0;
        while(std::abs(tof - tof_guess) > tof_maxErr && count < maxIts){
        	a_guess += a_step;

        	// Compute principal values of alpha and beta
        	double alpha = 2.0*asinh(sqrt(s/(-2.0*a_guess)));
        	double beta = 2.0*asinh(sqrt((s-c)/(-2.0*a_guess)));

        	if(type == 2)
        		beta *= -1;

        	tof_guess = sqrt(-pow(a,3)/pSys->getMu())*(beta - alpha + sinh(alpha) - sinh(beta));

        	if(count > 0 && sign(tof - last_tof_guess) != sign(tof - tof_guess))
        		a_step *= -0.1;

        	last_tof_guess = tof_guess;
        	count++;
        }

        if(std::abs(tof - tof_guess) > tof_maxErr)
        	throw DivergeException("LambertArcEngine::getLambertArc: No match for TOF on hyperbolic transfer");
        
        a = a_guess;
        double p_option1 = (-4*a*(s - r1_mag)*(s - r2_mag)/(c*c)) * (sinh((alpha + beta)/2.0)*sinh((alpha + beta)/2.0));
        double p_option2 = (-4*a*(s - r1_mag)*(s - r2_mag)/(c*c)) * (sinh((alpha + beta)/2.0)*sinh((alpha - beta)/2.0));

        p = type == 1 ? std::max(p_option1, p_option2) : std::min(p_option1, p_option2);
    }

    // Compute other orbit parameters
    e = sqrt(1.0 - p/a);
    double ta1 = acos( (p - r1_mag)/(r1_mag*e) );
    double ta2 = acos( (p - r2_mag)/(r2_mag*e) );

    /*  Four Possible transfer options
     *  1) ta_a > ta_d > 0
     *  2) ta_d > 0 > ta_a
     *  3) ta_a > 0 > ta_d
     *  4) 0 > ta_d > ta_a
     */
     double possibleTx[] = {ta2 - ta1,
                            2*PI - ta2 - ta1,
                            ta2 + ta1,
                            2*PI - ta2 + ta1};
    bool foundMatch = false;
    for(unsigned int i = 0; i < 4; i++){
        if(std::abs(possibleTx[i] - transAngle) < ta_maxErr){
            switch(i){
                case 1: ta2 = 2*PI - ta2; break;
                case 2: ta1 = 2*PI - ta1; break;
                case 3: 
                    ta2 = 2*PI - ta2;
                    ta1 = 2*PI - ta1;
                    break;
            }
            foundMatch = true;
            break;
        }
    }

    if(!foundMatch)
        throw Exception("r2bp_getLambertArc: Could not find a set of true anomalies to match the transfer angle");

    double v1_mag = sqrt(2*pSys->getMu()*(1.0/r1_mag - 1.0/(2*a)));
    double v2_mag = sqrt(2*pSys->getMu()*(1.0/r2_mag - 1.0/(2*a)));
    double fpa1 = acos(sqrt(pSys->getMu() * p)/(r1_mag*v1_mag));
    double fpa2 = acos(sqrt(pSys->getMu() * p)/(r2_mag*v2_mag));

    if(ta1 > PI)
        fpa1 *= -1;

    if(ta2 > PI)
        fpa2 *= -1;

    char subTypeChars[] {'A', 'B', 'H'};
    if(verbosity >= Verbosity_tp::ALL_MSG){
        printf("\nActual Transfer: Type %d%c\n", type, subTypeChars[subType]);
        printf("  inc = %.2f deg\n  raan = %.2f deg\n", inc*180/PI, raan*180/PI);
        printf("  a = %.4f km\n  tof = %.2f sec\n  alpha = %.4f deg\n  beta = %.4f deg\n", a, tof_guess, alpha*180/PI, beta*180/PI);
        printf("  p = %.4f km\n", p);
        printf("  e = %.4f\n  ta1 = %.4f deg\n  v1 = %.4f km/s\n  fpa1 = %.4f deg\n", e, ta1*180/PI, v1_mag, fpa1*180/PI);
        printf("  ta2 = %.4f deg\n  v2 = %.4f km/s\n  fpa2 = %.4f deg\n", ta2*180/PI, v2_mag, fpa2*180/PI);
        printf("  argPeri = %.2f deg\n", (theta1 - ta1)*180/PI);
	}

    node1.setExtraParam("sma", a);
    node1.setExtraParam("ecc", e);
    node1.setExtraParam("fpa", fpa1);
    node1.setExtraParam("ta", ta1);
    node1.setExtraParam("range", r1_mag);
    node1.setExtraParam("speed", v1_mag);
    node1.setExtraParam("argPeri", (theta1 - ta1));

    node2.setExtraParam("sma", a);
    node2.setExtraParam("ecc", e);
    node2.setExtraParam("fpa", fpa2);
    node2.setExtraParam("ta", ta2);
    node2.setExtraParam("range", r2_mag);
    node2.setExtraParam("speed", v2_mag);

    if(tof > tof_par){
    	// Compute Eccentric and Mean Anomaly
    	double E1 = acos((1.0 - r1_mag/a)/e);

    	// Correct for quadrant ambiguity
    	if(sign(ta1) < 0)
    		E1 = 2*PI - E1;

    	double E2 = E1 + (alpha - beta);

    	double M1 = E1 - e*sin(E1);
    	double M2 = E2 - e*sin(E2);

    	node1.setExtraParam("E", E1);
    	node1.setExtraParam("M", M1);

    	node2.setExtraParam("E", E2);
    	node2.setExtraParam("M", M2);
    }
    
    // Compute the velocity vector for each node
    double v1[] = {v1_mag*sin(fpa1), v1_mag*cos(fpa1), 0};	// in r-theta-h coordinates
    double v2[] = {v2_mag*sin(fpa2), v2_mag*cos(fpa2), 0};

    // DCM that converts r-theta-h coordinates to x-y-z coordinates
	// The first row is (x dot r, x dot theta), second row is (y dot r, y dot theta),
	// third row is (z dot r, z dot theta)
    double dcm1[][2] =	{{cos(raan)*cos(theta1) - sin(raan)*cos(inc)*sin(theta1), -cos(raan)*sin(theta1) - sin(raan)*cos(inc)*cos(theta1)},
                     	 {sin(raan)*cos(theta1) + cos(raan)*cos(inc)*sin(theta1), -sin(raan)*sin(theta1) + cos(raan)*cos(inc)*cos(theta1)},
                     	 {sin(inc)*sin(theta1), sin(inc)*cos(theta1)}
                    	};

    double dcm2[][2] =	{{cos(raan)*cos(theta2) - sin(raan)*cos(inc)*sin(theta2), -cos(raan)*sin(theta2) - sin(raan)*cos(inc)*cos(theta2)},
                     	 {sin(raan)*cos(theta2) + cos(raan)*cos(inc)*sin(theta2), -sin(raan)*sin(theta2) + cos(raan)*cos(inc)*cos(theta2)},
                     	 {sin(inc)*sin(theta2), sin(inc)*cos(theta2)}
                    	};

    state1[3] = v1[0]*dcm1[0][0] + v1[1]*dcm1[0][1];
	state1[4] = v1[0]*dcm1[1][0] + v1[1]*dcm1[1][1];
	state1[5] = v1[0]*dcm1[2][0] + v1[1]*dcm1[2][1];

	state2[3] = v2[0]*dcm2[0][0] + v2[1]*dcm2[0][1];
	state2[4] = v2[0]*dcm2[1][0] + v2[1]*dcm2[1][1];
	state2[5] = v2[0]*dcm2[2][0] + v2[1]*dcm2[2][1];

	node1.setState(state1);
	node2.setState(state2);

    int id1 = traj.addNode(node1);
    int id2 = traj.addNode(node2);
    traj.addSeg(Segment(id1, id2, tof));
    return traj;
}//====================================================


//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void LambertArcEngine::copyMe(const LambertArcEngine &e){
	Engine::copyBaseEngine(e);
	tof_maxErr = e.tof_maxErr;
	ta_maxErr = e.ta_maxErr;
	maxIts = e.maxIts;
}//====================================================

void LambertArcEngine::cleanEngine(){
	bIsClean = true;
}//====================================================

void LambertArcEngine::reset(){
	if(!bIsClean)
		cleanEngine();

	tof_maxErr = 1e-3;
	ta_maxErr = 0.01*PI/180;
	maxIts = 100;
}//====================================================


}// End of Astrohelion namespace