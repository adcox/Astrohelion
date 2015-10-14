/**
 *	@file tpat_linear_motion_engine.cpp
 *	@brief Uses linear EOMS near libration points to generate trajectories
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

#include "tpat_linear_motion_engine.hpp"

#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_traj.hpp"
#include "tpat_traj_step.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <complex>
#include <vector>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Default, do-nothing constructor
 */
tpat_linear_motion_engine::tpat_linear_motion_engine(){}

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the ratio between short- and long-period motion near L4 and L5
 *
 *	This ratio is only applied to mixed-period motion simulations
 *	@return the ratio between short- and long-period motion near L4 and L5
 */
double tpat_linear_motion_engine::getMPORatio() const { return nu; }

/**
 *	@brief Retrieve the number of revolutions to simulate for
 *
 *	Note that this rev count applies to the in-plane oscillations only; the out-of-plane
 *	motion will most likely have a different period and therefore perform a different
 *	number of revs during the same time period
 *
 *	@return the number of rotatins to simulate for
 */
int tpat_linear_motion_engine::getNumRevs() const { return rots; }

/**
 *	@brief Retrieve the step size (in non-dimensional units) for the time vector
 *	@return he step size (in non-dimensional units) for the time vector
 */
double tpat_linear_motion_engine::getTimeStep() const { return t_step; }

/**
 *	@brief Retrieve the acceptable numerical tolerance for computations
 *
 *	This tolerance is used to target Lagrange point locations and to determine
 *	where or not numbers are "equal" to zero (or any other value)
 *	@return the acceptable numerical tolerance for computations
 */
double tpat_linear_motion_engine::getTol() const { return tol; }

/**
 *	@brief Set the ratio for short- and long-period L4 and L5 motion
 *
 *	This ratio is only applied to mixed-period simulations
 *	@param ratio the ratio
 */
void tpat_linear_motion_engine::setMPORatio(double ratio) {
	nu = ratio;
}

/**
 *	@brief Set the number of revolutions to simulate for.
 *
 *	The length of one revolutions is determined by the in-plane period, not the
 *	out-of-plane period
 *	@param numRevs number of revolutions
 */
void tpat_linear_motion_engine::setNumRevs(int numRevs) { rots = numRevs; }

/**
 *	@brief Set the step size for the time vector
 *	@param dt non-dimensional time step
 */
void tpat_linear_motion_engine::setTimeStep(double dt) { t_step = dt; }

/**
 *	@brief set the tolerance to use
 *	@param t tolerance, non-dimensional units
 */
void tpat_linear_motion_engine::setTol(double t){ tol = t; }

/**
 *	@brief get a human-readable string for a motion type
 *	@param type the motion type
 *	@return a human-redable string
 */
const char* tpat_linear_motion_engine::getTypeStr(motion_t type) const{
	switch(type){
		case NONE: return "NONE";
		case HYP: return "HYPERBOLIC";
		case ELLIP: return "ELLIPTIC";
		case SPO: return "SHORT-period_xy-ORBIT";
		case LPO: return "LONG-period_xy-ORBIT";
		case MPO: return "MIXED-period_xy-ORBIT";
		case CONVERGE: return "CONVERGENT";
		case DIVERGE: return "DIVERGENT";
		default: return "Unrecognized type...";
	}
}//========================================

//-----------------------------------------------------
//      Linear Motion Generation Functions
//-----------------------------------------------------

/**
 * 	@brief Construct a trajectory from linear approximations of the CR3BP EOMs
 *
 *	This function chooses the default motion types: Elliptical for the collinear points and
 *	SPO for the triangular points. If the system has a mass ratio such that it falls into 
 *	Case III, then convergent behavior will be employed.
 *
 *	@param L the Lagrange point number [1,5]
 *	@param r0 a three-element vector that specifies the initial position of the arc relative
 *		to the chosen Lagrange point
 *	@param sys the CR3BP system data object
 */
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLinear(int L, double r0[3], tpat_sys_data_cr3bp *sys){

	return getCR3BPLinear(L, r0, 0, 0, NONE, sys);
}//=======================================================

/**
 *	@brief Compute a linear approximation for a Lissajous orbit
 *
 *	@param L Lagrange point number. Choose 1, 2, or 3
 *	@param Axy in-plane amplitude, non-dimensional units
 *	@param xAmp whether or not Axy describes the x-amplitude; true -> Ax = Axy, false -> Ay = Axy
 *	@param phi starting phase angle for the in-plane motion
 *	@param Az out-of-plane amplitude, non-dimensional units
 *	@param psi starting phase angle for out-of-plane motion
 *	@param sysData a pointer to a system data object
 */
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLiss(int L, double Axy, bool xAmp, double phi, double Az, double psi,
	tpat_sys_data_cr3bp *sysData){

	double mu = sysData->getMu();

	// Locate Lagrange point
	double LPtPos[3] = {0};
	cr3bp_getEquilibPt(*sysData, L, tol, LPtPos);

	// Get partial derivatives of pseudo potential
	double ddots[6] = {0};
	cr3bp_getUDDots(mu, LPtPos[0], LPtPos[1], LPtPos[2], ddots);

	tpat_traj_cr3bp linTraj(sysData);

	double xi;
	double eta;
	double zeta;
	double xi_dot;
	double eta_dot;
	double zeta_dot;

	// Out of plane frequency
	double w_z = sqrt(-ddots[2]);

	if(L < 4){
		// Compute eigenvalues analytically
		double beta1 = 2 - (ddots[0] + ddots[1])/2;
		double beta2 = sqrt(-1*ddots[1]*ddots[0]);
		std::complex<double> Lam1 = -1*beta1 + sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> Lam2 = -1*beta1 - sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> eigenval[] = {sqrt(Lam1), -1.0*sqrt(Lam1), sqrt(Lam2), -1.0*sqrt(Lam2)};

		// Motion for ELLIP type only
		double s = std::imag(eigenval[2]);
		double beta3 = (s*s + ddots[0])/(2*s);
		double period_xy = 2*PI/s;

		// Initial conditions computed from amplitude/phase angle information
		double xi0 = xAmp ? Axy*cos(phi) : -Axy/beta3 * cos(phi);
		double eta0 = xAmp ? -Axy*beta3*sin(phi) : Axy*sin(phi);

		for(double t = 0; t < rots*period_xy; t += t_step){
			xi = xi0*cos(s*t) + eta0/beta3*sin(s*t);
			eta = eta0*cos(s*t) - beta3*xi0*sin(s*t);
			zeta = Az*sin(w_z*t + psi);
			xi_dot = -s*xi0*sin(s*t) + eta0*s/beta3 * cos(s*t);
			eta_dot = -s*eta0*sin(s*t) - s*beta3*xi0*cos(s*t);
			zeta_dot = -w_z*Az*cos(w_z*t + psi);

			double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

			tpat_traj_step step(state, t);
			linTraj.appendStep(step);
		}
	}else{
		throw tpat_exception("tpat_linear_motion_engine::getCR3BPLiss: Cannot compute Lissajous motion for anything other than the collinear points");
	}

	// Compute Jacobi Constant for each step; won't be constant because non-linear dynamics are
	// not enforced, but is still useful information
	for(int i = 0; i < linTraj.getLength(); i++){
		std::vector<double> state = linTraj.getState(i);
		linTraj.setJacobi(i, cr3bp_getJacobi(&(state[0]), mu));
	}
	
	return linTraj;
}//========================================================


/**
 * 	@brief Construct a trajectory from linear approximations of the CR3BP EOMs
 *
 *	This function sets the out-of-plane motion to have zero amplitude
 *
 *	@param L the Lagrange point number [1,5]
 *	@param r0 a three-element vector that specifies the initial position of the arc relative
 *		to the chosen Lagrange point
 *	@param type the type of linearized motion desired
 *	@param sysData the CR3BP system data object
 *
 *	@return a trajectory object containing one revolution of the trajectory. Because this motion
 *	is generated from simplified dynamics, no information about the STM or Jacobi Constant is 
 *	computed. Accelerations are also not computed. These values are all stored as NAN
 */
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLinear(int L, double r0[3], motion_t type, tpat_sys_data_cr3bp *sysData){

	return getCR3BPLinear(L, r0, 0, 0, type, sysData);
}//====================================================

/**
 * 	@brief Construct a trajectory from linear approximations of the CR3BP EOMs
 *
 *	@param L the Lagrange point number [1,5]
 *	@param r0 a three-element vector that specifies the initial position of the arc relative
 *		to the chosen Lagrange point
 *	@param Az out-of-plane amplitude, in non-dimensional units
 *	@param psi defines the starting elevation angle of the out-of-plane motion, radians
 *	@param type the type of linearized motion desired
 *	@param sysData the CR3BP system data object
 *
 *	@return a trajectory object containing one revolution of the trajectory. Because this motion
 *	is generated from simplified dynamics, no information about the STM or Jacobi Constant is 
 *	computed. Accelerations are also not computed. These values are all stored as NAN
 */
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLinear(int L, double r0[3], double Az, double psi,
	motion_t type, tpat_sys_data_cr3bp *sysData){

	double mu = sysData->getMu();

	// Locate Lagrange point
	double LPtPos[3] = {0};
	cr3bp_getEquilibPt(*sysData, L, tol, LPtPos);

	// Get partial derivatives of pseudo potential
	double ddots[6] = {0};
	cr3bp_getUDDots(mu, LPtPos[0], LPtPos[1], LPtPos[2], ddots);

	double xi0 = r0[0];		// Initial x-variation
	double eta0 = r0[1];	// Initial y-variation

	tpat_traj_cr3bp linTraj(sysData);

	double xi;
	double eta;
	double zeta;
	double xi_dot;
	double eta_dot;
	double zeta_dot;
	double period_xy;

	// Out of plane frequency
	double w_z = sqrt(-ddots[2]);

	if(L < 4){
		// Compute eigenvalues analytically
		double beta1 = 2 - (ddots[0] + ddots[1])/2;
		double beta2 = sqrt(-1*ddots[1]*ddots[0]);
		std::complex<double> Lam1 = -1*beta1 + sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> Lam2 = -1*beta1 - sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> eigenval[] = {sqrt(Lam1), -1.0*sqrt(Lam1), sqrt(Lam2), -1.0*sqrt(Lam2)};

		double s;
		switch(type){
			case NONE:	// for default behavior
			case ELLIP:
			{
				s = std::imag(eigenval[2]);
				double beta3 = (s*s + ddots[0])/(2*s);
				period_xy = 2*PI/s;

				for(double t = 0; t < rots*period_xy; t += t_step){
					xi = xi0*cos(s*t) + eta0/beta3*sin(s*t);
					eta = eta0*cos(s*t) - beta3*xi0*sin(s*t);
					zeta = Az*sin(w_z*t + psi);
					xi_dot = -s*xi0*sin(s*t) + eta0*s/beta3 * cos(s*t);
					eta_dot = -s*eta0*sin(s*t) - s*beta3*xi0*cos(s*t);
					zeta_dot = -w_z*Az*cos(w_z*t + psi);

					double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

					tpat_traj_step step(state, t);
					linTraj.appendStep(step);
				}
				break;
			}
			case HYP:
			{
				s = std::real(eigenval[0]);
				double alpha = (s*s - ddots[0])/(2*s);
				period_xy = 2*PI/s;

				for(double t = 0; t < rots*period_xy; t += t_step){
					xi = xi0*cosh(s*t) + eta0/alpha*sinh(s*t);
					eta = eta0*cosh(s*t) + alpha*xi0*sinh(s*t);
					
					xi_dot = s*xi0*sinh(s*t) + s*eta0/alpha*cosh(s*t);
					eta_dot = s*eta0*sinh(s*t) + s*alpha*xi0*cosh(s*t);
					
					zeta = Az*sin(w_z*t + psi);
					zeta_dot = -w_z*Az*cos(w_z*t + psi);

					double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

					tpat_traj_step step(state, t);
					linTraj.appendStep(step);
				}
	            break;
			}
			default:
				throw tpat_exception("Invalid type for collinear points");
		}
	}else{ // L = 4 or 5
		// Compute eigenvalues analytically
		std::complex<double> gc = 1 - 27*mu*(1-mu);
		std::complex<double> Lam1 = -0.5 + 0.5*sqrt(gc);
		std::complex<double> Lam2 = -0.5 - 0.5*sqrt(gc);
		std::complex<double> eigenval[] = {sqrt(Lam1), -1.0*sqrt(Lam1), -1.0*sqrt(Lam2), sqrt(Lam2)};

		// Other variables that depend on eigenvalues
		std::complex<double> s1 = eigenval[0];
		std::complex<double> s3 = eigenval[2];
		std::complex<double> alpha1 = (s1*s1 - ddots[0])/(2.0*s1 + ddots[3]);
		std::complex<double> alpha3 = (s3*s3 - ddots[0])/(2.0*s3 + ddots[3]);
		double a1 = std::real(alpha1);
		double a3 = std::real(alpha3);
		double b1 = std::imag(alpha1);
		double b3 = std::imag(alpha3);

		

		// Will be a real number, so re-cast for easier use
		double g = std::real(gc);
		if(g > tol){	// Case I, g > 0 (allow for some tiny error)
			// Make these scalars now that we're done computing
			double s1d = std::imag(s1);
			double s3d = std::imag(s3);
			
			if(s1d < 0 || s3d < 0){
				printErr("tpat_linear_motion_engine :: Eigenvalue order appears to be shifting... please hardcode abs()\n");
			}

			switch(type){
				case LPO:
					period_xy = 2*PI/s1d;

					for(double t = 0; t < rots*period_xy; t+= t_step){
						xi = xi0*cos(s1d*t) + (eta0 - a1*xi0)/b1 * sin(s1d*t);
						eta = eta0*cos(s1d*t) - (b1*xi0 - a1*(eta0 - a1*xi0)/b1)*sin(s1d*t);
						xi_dot = -s1d*xi0*sin(s1d*t) + s1d*(eta0 - a1*xi0)/b1 * cos(s1d*t);
						eta_dot = -s1d*eta0*sin(s1d*t) - s1d*(b1*xi0 - a1*(eta0 - a1*xi0)/b1)*cos(s1d*t);

						zeta = Az*sin(w_z*t + psi);
						zeta_dot = -w_z*Az*cos(w_z*t + psi);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
                	break;
                case NONE: // for default behavior
				case SPO:
					period_xy = 2*PI/s3d;

					for(double t = 0; t < rots*period_xy; t+= t_step){
						xi = xi0*cos(s3d*t) + (eta0 - a3*xi0)/b3 * sin(s3d*t);
						eta = eta0*cos(s3d*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(s3d*t);
						xi_dot = -s3d*xi0*sin(s3d*t) + s3d*(eta0 - a3*xi0)/b3 * cos(s3d*t);
						eta_dot = -s3d*eta0*sin(s3d*t) - s3d*(b3*xi0 - a3*(eta0 - a3*xi0)/b3)*cos(s3d*t);

						zeta = Az*sin(w_z*t + psi);
						zeta_dot = -w_z*Az*cos(w_z*t + psi);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
                	break;
				case MPO:
				{
					period_xy = s1d < s3d ? 2*PI/s1d : 2*PI/s3d;
					double R1 = xi0/(1 + nu);
                	double C1 = xi0*(a1 + nu*a3)/(2*(1 + nu)*(b1 + nu*b3)) - eta0/(2*(b1 + nu*b3));

                	for(double t = 0; t < rots*period_xy; t+= t_step){
						xi = R1*cos(s1d*t) - 2*C1*sin(s1d*t) + nu*R1*cos(s3d*t) - 2*nu*C1*sin(s3d*t);
						eta = (a1*R1 - 2*b1*C1)*cos(s1d*t) - (b1*R1 + 2*a1*C1)*sin(s1d*t) +
                    		nu*(a3*R1 - 2*b3*C1)*cos(s3d*t) - nu*(b3*R1 + 2*a3*C1)*sin(s3d*t);

                    	xi_dot = -s1d*R1*sin(s1d*t) - s1d*2*C1*cos(s1d*t) - s3d*nu*R1*sin(s3d*t) - s3d*2*nu*C1*cos(s3d*t);
						eta_dot = -s1d*(a1*R1 - 2*b1*C1)*sin(s1d*t) - s1d*(b1*R1 + 2*a1*C1)*cos(s1d*t) +
                    		-s3d*nu*(a3*R1 - 2*b3*C1)*sin(s3d*t) - s3d*nu*(b3*R1 + 2*a3*C1)*cos(s3d*t);
						
						zeta = Az*sin(w_z*t + psi);
						zeta_dot = -w_z*Az*cos(w_z*t + psi);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
                	break;
                }
				default:
					throw tpat_exception("Invalid type for triangular points, Case I");
			}
		}else if(std::abs(g) < tol){
			// Case II, g is 0 (approx.)
			double s1d = s1.imag();
			period_xy = 2*PI/s1d;

			for(double t = 0; t < rots*period_xy; t+= t_step){
				xi = xi0*cos(s1d*t) - (a1*xi0 - eta0)/b1 * sin(s1d*t);
				eta = eta0*cos(s1d*t) - (xi0*std::abs(alpha1) - a1*eta0)/b1 * sin(s1d*t);

				xi_dot = -s1d*xi0*sin(s1d*t) - s1d*(a1*xi0 - eta0)/b1 * cos(s1d*t);
				eta_dot = -s1d*eta0*sin(s1d*t) - s1d*(xi0*std::abs(alpha1) - a1*eta0)/b1 * cos(s1d*t);

				zeta = Az*sin(w_z*t + psi);
				zeta_dot = -w_z*Az*cos(w_z*t + psi);

				double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

				tpat_traj_step step(state, t);
				linTraj.appendStep(step);
			}
		}else{
			// Case III
			double p = std::real(s1);
			double q = std::imag(s1);
			switch(type){
				case NONE: // for default behavior
				case CONVERGE:
					period_xy = 2*PI/q;

					for(double t = 0; t < rots*period_xy; t += t_step){
						xi = exp(-p*t)*(xi0*cos(q*t) + (eta0 - a3*xi0)/b3 * sin(q*t));
	                    eta = exp(-p*t)*(eta0*cos(q*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(q*t));
	                    
	                    xi_dot = -p*exp(-p*t)*(xi0*cos(q*t) + (eta0 - a3*xi0)/b3 * sin(q*t)) + 
	                    	exp(-p*t)*(-q*xi0*sin(q*t) + q*(eta0 - a3*xi0)/b3 * cos(q*t));
	                    eta_dot = -p*exp(-p*t)*(eta0*cos(q*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(q*t)) +
	                    	exp(-p*t)*(-q*eta0*sin(q*t) - q*(b3*xi0 - a3*(eta0 - a3*xi0)/b3)*cos(q*t));

	                    zeta = Az*sin(w_z*t + psi);
						zeta_dot = -w_z*Az*cos(w_z*t + psi);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], zeta + LPtPos[2], xi_dot, eta_dot, zeta_dot};

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
					break;	
				case DIVERGE:
					throw tpat_exception("Triangular points, Case III: Diverge not yet implemented");
					break;
				default:
					throw tpat_exception("Invalid type for triangular points, Case III");
			}
		}
	}

	// Compute Jacobi Constant for each step; won't be constant because non-linear dynamics are
	// not enforced, but is still useful information
	for(int i = 0; i < linTraj.getLength(); i++){
		std::vector<double> state = linTraj.getState(i);
		linTraj.setJacobi(i, cr3bp_getJacobi(&(state[0]), mu));
	}

	return linTraj;
}//====================================================

