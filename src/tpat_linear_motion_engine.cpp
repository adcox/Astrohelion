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
		case SPO: return "SHORT-PERIOD-ORBIT";
		case LPO: return "LONG-PERIOD-ORBIT";
		case MPO: return "MIXED-PERIOD-ORBIT";
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
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLinear(int L, double r0[3], tpat_sys_data_cr3bp sys){

	return getCR3BPLinear(L, r0, NONE, sys);
}

/**
 * 	@brief Construct a trajectory from linear approximations of the CR3BP EOMs
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
tpat_traj_cr3bp tpat_linear_motion_engine::getCR3BPLinear(int L, double r0[3], motion_t type, tpat_sys_data_cr3bp sysData){

	double mu = sysData.getMu();

	// Locate Lagrange point
	double LPtPos[3] = {0};
	cr3bp_getEquilibPt(sysData, L, tol, LPtPos);

	// Get partial derivatives of pseudo potential
	double ddots[6] = {0};
	cr3bp_getUDDots(sysData.getMu(), LPtPos[0], LPtPos[1], LPtPos[2], ddots);

	double xi0 = r0[0];		// Initial x-variation
	double eta0 = r0[1];	// Initial y-variation

	tpat_traj_cr3bp linTraj(&sysData);

	double xi;
	double eta;
	double xi_dot0;
	double eta_dot0;
	double period;

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
				s = imag(eigenval[2]);
				double beta3 = (s*s + ddots[0])/(2*s);
				period = 2*PI/s;
				xi_dot0 = s*eta0/beta3;
				eta_dot0 = -1*s*beta3*xi0;

				for(double t = 0; t < rots*period; t += t_step){
					xi = xi0*cos(s*t) + eta0/beta3*sin(s*t);
					eta = eta0*cos(s*t) - beta3*xi0*sin(s*t);

					double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

					if(t == 0){
						state[3] = xi_dot0;
						state[4] = eta_dot0;
					}

					tpat_traj_step step(state, t);
					linTraj.appendStep(step);
				}
				break;
			}
			case HYP:
			{
				s = std::real(eigenval[0]);
				double alpha = (s*s - ddots[0])/(2*s);
				period = 2*PI/s;
				xi_dot0 = s*eta0/alpha;
	            eta_dot0 = s*alpha*xi0;

				for(double t = 0; t < rots*period; t += t_step){
					xi = xi0*cosh(s*t) + eta0/alpha*sinh(s*t);
					eta = eta0*cosh(s*t) + alpha*xi0*sinh(s*t);
					
					double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

					if(t == 0){
						state[3] = xi_dot0;
						state[4] = eta_dot0;
					}

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
					period = 2*PI/s1d;
					xi_dot0 = s1d*(eta0 - a1*xi0)/b1;
                	eta_dot0 = -s1d*(b1*xi0 - a1*(eta0 - a1*xi0)/b1);
					for(double t = 0; t < rots*period; t+= t_step){
						xi = xi0*cos(s1d*t) + (eta0 - a1*xi0)/b1 * sin(s1d*t);
						eta = eta0*cos(s1d*t) - (b1*xi0 - a1*(eta0 - a1*xi0)/b1)*sin(s1d*t);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

						if(t == 0){
							state[3] = xi_dot0;
							state[4] = eta_dot0;
						}

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
                	break;
                case NONE: // for default behavior
				case SPO:
					period = 2*PI/s3d;
					xi_dot0 = s3d*(eta0 - a3*xi0)/b3;
                	eta_dot0 = -s3d*(b3*xi0 - a3*(eta0 - a3*xi0)/b3);
					for(double t = 0; t < rots*period; t+= t_step){
						xi = xi0*cos(s3d*t) + (eta0 - a3*xi0)/b3 * sin(s3d*t);
						eta = eta0*cos(s3d*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(s3d*t);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

						if(t == 0){
							state[3] = xi_dot0;
							state[4] = eta_dot0;
						}

						tpat_traj_step step(state, t);
						linTraj.appendStep(step);
					}
                	break;
				case MPO:
				{
					period = s1d < s3d ? 2*PI/s1d : 2*PI/s3d;
					double R1 = xi0/(1 + nu);
                	double C1 = xi0*(a1 + nu*a3)/(2*(1 + nu)*(b1 + nu*b3)) - eta0/(2*(b1 + nu*b3));
                	xi_dot0 = -2*C1*(s1d + nu*s3d);
                	eta_dot0 = -R1*(s1d*b1 + nu*s3d*b3) - 2*C1*(s1d*a1 + nu*s3d*a3);

                	for(double t = 0; t < rots*period; t+= t_step){
						xi = R1*cos(s1d*t) - 2*C1*sin(s1d*t) + nu*R1*cos(s3d*t) - 2*nu*C1*sin(s3d*t);
						eta = (a1*R1 - 2*b1*C1)*cos(s1d*t) - (b1*R1 + 2*a1*C1)*sin(s1d*t) +
                    		nu*(a3*R1 - 2*b3*C1)*cos(s3d*t) - nu*(b3*R1 + 2*a3*C1)*sin(s3d*t);

						double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

						if(t == 0){
							state[3] = xi_dot0;
							state[4] = eta_dot0;
						}

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
			period = 2*PI/s1d;
			xi_dot0 = -s1d*(a1*xi0 - eta0)/b1;
            eta_dot0 = -s1d*(xi0*std::abs(alpha1) - a1*eta0)/b1;

			for(double t = 0; t < rots*period; t+= t_step){
				xi = xi0*cos(s1d*t) - (a1*xi0 - eta0)/b1 * sin(s1d*t);
				eta = eta0*cos(s1d*t) - (xi0*std::abs(alpha1) - a1*eta0)/b1 * sin(s1d*t);

				double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};

				if(t == 0){
					state[3] = xi_dot0;
					state[4] = eta_dot0;
				}

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
					period = 2*PI/q;
					xi_dot0 = -xi0*p + q*(eta0 - a3*xi0)/b3;
					eta_dot0 = -p*eta0 - q*(b3*xi0 - a3*(eta0 - a3*xi0)/b3);

					for(double t = 0; t < rots*period; t += t_step){
						xi = exp(-p*t)*(xi0*cos(q*t) + (eta0 - a3*xi0)/b3 * sin(q*t));
	                    eta = exp(-p*t)*(eta0*cos(q*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(q*t));
	                    
	                    double state[] = {xi + LPtPos[0], eta + LPtPos[1], 0, 0, 0, 0};
						
						if(t == 0){
							state[3] = xi_dot0;
							state[4] = eta_dot0;
						}

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

	return linTraj;
}//====================================================

