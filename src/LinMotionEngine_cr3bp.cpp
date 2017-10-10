/**
 *	\file LinMotionEngine_cr3bp.cpp
 *	\brief Uses linear EOMS near libration points to generate trajectories
 *	
 *	\author Andrew Cox
 *	\version September 28, 2017
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

#include "LinMotionEngine_cr3bp.hpp"

#include "Arcset_cr3bp.hpp"
#include "DynamicsModel.hpp"
#include "Exceptions.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  \brief Default, do-nothing constructor
 */
LinMotionEngine_cr3bp::LinMotionEngine_cr3bp() : LinMotionEngine(){}

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *	\brief Retrieve the ratio between short- and long-period motion near L4 and L5
 *
 *	This ratio is only applied to mixed-period motion simulations
 *	\return the ratio between short- and long-period motion near L4 and L5
 */
double LinMotionEngine_cr3bp::getMPORatio() const { return nu; }

/**
 *	\brief Set the ratio for short- and long-period L4 and L5 motion
 *
 *	This ratio is only applied to mixed-period simulations
 *	\param ratio the ratio
 */
void LinMotionEngine_cr3bp::setMPORatio(double ratio) { nu = ratio; }

/**
 *	\brief get a human-readable string for a motion type
 *	\param type the motion type
 *	\return a human-redable string
 */
const char* LinMotionEngine_cr3bp::getTypeStr(unsigned int type) const{
	switch(type){
		case LinMotion_cr3bp_tp::SPO: return "SHORT-PERIOD-XY-ORBIT";
		case LinMotion_cr3bp_tp::LPO: return "LONG-PERIOD-XY-ORBIT";
		case LinMotion_cr3bp_tp::MPO: return "MIXED-PERIOD-XY-ORBIT";
		default: return LinMotionEngine::getTypeStr(type);
	}
}//========================================


//-----------------------------------------------------------------------------
//      Linear Motion Generation Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Save data from the linear model to an Arcset object
 *  \details [long description]
 * 
 *  \param model Describes the relevant dynamical model
 *  \param pArc Pointer to the storage arcset
 *  \param t current time on the solution
 *  \param pState current state on the solution
 *  \param nodeCount current number of nodes that have been created
 *  \param eomParams object that stores parameters used to evaluate the nonlinear equations of motion
 *  \param dt time step between nodes
 *  \param t_step step between points on the solution
 *  \param core_dim dimension of the core state
 *  \param full_dim dimension of the full state (core + control + stm + extras)
 */
void LinMotionEngine_cr3bp::storeData(const DynamicsModel *model, Arcset_cr3bp *pArc, const double &t, std::vector<double> *pState,
	unsigned int &nodeCount, EOM_ParamStruct *eomParams, const double& dt, const double &t_step,
	const unsigned int &core_dim, const unsigned int &full_dim){

	if(t >= dt*nodeCount){
		// Add a node
		Node node(&(pState->front()), core_dim, t);
		int ID = model->sim_addNode(node, &(pState->front()), t, pArc, eomParams, Event_tp::NONE);
		// Link the previous segment, if it exists
		if(nodeCount > 0){
			Segment &lastSeg = pArc->getSegRefByIx(-1);
			lastSeg.setTerminus(ID);

			lastSeg.appendState(&(pState->front()), full_dim);
			lastSeg.appendTime(t);
			lastSeg.updateTOF();
		}

		// Create a new segment
		Segment seg(ID, Linkable::INVALID_ID, t_step);
		seg.appendState(&(pState->front()), full_dim);
		seg.setStateWidth(full_dim);
		seg.appendTime(t);
		model->sim_addSeg(seg, &(pState->front()), t, pArc, eomParams);
		nodeCount++;
	}else{
		Segment &lastSeg = pArc->getSegRefByIx(-1);
		lastSeg.appendState(&(pState->front()), full_dim);
		lastSeg.appendTime(t);
	}
}//====================================================

/**
 *	\brief Compute a linear approximation for a Lissajous orbit
 *
 *	\param L Lagrange point number. Choose 1, 2, or 3
 *	\param Axy in-plane amplitude, non-dimensional units
 *	\param xAmp whether or not Axy describes the x-amplitude; true -> Ax = Axy, false -> Ay = Axy
 *	\param phi starting phase angle for the in-plane motion
 *	\param Az out-of-plane amplitude, non-dimensional units
 *	\param psi starting phase angle for out-of-plane motion
 *	\param pArc pointer to an arcset to store the linear trajectory in
 *	\param numNodes the number of nodes to place around the linear trajectory; 
 *	this is nodes for the entire propagation (period * num revs)
 *	
 *	\throws Exception if <tt>L</tt> is not 1, 2, or 3
 */
void LinMotionEngine_cr3bp::getLiss(int L, double Axy, bool xAmp, double phi, double Az, double psi, Arcset_cr3bp* pArc,
	unsigned int numNodes){

	if(numNodes < 2)
		throw Exception("LinMotionEngine_cr3bp::getLinear: Must specify at least 2 nodes");

	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp*>(pArc->getSysData());

	double mu = pSys->getMu();

	// Locate Lagrange point
	double LPtPos[3] = {0};
	DynamicsModel_cr3bp::getEquilibPt(pSys, L, tol, LPtPos);

	// Get partial derivatives of pseudo potential
	double ddots[6] = {0};
	DynamicsModel_cr3bp::getUDDots(mu, LPtPos[0], LPtPos[1], LPtPos[2], ddots);

	// Out of plane frequency
	double w_z = sqrt(-ddots[2]);

	if(L < 4){
		// Get state size
		const DynamicsModel *model = pSys->getDynamicsModel();
		const unsigned int core_dim = model->getCoreStateSize();
		const unsigned int ctrl_dim = 0;	// No control implemented for CR3BP model(s)
		const unsigned int extra_dim = model->getExtraStateSize();
		const unsigned int full_dim = core_dim + ctrl_dim + pow(core_dim + ctrl_dim, 2) + extra_dim;
		EOM_ParamStruct eomParams(pSys, nullptr);

		// Compute eigenvalues analytically
		double beta1 = 2 - (ddots[0] + ddots[1])/2;
		double beta2 = sqrt(-1*ddots[1]*ddots[0]);
		std::complex<double> Lam1 = -1*beta1 + sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> Lam2 = -1*beta1 - sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> eigenval[] = {sqrt(Lam1), -1.0*sqrt(Lam1), sqrt(Lam2), -1.0*sqrt(Lam2)};

		// Motion for OSC type only
		double s = std::imag(eigenval[2]);
		double beta3 = (s*s + ddots[0])/(2*s);
		double period_xy = 2*PI/s;

		// Initial conditions computed from amplitude/phase angle information
		double xi0 = xAmp ? Axy*cos(phi) : -Axy/beta3 * cos(phi);
		double eta0 = xAmp ? -Axy*beta3*sin(phi) : Axy*sin(phi);

		std::vector<double> state(full_dim, 0);
		double t = 0, dt = revs*period_xy/(numNodes - 1);
		unsigned int nodeCount = 0;
		for(t = 0; t < revs*period_xy; t += t_step){
			state[0] = xi0*cos(s*t) + eta0/beta3*sin(s*t);			// xi
			state[1] = eta0*cos(s*t) - beta3*xi0*sin(s*t);			// eta
			state[2] = Az*sin(w_z*t + psi);							// zeta
			state[3] = -s*xi0*sin(s*t) + eta0*s/beta3 * cos(s*t);	// xi_dot
			state[4] = -s*eta0*sin(s*t) - s*beta3*xi0*cos(s*t);		// eta_dot
			state[5] = -w_z*Az*cos(w_z*t + psi);					// zeta_dot

			state[0] += LPtPos[0];
			state[1] += LPtPos[1];
			state[2] += LPtPos[2];

			storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
		}

		// Add a final node
		Node nodeF(&(state.front()), core_dim, t);
		Segment &lastSeg = pArc->getSegRefByIx(-1);

		// Link things up and update the segment
		int idf = model->sim_addNode(nodeF, &(state.front()), t, pArc, &eomParams, Event_tp::SIM_TOF);
		pArc->getNodeRef(idf).addLink(lastSeg.getID());
		lastSeg.setTerminus(idf);
		lastSeg.appendState(&(state.front()), full_dim);
		lastSeg.appendTime(t);
		lastSeg.updateTOF();
	}else{
		throw Exception("LinMotionEngine::getCR3BPLiss: Cannot compute Lissajous motion for anything other than the collinear points");
	}

	// Compute Jacobi Constant for each step; won't be constant because non-linear dynamics are
	// not enforced, but is still useful information
	for(unsigned int i = 0; i < pArc->getNumNodes(); i++){
		std::vector<double> state = pArc->getStateByIx(i);
		pArc->setJacobiByIx(i, DynamicsModel_cr3bp::getJacobi(&(state[0]), mu));
	}
}//========================================================


/**
 * 	\brief Construct a trajectory from linear approximations of the CR3BP EOMs
 *
 *	This function sets the out-of-plane motion to have zero amplitude
 *
 *	\param L the Lagrange point number [1,5]
 *	\param r0 a three-element vector that specifies the initial position of the arc relative
 *		to the chosen Lagrange point
 *	\param type the type of linearized motion desired
 *	\param pArc pointer to an arcset to store the linear trajectory in
 *	\param numNodes the number of nodes to place around the linear trajectory; 
 *	this is nodes for the entire propagation (period * num revs)
 *
 *	\return a trajectory object containing one revolution of the trajectory. Because this motion
 *	is generated from simplified dynamics, no information about the STM or Jacobi Constant is 
 *	computed. Accelerations are also not computed. These values are all stored as NAN
 */
void LinMotionEngine_cr3bp::getLinear(int L, double r0[3], unsigned int type, Arcset_cr3bp* pArc, unsigned int numNodes){
	return getLinear(L, r0, 0, 0, type, pArc, numNodes);
}//====================================================

/**
 * 	\brief Construct a trajectory from linear approximations of the CR3BP EOMs
 *
 *	\param L the Lagrange point number [1,5]
 *	\param r0 a three-element vector that specifies the initial position of the arc relative
 *		to the chosen Lagrange point
 *	\param Az out-of-plane amplitude, in non-dimensional units
 *	\param psi defines the starting elevation angle of the out-of-plane motion, radians
 *	\param type the type of linearized motion desired
 *	\param pArc pointer to an arcset to store the linear trajectory in
 *	\param numNodes the number of nodes to place around the linear trajectory; 
 *	this is nodes for the entire propagation (period * num revs)
 *
 *	\return a trajectory object containing one revolution of the trajectory. Because this motion
 *	is generated from simplified dynamics, no information about the STM or Jacobi Constant is 
 *	computed. Accelerations are also not computed. These values are all stored as NAN
 *	\throws Exception if the <tt>type</tt> does not correspond with the specified Lagrange
 *	point <tt>L</tt>
 */
void LinMotionEngine_cr3bp::getLinear(int L, double r0[3], double Az, double psi, unsigned int type, Arcset_cr3bp* pArc,
	unsigned int numNodes){

	if(numNodes < 2)
		throw Exception("LinMotionEngine_cr3bp::getLinear: Must specify at least 2 nodes");

	const SysData_cr3bp *pSys = static_cast<const SysData_cr3bp*>(pArc->getSysData());
	double mu = pSys->getMu();

	// Locate Lagrange point
	double LPtPos[3] = {0};
	DynamicsModel_cr3bp::getEquilibPt(pSys, L, tol, LPtPos);

	// Get partial derivatives of pseudo potential
	double ddots[6] = {0};
	DynamicsModel_cr3bp::getUDDots(mu, LPtPos[0], LPtPos[1], LPtPos[2], ddots);

	double xi0 = r0[0];		// Initial x-variation
	double eta0 = r0[1];	// Initial y-variation

	// Out of plane frequency
	double w_z = sqrt(-ddots[2]);

	// Get state size
	const DynamicsModel *model = pSys->getDynamicsModel();
	const unsigned int core_dim = model->getCoreStateSize();
	const unsigned int ctrl_dim = 0;	// No control implemented for CR3BP model(s)
	const unsigned int extra_dim = model->getExtraStateSize();
	const unsigned int full_dim = core_dim + ctrl_dim + pow(core_dim + ctrl_dim, 2) + extra_dim;

	std::vector<double> state(full_dim, 0);
	double t = 0, period_xy = 0, dt = 0;
	unsigned int nodeCount = 0;
	EOM_ParamStruct eomParams(pSys, nullptr);

	if(L < 4){
		// Compute eigenvalues analytically
		double beta1 = 2 - (ddots[0] + ddots[1])/2;
		double beta2 = sqrt(-1*ddots[1]*ddots[0]);
		std::complex<double> Lam1 = -1*beta1 + sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> Lam2 = -1*beta1 - sqrt(beta1*beta1 + beta2*beta2);
		std::complex<double> eigenval[] = {sqrt(Lam1), -1.0*sqrt(Lam1), sqrt(Lam2), -1.0*sqrt(Lam2)};

		double s;
		switch(type){
			case LinMotion_tp::NONE:	// for default behavior
			case LinMotion_tp::OSC:
			{
				s = std::imag(eigenval[2]);
				double beta3 = (s*s + ddots[0])/(2*s);
				period_xy = 2*PI/s;
				dt = revs*period_xy/(numNodes - 1);
				for(t = 0; t < revs*period_xy; t += t_step){
					state[0] = xi0*cos(s*t) + eta0/beta3*sin(s*t) + LPtPos[0];
					state[1] = eta0*cos(s*t) - beta3*xi0*sin(s*t) + LPtPos[1];
					state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
					state[3] = -s*xi0*sin(s*t) + eta0*s/beta3 * cos(s*t);
					state[4] = -s*eta0*sin(s*t) - s*beta3*xi0*cos(s*t);
					state[5] = -w_z*Az*cos(w_z*t + psi);

					storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
				}
				break;
			}
			case LinMotion_tp::HYP:
			{
				s = std::real(eigenval[0]);
				double alpha = (s*s - ddots[0])/(2*s);
				period_xy = 2*PI/s;
				dt = revs*period_xy/(numNodes - 1);
				for(t = 0; t < revs*period_xy; t += t_step){
					state[0] = xi0*cosh(s*t) + eta0/alpha*sinh(s*t) + LPtPos[0];
					state[1] = eta0*cosh(s*t) + alpha*xi0*sinh(s*t) + LPtPos[1];
					
					state[3] = s*xi0*sinh(s*t) + s*eta0/alpha*cosh(s*t);
					state[4] = s*eta0*sinh(s*t) + s*alpha*xi0*cosh(s*t);
					
					state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
					state[5] = -w_z*Az*cos(w_z*t + psi);

					storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
				}
	            break;
			}
			default:
				throw Exception("Invalid type for collinear points");
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
				astrohelion::printErr("LinMotionEngine_cr3bp :: Eigenvalue order appears to be shifting... please hardcode abs()\n");
			}

			switch(type){
				case LinMotion_cr3bp_tp::LPO:
					period_xy = 2*PI/s1d;
					dt = revs*period_xy/(numNodes - 1);
					for(t = 0; t < revs*period_xy; t+= t_step){
						state[0] = xi0*cos(s1d*t) + (eta0 - a1*xi0)/b1 * sin(s1d*t) + LPtPos[0];
						state[1] = eta0*cos(s1d*t) - (b1*xi0 - a1*(eta0 - a1*xi0)/b1)*sin(s1d*t) + LPtPos[1];
						state[3] = -s1d*xi0*sin(s1d*t) + s1d*(eta0 - a1*xi0)/b1 * cos(s1d*t);
						state[4] = -s1d*eta0*sin(s1d*t) - s1d*(b1*xi0 - a1*(eta0 - a1*xi0)/b1)*cos(s1d*t);

						state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
						state[5] = -w_z*Az*cos(w_z*t + psi);

						storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
					}
                	break;
                case LinMotion_tp::NONE: // for default behavior
				case LinMotion_cr3bp_tp::SPO:
					period_xy = 2*PI/s3d;
					dt = revs*period_xy/(numNodes - 1);
					for(t = 0; t < revs*period_xy; t+= t_step){
						state[0] = xi0*cos(s3d*t) + (eta0 - a3*xi0)/b3 * sin(s3d*t) + LPtPos[0];
						state[1] = eta0*cos(s3d*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(s3d*t) + LPtPos[1];
						state[3] = -s3d*xi0*sin(s3d*t) + s3d*(eta0 - a3*xi0)/b3 * cos(s3d*t);
						state[4] = -s3d*eta0*sin(s3d*t) - s3d*(b3*xi0 - a3*(eta0 - a3*xi0)/b3)*cos(s3d*t);

						state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
						state[5] = -w_z*Az*cos(w_z*t + psi);

						storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
					}
                	break;
				case LinMotion_cr3bp_tp::MPO:
				{
					period_xy = s1d < s3d ? 2*PI/s1d : 2*PI/s3d;
					double R1 = xi0/(1 + nu);
                	double C1 = xi0*(a1 + nu*a3)/(2*(1 + nu)*(b1 + nu*b3)) - eta0/(2*(b1 + nu*b3));
                	dt = revs*period_xy/(numNodes - 1);
                	for(t = 0; t < revs*period_xy; t+= t_step){
						state[0] = R1*cos(s1d*t) - 2*C1*sin(s1d*t) + nu*R1*cos(s3d*t) - 2*nu*C1*sin(s3d*t) + LPtPos[0];
						state[1] = (a1*R1 - 2*b1*C1)*cos(s1d*t) - (b1*R1 + 2*a1*C1)*sin(s1d*t) +
                    		nu*(a3*R1 - 2*b3*C1)*cos(s3d*t) - nu*(b3*R1 + 2*a3*C1)*sin(s3d*t) + LPtPos[1];

                    	state[3] = -s1d*R1*sin(s1d*t) - s1d*2*C1*cos(s1d*t) - s3d*nu*R1*sin(s3d*t) - s3d*2*nu*C1*cos(s3d*t);
						state[4] = -s1d*(a1*R1 - 2*b1*C1)*sin(s1d*t) - s1d*(b1*R1 + 2*a1*C1)*cos(s1d*t) +
                    		-s3d*nu*(a3*R1 - 2*b3*C1)*sin(s3d*t) - s3d*nu*(b3*R1 + 2*a3*C1)*cos(s3d*t);
						
						state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
						state[5] = -w_z*Az*cos(w_z*t + psi);

						storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
					}
                	break;
                }
				default:
					throw Exception("Invalid type for triangular points, Case I");
			}
		}else if(std::abs(g) < tol){
			// Case II, g is 0 (approx.)
			double s1d = s1.imag();
			period_xy = 2*PI/s1d;
			dt = revs*period_xy/(numNodes - 1);
			for(t = 0; t < revs*period_xy; t+= t_step){
				state[0] = xi0*cos(s1d*t) - (a1*xi0 - eta0)/b1 * sin(s1d*t) + LPtPos[0];
				state[1] = eta0*cos(s1d*t) - (xi0*std::abs(alpha1) - a1*eta0)/b1 * sin(s1d*t) +  LPtPos[1];

				state[3] = -s1d*xi0*sin(s1d*t) - s1d*(a1*xi0 - eta0)/b1 * cos(s1d*t);
				state[4] = -s1d*eta0*sin(s1d*t) - s1d*(xi0*std::abs(alpha1) - a1*eta0)/b1 * cos(s1d*t);

				state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
				state[5] = -w_z*Az*cos(w_z*t + psi);

				storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
			}
		}else{
			// Case III
			double p = std::real(s1);
			double q = std::imag(s1);
			switch(type){
				case LinMotion_tp::NONE: // for default behavior
				case LinMotion_tp::STAB_OSC:
					period_xy = 2*PI/q;
					dt = revs*period_xy/(numNodes - 1);
					for(t = 0; t < revs*period_xy; t += t_step){
						state[0] = exp(-p*t)*(xi0*cos(q*t) + (eta0 - a3*xi0)/b3 * sin(q*t)) + LPtPos[0];
	                    state[1] = exp(-p*t)*(eta0*cos(q*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(q*t)) + LPtPos[1];
	                    
	                    state[3] = -p*exp(-p*t)*(xi0*cos(q*t) + (eta0 - a3*xi0)/b3 * sin(q*t)) + 
	                    	exp(-p*t)*(-q*xi0*sin(q*t) + q*(eta0 - a3*xi0)/b3 * cos(q*t));
	                    state[4] = -p*exp(-p*t)*(eta0*cos(q*t) - (b3*xi0 - a3*(eta0 - a3*xi0)/b3)*sin(q*t)) +
	                    	exp(-p*t)*(-q*eta0*sin(q*t) - q*(b3*xi0 - a3*(eta0 - a3*xi0)/b3)*cos(q*t));

	                    state[2] = Az*sin(w_z*t + psi) + LPtPos[2];
						state[5] = -w_z*Az*cos(w_z*t + psi);

						storeData(model, pArc, t, &state, nodeCount, &eomParams, dt, t_step, core_dim, full_dim);
					}
					break;	
				case LinMotion_tp::UNSTAB_OSC:
					throw Exception("Triangular points, Case III: UNSTAB_OSC not yet implemented");
					break;
				default:
					throw Exception("Invalid type for triangular points, Case III");
			}
		}
	}

	// Add a final node
	Node nodeF(&(state.front()), core_dim, t);
	Segment &lastSeg = pArc->getSegRefByIx(-1);

	// Link things up and update the segment
	int idf = model->sim_addNode(nodeF, &(state.front()), t, pArc, &eomParams, Event_tp::SIM_TOF);
	pArc->getNodeRef(idf).addLink(lastSeg.getID());
	lastSeg.setTerminus(idf);
	lastSeg.appendState(&(state.front()), full_dim);
	lastSeg.appendTime(t);
	lastSeg.updateTOF();
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

void LinMotionEngine_cr3bp::reset(){
	LinMotionEngine::reset();

	nu = 1;
}//====================================================

}// END of Astrohelion namespace