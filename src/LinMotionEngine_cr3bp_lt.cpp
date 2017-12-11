/**
 *	\file LinMotionEngine_cr3bp_lt.cpp
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

#include "LinMotionEngine_cr3bp_lt.hpp"

#include <Eigen/Eigenvalues>

#include "Arcset_cr3bp_lt.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "DynamicsModel_cr3bp_lt.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  \brief Default, do-nothing constructor
 */
LinMotionEngine_cr3bp_lt::LinMotionEngine_cr3bp_lt() : LinMotionEngine(){}

//-----------------------------------------------------------------------------
//      Linear Motion Generation Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Compute the a trajectory in the linear dynamics
 *  \details [long description]
 * 
 *  \param eqPt 3D position vector of the equilibrium solution. This may be obtained 
 *  from DynamicsModel_cr3bp_lt::getEquilibPt
 *  \param f thrust magnitude (nondimensional) that corresponds to the equilibrium solution
 *  \param alpha thrust pointing angle (in plane, radians) that corresponds to the 
 *  equilibrium solution
 *  \param x0 Initial variational state (nondimensional, relative to eqPt)
 *  \param motionTp describes the type of motion
 *  \param pArc pointer to an arcset in which to store the linear trajectory
 *  \param pLaw pointer to a control law object; the control used on the linear trajectory is 
 *  stored in this object.
 *  \param numNodes How many nodes to create on the trajectory. Must be >= 2
 */
void LinMotionEngine_cr3bp_lt::getLinear(double eqPt[3], double f, double alpha, double x0[3],
	unsigned int motionTp, Arcset_cr3bp_lt *pArc, ControlLaw_cr3bp_lt *pLaw, unsigned int numNodes){

	if(numNodes < 2)
		throw Exception("LinMotionEngine_cr3bp_lt::getLinear: Must specify at least 2 nodes");

	const SysData_cr3bp_lt *pSys = static_cast<const SysData_cr3bp_lt*>(pArc->getSysData());
	double mu = pSys->getMu();

	// Get the partial derivatives of the pseudo potential
	double ddots[6] = {0};
	DynamicsModel_cr3bp::getUDDots(mu, eqPt[0], eqPt[1], eqPt[2], ddots);

	// Construct the A matrix (planar)
	MatrixXRd A = MatrixXRd::Zero(4,4);
	A(0,2) = 1;	A(1,3) = 1;		// Identity in top-right
	A(2,0) = ddots[0];	A(2,1) = ddots[3];	// partial derivatives of pseudo-potential
	A(3,0) = ddots[3];	A(3,1) = ddots[1];
	A(2,3) = 2;			A(3,2) = -2;		// Velocity terms

	// toCSV(A, "data/temp_A.csv");

	// Compute the eigenvalues and eigenvectors
	Eigen::EigenSolver<MatrixXRd> eigensolver(A);
	if(eigensolver.info() != Eigen::Success)
		throw Exception("LinMotionEngine_cr3bp_lt::getLinear: Unable to compute eigenvalues of A matrix");

	Eigen::VectorXcd vals = eigensolver.eigenvalues();
	std::vector<cdouble> eigVals = std::vector<cdouble>(vals.data(), vals.data()+6);
	MatrixXRcd eigVecs = eigensolver.eigenvectors();

	// toCSV(eigVecs, "data/temp_eigVecs.csv");
	// std::cout << "Eigenvectors:\n" << eigVecs << std::endl;

	bool pureReal[4], pureImag[4], mixed[4];
	int numPureReal = 0, numPureImag = 0, numMixed = 0;
	for(unsigned int i = 0; i < 4; i++){
		if(verbosity >= Verbosity_tp::ALL_MSG)
			printf("Eigenvalue: %f + %fj\n", real(eigVals[i]), imag(eigVals[i]));
		
		pureReal[i] = std::abs(imag(eigVals[i])) < 1e-12;
		pureImag[i] = std::abs(real(eigVals[i])) < 1e-12;
		mixed[i] = !pureReal[i] & !pureImag[i];

		numPureReal += pureReal[i];
		numPureImag += pureImag[i];
		numMixed += mixed[i];
	}

	// Create the control law object that describes the constant thrust vector with a very high Isp
	std::vector<double> params {f, 10000};
	*pLaw = ControlLaw_cr3bp_lt(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, params);
	EOM_ParamStruct eomParams(pSys, pLaw);

	// Get state size
	const DynamicsModel *model = pSys->getDynamicsModel();
	const unsigned int core_dim = model->getCoreStateSize();
	const unsigned int ctrl_dim = pLaw->getNumStates();
	const unsigned int extra_dim = model->getExtraStateSize();
	const unsigned int full_dim = core_dim + ctrl_dim + pow(core_dim + ctrl_dim, 2) + extra_dim;

	switch(motionTp){
		case LinMotion_tp::HYP:
			if(verbosity > Verbosity_tp::NO_MSG)
				printErr("LinMotionEngine_cr3bp_lt::getLinear: Hyperbolic mode is not implemented\n");
			break;
		case LinMotion_tp::NONE:
		case LinMotion_tp::OSC:
		case LinMotion_tp::STAB_OSC:
		case LinMotion_tp::UNSTAB_OSC:
		{
			if(numPureImag < 2){
				if(verbosity > Verbosity_tp::NO_MSG)
					printErr("LinMotionEngine_cr3bp_lt::getLinear: There is not a center subspace\n");
			}else{
				// Find the largest (i.e., positive) frequency
				double freq = 0;
				unsigned int ix = 0;
				for(unsigned int i = 0; i < 4; i++){
					if(!pureReal[i] && imag(eigVals[i]) > freq){
						freq = imag(eigVals[i]);
						ix = i;
					}
				}

				Eigen::VectorXcd eigVec = eigVecs.col(ix);
				std::vector<cdouble> vec(eigVec.data(), eigVec.data()+4);
				std::vector<double> u = real(vec);
				std::vector<double> w = imag(vec);
				double c4 = (x0[1]*u[0] - x0[0]*u[1])/(w[1]*u[0] - w[0]*u[1]);
				double c3 = (x0[0] - c4*w[0])/u[0];
				double rate = std::abs(real(eigVals[ix]))*(motionTp == LinMotion_tp::STAB_OSC ? -1 : 1);
				if(std::abs(rate) < 1e-14)
					rate = 0;

				if(motionTp == LinMotion_tp::OSC && std::abs(rate) > 1e-12){
					if(verbosity > Verbosity_tp::NO_MSG)
						printErr("LinMotionEngine_cr3bp_lt::getLinear: Requested oscillatory motion, but the oscillatory modes are spirals; aborting\n");
					return;
				}
				if((motionTp == LinMotion_tp::STAB_OSC || motionTp == LinMotion_tp::UNSTAB_OSC) && std::abs(rate) < 1e-12){
					if(verbosity > Verbosity_tp::NO_MSG)
						printWarn("LinMotionEngine_cr3bp_lt::getLinear: No stable/unstable component to oscillatory motion; proceeding with pure oscillatory\n");
				}

				double period = 2*PI/freq;

				std::vector<double> state(full_dim, 0);
				state[6] = 1;					// Set mass = 1 for all time
				state[core_dim] = alpha;		// Set alpha, beta (core_dim + 1) is always zero
				double t = 0, dt = revs*period/(numNodes - 1);
				unsigned int nodeCount = 0;
				for(t = 0; t < revs*period; t += t_step){
					for(unsigned int i = 0; i < 6; i++){
						// Leave zeta and zeta_dot as zeros
						int ix = i - std::floor(i/3);
						if(i != 2 && i < 5)
							state[i] = c3*exp(rate*t)*(u[ix]*cos(freq*t) - w[ix]*sin(freq*t)) + c4*exp(rate*t)*(w[ix]*cos(freq*t) + u[ix]*sin(freq*t));
						else
							state[i] = 0;
					}
					state[0] += eqPt[0];	// Shift to be centered around equilibrium solution
					state[1] += eqPt[1];
					state[2] += eqPt[2];
					

					if(t >= dt*nodeCount){
						// Add a node
						Node node0(&(state.front()), core_dim, t);
						int ID = model->sim_addNode(node0, &(state.front()), t, pArc, &eomParams, Event_tp::NONE);
						// Link the previous segment, if it exists
						if(nodeCount > 0){
							Segment &lastSeg = pArc->getSegRefByIx(-1);
							lastSeg.setTerminus(ID);

							lastSeg.appendState(&(state.front()), full_dim);
							lastSeg.appendTime(t);
							lastSeg.updateTOF();
						}

						// Create a new segment
						Segment seg(ID, Linkable::INVALID_ID, t_step);
						seg.appendState(&(state.front()), full_dim);
						seg.setStateWidth(full_dim);
						seg.appendTime(t);
						model->sim_addSeg(seg, &(state.front()), t, pArc, &eomParams);
						nodeCount++;
					}else{
						Segment &lastSeg = pArc->getSegRefByIx(-1);
						lastSeg.appendState(&(state.front()), full_dim);
						lastSeg.appendTime(t);
					}
				}// End of loop to construct trajectory

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
			}

			break;
		}
		default:
			throw Exception("LinMotionEngine_cr3bp_lt::getLinear: Unrecognized motion type");
	}
}//====================================================

}// END of astrohelion namespace






