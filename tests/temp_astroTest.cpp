#include "AllIncludes.hpp"

#include <algorithm>
#include <chrono>
#include <exception>
#include <iostream>
#include "matio.h"
#include <random>
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

std::vector<double> concat(std::vector<double> &a, std::vector<double> &b){
	a.insert(a.end(), b.begin(), b.end());
	return a;
}

#pragma omp declare reduction \
(concatVecs : std::vector<double> : omp_out=concat(omp_out, omp_in))\
initializer(omp_priv = std::vector<double>{} )

int main(int argc, char** argv){
	
	if(argc != 3){
		printErr("require 2 inputs: <num segs> <num iterations>");
		return EXIT_SUCCESS;
	}

	const unsigned int numSegs = atoi(argv[1]);
	const unsigned int numIts = atoi(argv[2]);
	int onePC = numIts > 100 ? std::floor(numIts/100) : 1;
	const int printAtIt = std::min(onePC, 500);

	const SysData_cr3bp_lt sys("earth", "moon", 100);
	const double f = 1e-2, Isp = 1500;
	const double tof = 2*PI;

	const double seg_tof = tof/static_cast<double>(numSegs);

	// Law for inertially-fixed thrust vector
	unsigned int lawID = ltlaw::GEN_INERT | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	std::vector<double> params {0, f, Isp};
	ControlLaw_cr3bp_lt inertLaw(lawID, params);

	// Law for rotating-frame-fixed thrust vector
	lawID = ltlaw::GENERAL | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	params = std::vector<double> {f, Isp};
 	ControlLaw_cr3bp_lt rotLaw(lawID, params);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	std::vector<double> data;

	// obtain a seed from the system clock:
  	const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	int prog = 0, threadIt = 0;

  	// Random number generators with uniform distributions
	std::default_random_engine rand(seed);
	std::uniform_real_distribution<double> angleGen(-PI, PI);
	std::uniform_real_distribution<double> posGen(-1.2, 1.2);
	std::uniform_real_distribution<double> velGen(-0.75, 0.75);

	#pragma omp parallel for firstprivate(sim, rand, angleGen, posGen, velGen,\
		inertLaw, rotLaw, threadIt) schedule(dynamic) reduction(concatVecs:data)
	for(unsigned int i = 0; i < numIts; i++){
		std::vector<double> q0(7, 0), rotCtrl0(2, 0);

		if(threadIt == 0)
			rand.seed(seed + omp_get_thread_num());

		// Generate random initial states
		q0[0] = posGen(rand);
		q0[1] = posGen(rand);
		q0[3] = velGen(rand);
		q0[4] = velGen(rand);
		q0[6] = 1;	// always begin with 100% mass

		rotCtrl0[0] = angleGen(rand);
		double theta0 = angleGen(rand);
		double t0 = velGen(rand);

		std::vector<double> p {theta0, f, Isp};
		inertLaw.setParams(p);

		Arcset_cr3bp_lt arc_rotFixed(&sys), arc_inertFixed(&sys);

		// Single propagation for arc with thrust fixed in rotating frame
		try{
			sim.runSim_manyNodes(q0, rotCtrl0, t0, tof, numSegs+1, &arc_rotFixed, &rotLaw);
		}catch(const Exception &e){
			#pragma omp atomic
			prog++;

			continue;
		}

		if(arc_rotFixed.getNodeRefByIx(-1).getTriggerEvent() != Event_tp::SIM_TOF){
			// printErr("  Simulation with rot-fixed did not end as expected\n");

			#pragma omp atomic
			prog++;

			continue;	// don't proceed with solutions that crash or have other errors
		}

		// Generate an arcset with several propagations with inertially-fixed
		// thrust to approximate the rotating solution
		Arcset_cr3bp_lt temp(&sys);
		std::vector<double> qi = q0;
		bool completedInertArc = true;
		for(unsigned int s = 0; s < numSegs; s++){
			std::vector<double> inertCtrl0 {rotCtrl0[0] + theta0 + 
				t0 + s*seg_tof + 0.5*seg_tof, 0};

			temp.reset();
			try{
				sim.runSim(qi, inertCtrl0, t0 + s*seg_tof, seg_tof, &temp, &inertLaw);
			}catch(const Exception &e){
				#pragma omp atomic
				prog++;

				completedInertArc = false;
				break;
			}
			
			if(temp.getNodeRefByIx(-1).getTriggerEvent() != Event_tp::SIM_TOF){
				// printErr("  Simulation with inert-fixed did not end as expected\n");
				completedInertArc = false;
				break;
			}
			qi = temp.getStateByIx(-1);
			if(s == 0){
				arc_inertFixed = temp;
			}else{
				arc_inertFixed.appendSetAtNode(&temp, 
					arc_inertFixed.getNodeRefByIx(-1).getID(),
					temp.getNodeRefByIx(0).getID(), 0);
			}
		}

		if(!completedInertArc){

			#pragma omp atomic
			prog++;

			continue;
		}

		arc_inertFixed.putInChronoOrder();
		// arc_inertFixed.print();

		// Compute differences between the two solutions
		std::vector<double> diffq(6*numSegs, 0);
		std::vector<double> diffH(numSegs+1, 0);
		for(unsigned int s = 0; s < numSegs; s++){
			std::vector<double> q_rot = arc_rotFixed.getSegRefByIx(s).getStateByRow(-1);
			std::vector<double> q_inert = arc_inertFixed.getSegRefByIx(s).getStateByRow(-1);

			for(unsigned int i = 0; i < 6; i++){
				diffq[6*s + i] = q_inert[i] - q_rot[i];
			}

			double H_rot = DynamicsModel_cr3bp_lt::getHamiltonian(
				arc_rotFixed.getSegRefByIx(s).getTimeByIx(-1), &(q_rot[0]),
				&sys, &rotLaw);
			double H_inert = DynamicsModel_cr3bp_lt::getHamiltonian(
				arc_inertFixed.getSegRefByIx(s).getTimeByIx(-1), &(q_inert[0]),
				&sys, &inertLaw);

			diffH[s+1] = H_inert - H_rot;
		}

		// additionally compute the Hamiltonian difference at the initial state
		std::vector<double> q_rot = arc_rotFixed.getSegRefByIx(0).getStateByRow(0);
		std::vector<double> q_inert = arc_inertFixed.getSegRefByIx(0).getStateByRow(0);
		double H_rot = DynamicsModel_cr3bp_lt::getHamiltonian(
			arc_rotFixed.getSegRefByIx(0).getTimeByIx(0), &(q_rot[0]),
			&sys, &rotLaw);
		double H_inert = DynamicsModel_cr3bp_lt::getHamiltonian(
			arc_inertFixed.getSegRefByIx(0).getTimeByIx(0), &(q_inert[0]),
			&sys, &inertLaw);
		diffH[0] = H_inert - H_rot;

		// save data to output vector
		data.insert(data.end(), q0.begin(), q0.end());			// 7
		data.push_back(t0);										// 1
		data.push_back(theta0);									// 1
		data.push_back(numSegs);								// 1
		data.push_back(seg_tof);								// 1
		data.insert(data.end(), diffq.begin(), diffq.end());	// 6*numSegs
		data.insert(data.end(), diffH.begin(), diffH.end());	// numSegs+1

		arc_rotFixed.saveToMat("arc_rotFixed.mat");
		arc_inertFixed.saveToMat("arc_inertFixed.mat");

		#pragma omp atomic
		prog++;

		if(prog % printAtIt == 0){
			std::cout << 100*(static_cast<double>(prog)/static_cast<double>(numIts)) << "%\n";
		}

		threadIt++;
	}

	double rowWid = 7 + 4 + 6*numSegs + numSegs+1;
	char filename[128];
	sprintf(filename, "montecarlo_results_f%02.1e_Isp%d_segs%u.mat", f, 
		static_cast<int>(Isp), numSegs);
	saveMatrixToFile(filename, "Results", data, data.size()/rowWid, rowWid);

	return EXIT_SUCCESS;
}