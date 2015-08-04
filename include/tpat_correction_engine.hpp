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

#ifndef H_CORRECTIONS
#define H_CORRECTIONS

#include "tpat_nodeset.hpp"
#include "tpat_sys_data.hpp"

#include <vector>

// Forward declarations
class tpat_nodeset_bcr4bpr;
class tpat_sys_data_bcr4bpr;
class tpat_constraint;
class tpat_node;
class tpat_nodeset_cr3bp;
class tpat_matrix;
class tpat_simulation_engine;

/**
 *	@brief a custom data class to encapsulate data used in each iteration
 *	of the corrections process.
 *
 *	This data object can be passed to other functions, allowing us to break 
 *	the master corrections loop into smaller functions without requiring an
 *	obscene amount of arguments to be passed in.
 */
struct iterationData{
	public:
		tpat_sys_data *sysData;					//!< A pointer to the system data object used for this corrections process
		std::vector<double> X;					//!< Free-Variable Vector
		std::vector<double> FX;					//!< Constraint Function Vector
		std::vector<double> DF;					//!< Jacobian Matrix
		std::vector<double> deltaVs;			//!< nx3 vector of non-dim delta-Vs
		std::vector<double> primPos;			//!< Store the positions of the primaries
		std::vector<double> primVel;			//!< Store the velocities of the primaries
		std::vector<tpat_node> origNodes;		//!< Store the original nodes that birthed this correction process
		std::vector<tpat_traj> allSegs;			//!< A collection of all integrated segments
		std::vector<tpat_constraint> allCons;	//!< A list of all constraints
		std::vector<int> slackAssignCon;		//!< Indices of constraints, index of entry corresponds to a slack variable
		std::vector<int> conRows;				//!< Each entry holds the row # for the constraint; i.e. 0th element holds row # for 0th constraint

		int numNodes = 0;			//!< Number of nodes in the entire nodeset
		int count = 0;				//!< Count of number of iterations through corrections process

		int numSlack = 0;			//!< # slack variables
		int totalCons = 0;			//!< Total # constraints -> # rows of DF
		int totalFree = 0;			//!< Total # free var. -> # cols of DF

		bool varTime = true;		//!< Whether or not the simulation is using variable time
};

/**
 *	@brief An engine object to perform corrections, such as single-shooting
 *	and multiple shooting algorithms
 *
 *	@author Andrew Cox
 *	@version August 3, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_correction_engine{
	public:
		// *structors
		/** Default, do-nothing constructor */
		tpat_correction_engine(){}
		tpat_correction_engine(const tpat_correction_engine&);
		~tpat_correction_engine();

		// Operators
		tpat_correction_engine& operator =(const tpat_correction_engine &e);

		// Set and get functions
		bool usesVarTime() const;
		int getMaxIts() const;
		double getTol() const;
		bool isVerbose() const;
		bool isFindingEvent() const;
		tpat_nodeset_cr3bp getCR3BP_Output();
		tpat_nodeset_bcr4bpr getBCR4BPR_Output();

		void setVarTime(bool);
		void setMaxIts(int);
		void setTol(double);
		void setVerbose(bool);
		void setFindEvent(bool);

		// Utility/Action functions
		void correct_cr3bp(tpat_nodeset_cr3bp*);
		void correct_bcr4bpr(tpat_nodeset_bcr4bpr*);

	private:
		/** Whether or not to spit out lots of messages */
		bool verbose = false;

		/** Whether or not to use variable time in the corrections process */
		bool varTime = true;

		/** Maximum number of iterations before giving up */
		int maxIts = 20;

		/** Maximum acceptable error value, non-dimensional units */
		double tol = 1e-12;

		/** Whether or not an input nodeset was supplied and corrected */
		bool receivedNodesetIn = false;

		/** Whether or not space has been dynamically allocated for nodeset_out */
		bool createdNodesetOut = false;

		/** Flag to turn on when this algorithm is being used to locate an event */
		bool findEvent = false;

		/** The input nodeset, not modified*/
		tpat_nodeset *nodeset_in = 0;

		/** The output nodeset, constructed from the corrected arcs */
		tpat_nodeset *nodeset_out = 0;

		/** Whether or not the engine is ready to be cleaned and/or deconstructed */
		bool isClean = true;

		void copyEngine(const tpat_correction_engine&);
		void correct(tpat_nodeset*);

		// void createPosVelCons(iterationData*, tpat_sys_data::system_t, int);
		// void targetState(iterationData*, tpat_constraint, int, int);
		// void targetMatchAll(iterationData*, tpat_constraint, int, int);
		// void targetMatchCust(iterationData*, tpat_constraint, int, int);
		// void targetDist(iterationData*, tpat_constraint, tpat_sys_data*, int, int);
		// void targetSP(iterationData*, tpat_sys_data_bcr4bpr*, int, int);
		// void targetJC(iterationData*, tpat_constraint, tpat_sys_data*, int, int);
		// void updateDeltaVCon(iterationData*, tpat_sys_data*, int, int);
		// void updatePrimPos(iterationData*, tpat_sys_data*, double);
		// void updatePrimVel(iterationData*, tpat_sys_data*, double);

		tpat_matrix solveUpdateEq(iterationData*);

		void cleanEngine();
};

#endif