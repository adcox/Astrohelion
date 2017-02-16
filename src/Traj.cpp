/**
 *  @file Traj.cpp
 *	@brief Stores information about a trajectory
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "Traj.hpp"

#include "Exceptions.hpp"
#include "Node.hpp"
#include "Nodeset.hpp"
#include "SimEngine.hpp"
#include "Utilities.hpp"

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param data a pointer to a system data object
 */
Traj::Traj(const SysData *data) : BaseArcset(data) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
Traj::Traj(const Traj &t) : BaseArcset(t) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
Traj::Traj(const BaseArcset &a) : BaseArcset(a) {
	initExtraParam();
}//====================================================

/**
 *  @brief Default destructor
 */
Traj::~Traj(){}

/**
 *	@brief Create a trajectory from a nodeset
 *
 *	This algorithm will concatenate trajectories integrated from each node in 
 *	the nodeset. It does not check to make sure the arcs are continuous; that
 *	is up to you. The trajectory is constructed via a simulation engine that ignores
 *	crashes as we assume the initial nodeset has been propagated to either ignore
 *	or avoid the primaries; will not challenge that behavior. Each node is integrated
 *	for the associated time-of-flight and added (via operator +()) to a trajectory object.
 *
 *	@param set a nodeset
 *	@return a trajectory formed from the integrated nodeset
 */
Traj Traj::fromNodeset(Nodeset set){
	SimEngine simEngine;
	simEngine.setMakeCrashEvents(false);	// don't trigger crashes; assume this has been taken care of already
	Traj totalTraj(set.getSysData());

	if(!set.isInChronoOrder())
		set.putInChronoOrder();

	int coreSize = set.getSysData()->getDynamicsModel()->getCoreStateSize();
	MatrixXRd prevSTM;

	for(int s = 0; s < set.getNumSegs(); s++){
		double tof = set.getSegByIx(s).getTOF();
		simEngine.setRevTime(tof < 0);
		Node origin = set.getNode(set.getSegByIx(s).getOrigin());
		Traj temp(set.getSysData());

		bool cont = true;
		std::vector<bool> velCon = set.getSegByIx(s).getVelCon();
		for(unsigned int i = 0; i < velCon.size(); i++){
			if(!velCon[i]){
				cont = false;
				break;
			}
		}


		// Initialize the STM to identity if this is the first segment or if there is a discontinuity
		// Otherwise, maintain continuity by beginning with the previous STM
		MatrixXRd stm0 = (s == 0 || !cont) ? MatrixXRd::Identity(coreSize, coreSize) : prevSTM;
		// printf("(Cont = %c) STM0 = \n", cont ? 'Y' : 'N');
		// for(int r = 0; r < 6; r++){
		// 	for(int c = 0; c < 6; c++){
		// 		printf("%10.2f", stm0(r,c));
		// 	}
		// 	printf("\n");
		// }
		simEngine.runSim(origin.getState(), stm0, origin.getEpoch(), tof, &temp);

		if(s == 0){
			totalTraj = temp;
		}else{
			// Shift the time on the newly propagated segment so it starts where the previous segment left off
			temp.shiftAllTimes(totalTraj.getEpochByIx(-1));

			// Use += so that each piece is put into chronological order, even though this significantly increases run time
			totalTraj += temp;
			// totalTraj.appendSetAtNode(&temp, totalTraj.getNodeByIx(-1).getID(), 0, 0);
		}

		// temp.saveToMat("fromNodeset_tempTraj.mat");
		// totalTraj.saveToMat("fromNodeset_totalTraj.mat");
		// waitForUser();

		prevSTM = set.getSegByIx(s).getSTM();
	}

	// Summation of trajectories jumbles up the order a bit, so fix it before returning
	totalTraj.putInChronoOrder();
	return totalTraj;
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object
 *  @return a pointer to the newly created trajectory
 */
baseArcsetPtr Traj::create( const SysData *sys) const{
	return baseArcsetPtr(new Traj(sys));
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned trajectory
 */
baseArcsetPtr Traj::clone() const{
	return baseArcsetPtr(new Traj(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  @brief Combine two trajectories.
 *  @details This function concatenates two trajectory objects. It is assumed
 *  that the first state on <tt>rhs</tt> is identical to the final state on
 *  <tt>rhs</tt>. The <tt>rhs</tt> object is also assumed to occur after
 *  (chronologically) <tt>lhs</tt>
 * 
 *  @param lhs reference to a trajectory object
 *  @param rhs reference to a trajectory object
 * 
 *  @return the concatenation of lhs + rhs.
 */
Traj operator +(const Traj &lhs, const Traj &rhs){
	const Traj lhs_cpy(lhs);
	const Traj rhs_cpy(rhs);
	Traj result(lhs.pSysData);

	BaseArcset::sum(&lhs, &rhs, &result);

	return result;
}//====================================================

/**
 *  @brief Concatenate this object with another trajectory
 * 
 *  @param rhs reference to a trajectory object
 *  @return the concatenation of this and <tt>rhs</tt>
 *  @see operator +()
 */
Traj& Traj::operator +=(const Traj &rhs){
	Traj temp = *this + rhs;
	copyMe(temp);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Compute the total TOF for this trajectory
 *  @return the total TOF for this trajectory, units consistent
 *  with time units of the SysData object for this trajectory.
 *  If the function is unable to determine the initial or final epoch,
 *  a value of NAN is returned.
 */
double Traj::getTotalTOF() const{
	// Sort the trajectory into chronological order if it isn't already
	// (unlikely to be out of order, but possible)
	if(!bInChronoOrder){
		std::vector<ArcPiece> pieces = getChronoOrder();
		
		// Find the first and final nodes, get their epochs
		double T0 = NAN, Tf = NAN;
		for(unsigned int i = 0; i < pieces.size(); i++){
			if(pieces[i].type == ArcPiece::Piece_tp::NODE){
				T0 = getEpoch(pieces[i].id);
				// printf("T0 = %.4f at node id %d\n", T0, pieces[i].id);
				break;
			}
		}
		for(int i = static_cast<int>(pieces.size()) - 1; i >= 0; i--){
			if(pieces[i].type == ArcPiece::Piece_tp::NODE){
				Tf = getEpoch(pieces[i].id);
				// printf("Tf = %.4f at node id %d\n", Tf, pieces[i].id);
				break;
			}
		}

		// Return the difference so long as an initial and final 
		// epoch have been identified
		if(!std::isnan(T0) && !std::isnan(Tf))
			return Tf - T0;
		else
			return NAN;
	}else{
		return getEpochByIx(-1) - getEpochByIx(0);
	}
}//====================================================

/**
 *	@brief Retrieve the time along the trajectory at a specific step
 *	@param ix node index; if < 0, it will count backwards from end of trajectory
 *	@return the non-dimensional time along the trajectory at the specified step
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
double Traj::getTimeByIx(int ix) const {
	if(ix < 0)
		ix += nodes.size();
	
	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj::getTimeByIx: invalid index");

	return nodes[ix].getEpoch();
}//====================================================

/**
 *  @brief Set the time associated with a node
 *  @details [long description]
 * 
 *	@param ix node index; if < 0, it will count backwards from end of trajectory
 *  @param t time associated with the node
 *  @throws Exception if <tt>ix</tt> is out of bounds
 */
void Traj::setTimeByIx(int ix, double t){
	if(ix < 0)
		ix += nodes.size();
	
	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj::setTimeByIx: invalid index");

	nodes[ix].setEpoch(t);
}//====================================================

/**
 *  @brief Shift all time values by a constant amount
 *  @details This can be useful for use with the EM2SE and SE2EM functions
 * 
 *  @param amount a constant, non-dimensional time shift to apply to 
 *  all time values for points on this trajectory
 */
void Traj::shiftAllTimes(double amount){
	for(unsigned int i = 0; i < nodes.size(); i++){
		nodes[i].setEpoch(nodes[i].getEpoch() + amount);
	}
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Discretize a trajectory into a set of nodes without using integration
 *
 *	This method uses the existing steps in a trajectory to create nodes. Arc segments
 *	are created with equal number of steps regardless of the time or arc length separating
 *	the nodes.
 *
 *	@param numNodes number of nodes, including both the initial and final states on 
 *	the trajectory, which are always included
 *	@return a nodeset with the specified number of nodes
 *	@throws Exception if <tt>numNodes</tt> is less than two
 */
Nodeset Traj::discretize(int numNodes) const{
	if(numNodes < 2)
		throw Exception("Traj::discretize: Cannot split a trajectory into fewer than 2 nodes");

	if(numNodes > static_cast<int>(nodes.size())){
		astrohelion::printWarn("Traj::discretize: User requested more nodes than there are states; returning one node per step, will not meet requested number of nodes\n");
		numNodes = nodes.size();
	}

	double stepSize = static_cast<double>(nodes.size()-1)/static_cast<double>(numNodes - 1.0);

	Nodeset nodeset(pSysData);
	int n = 0;
	while(n < numNodes){
		// Round the step number
		int ix = std::floor(n*stepSize);
		int prevIx = std::floor((n-1)*stepSize);

		// Create a node from this step
		nodeset.addNode(nodes[ix]);

		if(n > 0){
			double tof = getEpochByIx(ix) - getEpochByIx(prevIx);
			nodeset.addSeg(Segment(n-1, n, tof));
		}

		n++;
	}

	return nodeset;
}//=================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void Traj::saveToMat(const char* filename) const{
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		astrohelion::printErr("Traj::saveToMat: Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveEpoch(matfp, "Time");
		saveTOF(matfp, "TOFs");
		saveSTMs(matfp);
		pSysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Print a useful message describing this trajectory to the standard output
 */
void Traj::print() const {
	printf("%s Trajectory\n", pSysData->getTypeStr().c_str());
	for(int i = 0; i < pSysData->getNumPrimaries(); i++)
		printf("  Primary %d: %s\n", i, pSysData->getPrimary(i).c_str());

	printf("  %d Nodes\n  %d Segments\n", getNumNodes(), getNumSegs());
	printf("  TOF = %.4f\n", getTotalTOF());
}//====================================================

/**
 *	@brief Initialize the extra param vector for trajectory-specific info
 */
void Traj::initExtraParam(){
	// Nothing to do here!
}//====================================================

/**
 *  @brief Populate data in this trajectory from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 *  @throws Exception if the data file cannot be opened
 */
void Traj::readFromMat(const char *filepath){

	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Traj: Could not load data from file");
	}
	initNodesSegsFromMat(matfp, "State");
	readStateFromMat(matfp, "State");
	readAccelFromMat(matfp);
	readSTMFromMat(matfp);
	readEpochFromMat(matfp, "Time");
	readTOFFromMat(matfp, "TOFs");

	Mat_Close(matfp);
}//====================================================



}// END of Astrohelion namespace