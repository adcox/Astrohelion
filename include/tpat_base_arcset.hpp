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

#ifndef H_TPAT_BASE_ARCSET
#define H_TPAT_BASE_ARCSET

#include "tpat.hpp"

#include "tpat_constraint.hpp"
#include "tpat_node.hpp"
#include "tpat_segment.hpp"
#include "tpat_sys_data.hpp"

#include "matio.h"
#include <vector>

// Forward Declarations
struct tpat_arc_piece;

/**
 *	@brief Abstract class that provides the framework for trajectories and nodesets
 *	
 *	The arcset object specifies default and mandatory behaviors for all derivative
 *	classes (i.e. tpat_traj and tpat_nodeset). All variables and data for an arc or 
 *	one of its derivative classes are declared and stored here; in other words, no 
 *	derivative classes declare class-specific data objects. This architecture has been
 *	chosen to facilitate easy casting between model-specific derivative classes with
 *	the added bonus of being able to cast easily between, say, a trajectory and a
 *	nodeset.
 *	
 *	This class contains all information about any trajectory or nodeset in a few objects:
 *	* steps - a vector of tpat_arc_step objects, each of which contains information about
 *		the state, acceleration, STM, any any other parameter values at one single step
 *	* sysData - a pointer to a system data object that describes the system this arc
 *		has been generated in
 *	* numExtraParam and extraParamRowSize describe the number of extra parameters and 
 *		the number of elements in each parameter; this various for different dynamical
 *		model-specific derivative classes
 *	* tol - The maximum numerical tolerance with which the data in this object has been 
 *		computed
 *	
 *	The following behavior is mandatory for all derivative classes:
 *	* Ability to add two arcs together via operator +()
 *	* Ability to save to a matlab file via saveToMat()
 *	* Ability to display a textual representation of the object via print()
 *	
 *	Additionally, the following behavior is defined for all derivative classes, though
 *	they may override the default:
 *	* Assignment operator
 *	* Access to position and velocity values at each step via getState()
 *	* Access to acceleration values at each step via getAccel()
 *	* Access to any extra parameters that evolve each step via getExtraParam()
 *	* Access to the STM at each step via getSTM()
 *	* Access to individual step objects via getStep()
 *	* Access to the system data object pointer that describes the system this arc was integrated in
 *
 *	@author Andrew Cox
 *	@version April 28, 2016
 *	@copyright GNU GPL v3.0	
 */
class tpat_base_arcset : public tpat{

public:
	// *structors
	tpat_base_arcset(const tpat_sys_data*);
	tpat_base_arcset(const tpat_base_arcset&);
	virtual tpat_base_arcset* create( const tpat_sys_data* ) const = 0;	//!< Virtual constructor for creation
	virtual tpat_base_arcset* clone() const = 0;							//!< Virtual constructor for copying

	virtual ~tpat_base_arcset();

	// Operators
	tpat_base_arcset& operator =(const tpat_base_arcset&);
	static void sum(const tpat_base_arcset*, const tpat_base_arcset*, tpat_base_arcset*);

	// Set and Get functions
	void addConstraint(tpat_constraint);
	int addNode(tpat_node);
	int addSeg(tpat_segment);
	int appendSetAtNode(const tpat_base_arcset*, int, int, double);
	void clearArcConstraints();
	void clearAllConstraints();
	void deleteNode(int);
	void deleteSeg(int);
	std::vector<double> getAccel(int) const;
	std::vector<double> getAccelByIx(int) const;
	std::vector<tpat_constraint> getArcConstraints() const;
	std::vector<tpat_arc_piece> getChronoOrder() const;
	std::vector<double> getCoord(int) const;
	double getEpoch(int) const;
	double getEpochByIx(int) const;
	std::vector<double> getExtraParam(int, int) const;
	int getNextNodeID() const;
	int getNextSegID() const;
	tpat_node getNode(int) const;
	tpat_node getNodeByIx(int) const;
	int getNumCons() const;
	int getNumNodes() const;
	int getNumSegs() const;
	double getTOF(int) const;
	double getTOFByIx(int) const;
	virtual double getTotalTOF() const;
	tpat_segment getSeg(int) const;
	tpat_segment getSegByIx(int) const;
	std::vector<double> getState(int) const;
	std::vector<double> getStateByIx(int) const;
	MatrixXRd getSTM(int) const;
	MatrixXRd getSTMByIx(int) const;
	const tpat_sys_data* getSysData() const;
	double getTol() const;
	void putInChronoOrder();
	void setAccel(int, std::vector<double>);
	void setAccelByIx(int, std::vector<double>);
	void setState(int, std::vector<double>);
	void setStateByIx(int, std::vector<double>);
	void setSTMByIx(int, MatrixXRd);
	void setSTM(int, MatrixXRd);

	void setTol(double);

	// Utility Functions
	void printInChrono() const;

	/**
	 *  @brief Loads the object from a Matlab binary file
	 *  @param filepath the filepath to the Matlab file
	 */
	virtual void readFromMat(const char *filepath) = 0;

	/**
	 *	@brief Saves the object to a Matlab binary file
	 *	@param filepath the filepath to the Matlab file
	 */
	virtual void saveToMat(const char *filepath) const = 0;

	/**
	 *	@brief Displays a useful messages about the object
	 */
	virtual void print() const = 0;

protected:
	/** A pointer to the system data object that the describes the system this arc exists in */
	const tpat_sys_data *sysData;

	/** Contains all nodes or integration steps along an arc data object */
	std::vector<tpat_node> nodes {};

	/** Contains all segments that link the nodes of this object */
	std::vector<tpat_segment> segs {};

	/** Each entry corresponds to one node ID. The value of the entry is 
	 * the index of the node in the <tt>nodes</tt> array. If the value is
	 * equal to linkable::INVALID_ID, then the node no longer exists.
	 * 
	 * The current implementation requires that nextNodeID begin at 0  
	 * and increment by one through all integers.
	 */
	std::vector<int> nodeIDMap {};

	/**
	 * Each entry corresponds to one segment ID. The value of the entry is
	 * the index of the segment in the <tt>segs</tt> array. If the value is equal
	 * to linkable::INVALID_ID, then the segment no longer exists.
	 * 
	 * The current implementation requires that nextSegID begin at 0 and increment
	 * by one through all integers.
	 */
	std::vector<int> segIDMap {};

	/** A set of constraints that apply to the arc data object as a whole */
	std::vector<tpat_constraint> cons {};

	/** 
	 *	Number of variables stored in the extraParam vector. This
	 *	parameter should be set by the constructor of all derived
	 *	classes
	 */
	int numExtraParam = 0;

	/** 
	 *	Number of elements in each extra parameter. This parameter
	 *	should be set by the constructor of all derived classes
	 */
	std::vector<int> extraParamRowSize {};

	double tol = 0;		//!< Tolerance used to compute this data
	int nextNodeID = 0;	//!< A counter that stores the next available node ID
	int nextSegID = 0;	//!< A counter that stores the next available segment ID

	void copyMe(const tpat_base_arcset&);

	void initNodesSegsFromMat(mat_t *, const char*);
	void readStateFromMat(mat_t*, const char*);
	void readAccelFromMat(mat_t*);
	void readEpochFromMat(mat_t*, const char*);
	void readExtraParamFromMat(mat_t*, int, const char*);
	void readSTMFromMat(mat_t*);
	void readTOFFromMat(mat_t*, const char*);
	void saveAccel(mat_t*) const;
	void saveEpoch(mat_t*) const;
	void saveEpoch(mat_t*, const char*) const;
	void saveExtraParam(mat_t*, int, const char*) const;
	void saveState(mat_t*) const;
	void saveState(mat_t*, const char*) const;
	void saveSTMs(mat_t*) const;
	void saveTOF(mat_t*, const char*) const;
};//END OF tpat_base_arcset//--//--//--//--//--//--//--//--//--//--//--//--//

/**
 *  @brief A structure used to represent nodes and segments.
 *  @details This structure is used when putting the arcset
 *  object in chronological order
 */
struct tpat_arc_piece{

	/**
	 * @brief Enumerated type to describe the types of objects represented
	 * by the piece
	 * @details There are several types:
	 * 	* NODE - represent a tpat_node object
	 * 	* SEG - represent a tpat_segment object
	 */
	enum piece_tp {NODE, SEG};
	piece_tp type;	//!< The type of object represented by this piece
	int ID;	//!< The ID of the object represented by this piece

	/**
	 *  @brief Constructor
	 * 
	 *  @param tp the type of object represented by this piece
	 *  @param id the ID of the object represented by this piece
	 */
	tpat_arc_piece(piece_tp tp, int id) : type(tp), ID(id){};

	/**
	 *  @brief Comparison operator
	 * 
	 *  @param lhs piece reference
	 *  @param rhs piece reference
	 * 
	 *  @return true if the type and ID of the two pieces match
	 */
	friend bool operator ==(const tpat_arc_piece &lhs, const tpat_arc_piece &rhs){
		return lhs.type == rhs.type && lhs.ID == rhs.ID;
	}//================================================

	/**
	 *  @brief Comparison operator
	 * 
	 *  @param lhs piece reference
	 *  @param rhs piece reference
	 * 
	 *  @return true if the type and ID of the two pieces do NOT match
	 */
	friend bool operator !=(const tpat_arc_piece &lhs, const tpat_arc_piece &rhs){
		return !(lhs == rhs);
	}//================================================
};//END OF TPAT_ARC_PIECE//--//--//--//--//--//--//--//--//--//--//--//--//
#endif