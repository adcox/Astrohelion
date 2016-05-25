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
#include <memory>
#include <vector>

// Forward Declarations
struct TPAT_Arc_Piece;
class TPAT_Base_Arcset;

/**
 * @brief Smart pointer to a TPAT_Base_Arcset object
 */
typedef std::shared_ptr<TPAT_Base_Arcset> baseArcsetPtr;

/**
 *	@brief Abstract class that provides the framework for trajectories and nodesets
 *	
 *	The arcset object specifies default and mandatory behaviors for all derivative
 *	classes (i.e. TPAT_Traj and TPAT_Nodeset). All variables and data for an arc or 
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
class TPAT_Base_Arcset : public TPAT{

public:
	// *structors
	TPAT_Base_Arcset(const TPAT_Sys_Data*);
	TPAT_Base_Arcset(const TPAT_Base_Arcset&);
	virtual baseArcsetPtr create( const TPAT_Sys_Data* ) const = 0;		//!< Virtual constructor for creation
	virtual baseArcsetPtr clone() const = 0;							//!< Virtual constructor for copying

	virtual ~TPAT_Base_Arcset();

	// Operators
	TPAT_Base_Arcset& operator =(const TPAT_Base_Arcset&);
	static void sum(const TPAT_Base_Arcset*, const TPAT_Base_Arcset*, TPAT_Base_Arcset*);

	// Set and Get functions
	void addConstraint(TPAT_Constraint);
	int addNode(TPAT_Node);
	int addSeg(TPAT_Segment);
	int appendSetAtNode(const TPAT_Base_Arcset*, int, int, double);
	void clearArcConstraints();
	void clearAllConstraints();
	std::vector<int> concatArcset(const TPAT_Base_Arcset*);
	void deleteNode(int);
	void deleteSeg(int);
	std::vector<double> getAccel(int) const;
	std::vector<double> getAccelByIx(int) const;
	std::vector<TPAT_Constraint> getArcConstraints() const;
	std::vector<TPAT_Arc_Piece> getChronoOrder() const;
	std::vector<double> getCoord(int) const;
	double getEpoch(int) const;
	double getEpochByIx(int) const;
	std::vector<double> getExtraParam(int, int) const;
	int getNextNodeID() const;
	int getNextSegID() const;
	TPAT_Node getNode(int) const;
	TPAT_Node getNodeByIx(int) const;
	int getNodeIx(int) const;
	int getNumCons() const;
	int getNumNodes() const;
	int getNumSegs() const;
	double getTOF(int) const;
	double getTOFByIx(int) const;
	virtual double getTotalTOF() const;
	TPAT_Segment getSeg(int) const;
	TPAT_Segment getSegByIx(int) const;
	int getSegIx(int) const;
	std::vector<double> getState(int) const;
	std::vector<double> getStateByIx(int) const;
	MatrixXRd getSTM(int) const;
	MatrixXRd getSTMByIx(int) const;
	const TPAT_Sys_Data* getSysData() const;
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
	void printNodeIDMap() const;
	void printSegIDMap() const;
	
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
	const TPAT_Sys_Data *sysData;

	/** Contains all nodes or integration steps along an arc data object */
	std::vector<TPAT_Node> nodes {};

	/** Contains all segments that link the nodes of this object */
	std::vector<TPAT_Segment> segs {};

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
	std::vector<TPAT_Constraint> cons {};

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

	void copyMe(const TPAT_Base_Arcset&);

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

	std::vector<TPAT_Arc_Piece> sortArcset(int, std::vector<TPAT_Arc_Piece>) const;
};//END OF TPAT_Base_Arcset//--//--//--//--//--//--//--//--//--//--//--//--//

/**
 *  @brief A structure used to represent nodes and segments.
 *  @details This structure is used when putting the arcset
 *  object in chronological order
 */
struct TPAT_Arc_Piece{

	/**
	 * @brief Enumerated type to describe the types of objects represented
	 * by the piece
	 * @details There are several types:
	 * 	* NODE - represent a TPAT_Node object
	 * 	* SEG - represent a TPAT_Segment object
	 */
	enum class Piece_Tp {NODE, SEG};
	Piece_Tp type;	//!< The type of object represented by this piece
	int id;	//!< The ID of the object represented by this piece

	/**
	 *  @brief Constructor
	 * 
	 *  @param tp the type of object represented by this piece
	 *  @param i the ID of the object represented by this piece
	 */
	TPAT_Arc_Piece(Piece_Tp tp, int i) : type(tp), id(i){};

	/**
	 *  @brief Comparison operator
	 * 
	 *  @param lhs piece reference
	 *  @param rhs piece reference
	 * 
	 *  @return true if the type and ID of the two pieces match
	 */
	friend bool operator ==(const TPAT_Arc_Piece &lhs, const TPAT_Arc_Piece &rhs){
		return lhs.type == rhs.type && lhs.id == rhs.id;
	}//================================================

	/**
	 *  @brief Comparison operator
	 * 
	 *  @param lhs piece reference
	 *  @param rhs piece reference
	 * 
	 *  @return true if the type and ID of the two pieces do NOT match
	 */
	friend bool operator !=(const TPAT_Arc_Piece &lhs, const TPAT_Arc_Piece &rhs){
		return !(lhs == rhs);
	}//================================================
};//END OF TPAT_Arc_Piece//--//--//--//--//--//--//--//--//--//--//--//--//
#endif