/**
 *  @file BaseArcset.hpp
 *	@brief Basic arcset class (abstract)
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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

#pragma once

#include "Core.hpp"

#include "Constraint.hpp"
#include "Node.hpp"
#include "Segment.hpp"
#include "SysData.hpp"

#include "matio.h"
#include <map>
#include <memory>
#include <vector>

namespace astrohelion{

// Forward Declarations
struct ArcPiece;
class BaseArcset;

/**
 * @brief Smart pointer to a BaseArcset object
 */
typedef std::shared_ptr<BaseArcset> baseArcsetPtr;

/**
 *	@brief Abstract class that provides the framework for trajectories and nodesets
 *	
 *	The arcset object specifies default and mandatory behaviors for all derivative
 *	classes (i.e. Traj and Nodeset). All variables and data for an arc or 
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
class BaseArcset : public Core{

public:
	/**
	 *  @name *structors
	 *  @{
	 */
	BaseArcset(const SysData*);
	BaseArcset(const BaseArcset&);
	virtual baseArcsetPtr create( const SysData* ) const = 0;		//!< Virtual constructor for creation
	virtual baseArcsetPtr clone() const = 0;							//!< Virtual constructor for copying

	virtual ~BaseArcset();
	//@}

	/**
	 *  @name Operators
	 *  @{
	 */
	BaseArcset& operator =(const BaseArcset&);
	static void sum(const BaseArcset*, const BaseArcset*, BaseArcset*);
	//@}
	
	/**
	 *  @name Set and Get Functions
	 *  @{
	 */
	void addConstraint(Constraint);
	int addNode(Node);
	int addSeg(Segment);
	int appendSetAtNode(const BaseArcset*, int, int, double);
	void clearArcConstraints();
	void clearAllConstraints();
	std::vector<int> concatArcset(const BaseArcset*);
	void deleteNode(int);
	void deleteSeg(int);
	std::vector<double> getAccel(int);
	std::vector<double> getAccelByIx(int);
	std::vector<Constraint> getArcConstraints() const;
	std::vector<ArcPiece> getChronoOrder() const;
	std::vector<double> getCoord(unsigned int) const;
	double getEpoch(int) const;
	double getEpochByIx(int) const;
	std::vector<double> getEpochs() const;
	double getExtraParamByIx(int, std::string) const;
	std::vector<double> getExtraParamVecByIx(int, std::string) const;
	int getNextNodeID() const;
	int getNextSegID() const;
	Node getNode(int) const;
	Node getNodeByIx(int) const;
	Node& getNodeRef(int);
	Node& getNodeRefByIx(int);
	int getNodeIx(int) const;
	int getNumCons() const;
	int getNumNodes() const;
	int getNumSegs() const;
	double getTOF(int) const;
	double getTOFByIx(int) const;
	virtual double getTotalTOF() const;
	Segment getSeg(int) const;
	Segment getSegByIx(int) const;
	int getSegIx(int) const;
	std::vector<double> getState(int) const;
	std::vector<double> getStateByIx(int) const;
	MatrixXRd getSTM(int) const;
	MatrixXRd getSTMByIx(int) const;
	const SysData* getSysData() const;
	double getTol() const;
	bool isInChronoOrder() const;
	
	void putInChronoOrder(bool force = false);
	void setAccel(int, std::vector<double>);
	void setAccelByIx(int, std::vector<double>);
	void setState(int, std::vector<double>);
	void setStateByIx(int, std::vector<double>);
	void setSTMByIx(int, MatrixXRd);
	void setSTM(int, MatrixXRd);
	void setTol(double);
	void updateEpochs(int, double);
	//@}

	/**
	 *  @name Utility Functions
	 *  @{
	 */
	void printInChrono() const;
	void printNodeIDMap() const;
	void printSegIDMap() const;
	/**
	 *	@brief Displays a useful messages about the object
	 */
	virtual void print() const = 0;
	//@}

	/**
	 *  @name *File I/O
	 *  @{
	 */
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

	//@}
protected:
	/** A pointer to the system data object that the describes the system this arc exists in */
	const SysData *pSysData;

	/** Contains all nodes or integration steps along an arc data object */
	std::vector<Node> nodes {};

	/** Contains all segments that link the nodes of this object */
	std::vector<Segment> segs {};
	
	/**
	 * The key is the ID of the node, the value is the index in 
	 * node storage vector.
	 * 
	 * The current implementation requires that nextNodeID begin at 0  
	 * and increment by one through all integers.
	 */
	std::map<int, int> nodeIDMap {};

	/**
	 * The key is the ID of the segment, the value is the index in 
	 * segment storage vector.
	 * 
	 * The current implementation requires that nextSegID begin at 0  
	 * and increment by one through all integers.
	 */
	std::map<int, int> segIDMap {};

	/** A set of constraints that apply to the arc data object as a whole */
	std::vector<Constraint> cons {};

	double tol = 0;		//!< Tolerance used to compute this data
	int nextNodeID = 0;	//!< A counter that stores the next available node ID
	int nextSegID = 0;	//!< A counter that stores the next available segment ID

	bool bInChronoOrder = false;	//!< Whether or not the arcset is in chronological order

	void copyMe(const BaseArcset&);

	/**
	 *  @name File I/O
	 *  @{
	 */
	void initNodesSegsFromMat(mat_t *, const char*);
	void readStateFromMat(mat_t*, const char*);
	void readAccelFromMat(mat_t*);
	void readEpochFromMat(mat_t*, const char*);
	void readExtraParamFromMat(mat_t*, std::string, const char*);
	void readExtraParamVecFromMat(mat_t*, std::string, size_t, const char*);
	void readSTMFromMat(mat_t*);
	void readTOFFromMat(mat_t*, const char*);
	void saveAccel(mat_t*) const;
	void saveEpoch(mat_t*) const;
	void saveEpoch(mat_t*, const char*) const;
	void saveExtraParam(mat_t*, std::string, const char*) const;
	void saveExtraParamVec(mat_t*, std::string, size_t len, const char*) const;
	void saveState(mat_t*, const char* pVarName = "State") const;
	void saveSTMs(mat_t*) const;
	void saveTOF(mat_t*, const char*) const;
	//@}

	std::vector<ArcPiece> sortArcset(int, std::vector<ArcPiece>) const;
};//END OF BaseArcset//--//--//--//--//--//--//--//--//--//--//--//--//

/**
 *  @brief A structure used to represent nodes and segments.
 *  @details This structure is used when putting the arcset
 *  object in chronological order
 */
struct ArcPiece{

	/**
	 * @brief Enumerated type to describe the types of objects represented
	 * by the piece
	 */
	enum class Piece_tp {
		NODE,	//!< Represent a Node object
		SEG 	//!< Represent a Segment object
	};

	Piece_tp type;	//!< The type of object represented by this piece
	int id;	//!< The ID of the object represented by this piece

	/**
	 *  @brief Constructor
	 * 
	 *  @param tp the type of object represented by this piece
	 *  @param i the ID of the object represented by this piece
	 */
	ArcPiece(Piece_tp tp, int i) : type(tp), id(i) {}

	/**
	 *  @brief Comparison operator
	 * 
	 *  @param lhs piece reference
	 *  @param rhs piece reference
	 * 
	 *  @return true if the type and ID of the two pieces match
	 */
	friend bool operator ==(const ArcPiece &lhs, const ArcPiece &rhs){
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
	friend bool operator !=(const ArcPiece &lhs, const ArcPiece &rhs){
		return !(lhs == rhs);
	}//================================================
};//END OF ArcPiece//--//--//--//--//--//--//--//--//--//--//--//--//

}// END of Astrohelion namespace