/**
 *  @file tpat_segment.cpp
 *	@brief Stores information about how nodes are linked
 *
 *	@author Andrew Cox
 *	@version April 28, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_segment.hpp"

#include "tpat_exceptions.hpp"
 
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Default constructor
 */
tpat_segment::tpat_segment(){ initArrays(); }

/**
 *  @brief Construct a new segment
 * 	@details This constructor creates an STM equal to the 6x6 Identity matrix
 *  @param originID ID associated with the origin node
 *  @param terminusID ID associated with the terminal node (if no terminal node exists, 
 *  enter a value of NAN)
 *  @param tof time-of-flight along the segment
 */
tpat_segment::tpat_segment(int originID, int terminusID, double tof){
	initArrays();
	addLink(originID);
	addLink(terminusID);
	this->tof = tof;
}//====================================================

/**
 *  @brief Construct a new segment
 *
 *  @param originID ID associated with the origin node
 *  @param terminusID ID associated with the terminal node (if no terminal node exists, 
 *  enter a value of NAN)
 *  @param tof time-of-flight along the segment
 *  @param stm_data a 36-element row-major order array containing a 6x6 state transition matrix
 *  associated with this segment
 */
tpat_segment::tpat_segment(int originID, int terminusID, double tof, const double stm_data[36]){
	initArrays();
	addLink(originID);
	addLink(terminusID);
	this->tof = tof;
	std::copy(stm_data, stm_data+36, stm);
}//====================================================

/**
 *  @brief Copy constructor
 *  @param s segment reference
 */
tpat_segment::tpat_segment(const tpat_segment &s) : tpat_linkable(s){
	copyMe(s);
}//====================================================

/**
 *  @brief Destructor
 *  @details Clears all data vectors
 */
tpat_segment::~tpat_segment(){
	cons.clear();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param s an arc segment object
 *	@return set this object equal to s and return *this
 */
tpat_segment& tpat_segment::operator =(const tpat_segment &s){
	tpat_linkable::operator =(s);
	copyMe(s);
	return *this;
}//====================================================

/**
 *	@brief Determine if two segments are identical
 *
 *	Conditions for identicalness:
 *	* Same origin, terminus, ID, and time-of-flight
 *	* (Not Active) Exact same constraint vector
 *
 *	@return whether or not two segments are identical
 */
bool operator ==(const tpat_segment &lhs, const tpat_segment &rhs){
	
	if(lhs.tof != rhs.tof)
		return false;

	const tpat_linkable link_lhs(lhs);
	const tpat_linkable link_rhs(rhs);
	return link_lhs == link_rhs;
}//====================================================

/**
 *	@brief Determine if two segments are different
 *	@return whether two segments are different
 *	@see operator==
 */
bool operator != (const tpat_segment &lhs, const tpat_segment &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a constraint to the current set for this segment
 *	@param c a new constraint
 */
void tpat_segment::addConstraint(tpat_constraint c){
	cons.push_back(c);
}//====================================================

/**
 *	@brief Clear all constraints associated with this segment
 */
void tpat_segment::clearConstraints(){ cons.clear(); }

/**
 *	@brief Get all constraints for this segment
 *	@return a vector containing all constraints applied to this segment
 */
std::vector<tpat_constraint> tpat_segment::getConstraints() const{
	return cons;
}//====================================================

/**
 *  @brief Retrieve the number of constraints stored by this object
 *  @return the number of constraints stored by this object
 */
int tpat_segment::getNumCons() const { return (int)(cons.size()); }

/**
 *  @brief Retrieve the ID of the origin node (chronologically)
 *  @return the ID of the origin node (chronologically)
 */
int tpat_segment::getOrigin() const { return links[ORIG_IX]; }

/**
 *  @brief Retrieve the ID of the terminal node (chronologically)
 *  @return the ID of the terminal node (chronologically)
 */
int tpat_segment::getTerminus() const { return links[TERM_IX]; }

/**
 *  @brief Retrieve the time-of-flight along this segment
 *  @return time-of-flight along this segment, units consistent with the parent system
 */
double tpat_segment::getTOF() const { return tof; }

/**
 *  @brief Retrieve the STM associated with this segment
 *  @return the STM associated with this segment
 */
MatrixXRd tpat_segment::getSTM() const{
	double stmData[36];
	std::copy(stm, stm+36, stmData);
	MatrixXRd temp = Eigen::Map<MatrixXRd>(stmData, 6, 6);
	return temp;
}//====================================================

/**
 *  @brief Retrieve the STM elements in row-major order
 *  @return the STM elements in row-major order
 */
std::vector<double> tpat_segment::getSTMElements() const{
	std::vector<double> stmVec;
	stmVec.insert(stmVec.begin(), stm, stm+36);
	return stmVec;
}//====================================================

/**
 *	@brief Remove the specified constraint
 *	@param ix the index of the constraint. If the ix < 0, it will
 *	count backwards from the end of the set
 */
void tpat_segment::removeConstraint(int ix){
	if(ix < 0)
		ix += cons.size();
	
	cons.erase(cons.begin() + ix);
}//====================================================

/**
 *	@brief Set the list of constraints for this arc segment
 *	@param constraints a vector of constraints
 */
void tpat_segment::setConstraints(std::vector<tpat_constraint> constraints){
	cons = constraints;
}//====================================================

/**
 *  @brief Set the ID of the node at the origin (chronologically) of this segment
 *  @param o the ID of the node at the origin (chronologically) of this segment
 */
void tpat_segment::setOrigin(int o){ links[ORIG_IX] = o; }

/**
 *	@brief Set the STM for this step
 *	@param m a 6x6 state transition matrix (non-dim)
 *	@throw tpat_exception if <tt>m</tt> is not 6x6
 */
void tpat_segment::setSTM(MatrixXRd m){
	if(m.rows() != 6 || m.cols() != 6)
		throw tpat_exception("tpat_segment::setSTM: STM must be 6x6");
	
	std::copy(m.data(), m.data()+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements an array of stm elements in 
 *	row-major order. Note that the array must have at least 36
 *	elements, or un-initialized memory may be accessed
 */
void tpat_segment::setSTM(const double *elements){
	std::copy(elements, elements+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements a 36-element vector of STM elements
 *	(non-dimensional) in row-major order
 *	@throw tpat_exception if <tt>elements</tt> does not have 36 elements
 */
void tpat_segment::setSTM(std::vector<double> elements){
	if(elements.size() != 36)
		throw tpat_exception("tpat_segment::setSTM: input vector must have 36 elements");
	std::copy(elements.begin(), elements.begin()+36, stm);
}//====================================================

/**
 *  @brief Set the ID of the node at the terminus (chronologically) of this segment
 *  @param t the ID of the node at the terminus (chronologically) of this segment
 */
void tpat_segment::setTerminus(int t){ links[TERM_IX] = t; }

/**
 *  @brief Set the TOF along this segment
 *  @param t the TOF along this segment, units consistent with the parent system
 */
void tpat_segment::setTOF(double t){ tof = t; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

 /**
 *	@brief Copy a segment into this one
 *	@param s a segment reference
 */
void tpat_segment::copyMe(const tpat_segment &s){
	// initArrays();
	std::copy(s.stm, s.stm+36, stm);
	cons = s.cons;
	tof = s.tof;
	tpat_linkable::copyMe(s);
}//====================================================

/**
 *	@brief Initialize all storage arrays
 */
void tpat_segment::initArrays(){
	for(int i = 0; i < 36; i++)
		stm[i] = i % 7 == 0 ? 1 : 0;
}//====================================================