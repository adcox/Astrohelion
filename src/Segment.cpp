/**
 *  @file Segment.cpp
 *	@brief Stores information about how nodes are linked
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

#include "Segment.hpp"

#include "Exceptions.hpp"
 

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Default constructor
 */
Segment::Segment(){}

/**
 *  @brief Construct a new segment
 * 	@details This constructor creates an STM equal to the 6x6 Identity matrix
 *  @param originID ID associated with the origin node
 *  @param terminusID ID associated with the terminal node (if no terminal node exists, 
 *  enter a value of NAN)
 *  @param tof time-of-flight along the segment
 */
Segment::Segment(int originID, int terminusID, double tof){
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
Segment::Segment(int originID, int terminusID, double tof, const double stm_data[36]){
	addLink(originID);
	addLink(terminusID);
	this->tof = tof;
	std::copy(stm_data, stm_data+36, stm);
}//====================================================

/**
 *  @brief Copy constructor
 *  @param s segment reference
 */
Segment::Segment(const Segment &s) : Linkable(s){
	copyMe(s);
}//====================================================

/**
 *  @brief Destructor
 */
// Segment::~Segment(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator
 *	@param s an arc segment object
 *	@return set this object equal to s and return *this
 */
Segment& Segment::operator =(const Segment &s){
	Linkable::operator =(s);
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
bool operator ==(const Segment &lhs, const Segment &rhs){
	
	if(lhs.tof != rhs.tof)
		return false;

	const Linkable link_lhs(lhs);
	const Linkable link_rhs(rhs);
	return link_lhs == link_rhs;
}//====================================================

/**
 *	@brief Determine if two segments are different
 *	@return whether two segments are different
 *	@see operator==
 */
bool operator != (const Segment &lhs, const Segment &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Add a constraint to the current set for this segment
 *	@param c a new constraint
 */
void Segment::addConstraint(Constraint c){
	cons.push_back(c);
}//====================================================

/**
 *	@brief Clear all constraints associated with this segment
 */
void Segment::clearConstraints(){ cons.clear(); }

/**
 *	@brief Get all constraints for this segment
 *	@return a vector containing all constraints applied to this segment
 */
std::vector<Constraint> Segment::getConstraints() const{
	return cons;
}//====================================================

/**
 *  @brief Retrieve the number of constraints stored by this object
 *  @return the number of constraints stored by this object
 */
int Segment::getNumCons() const { return (int)(cons.size()); }

/**
 *  @brief Retrieve the ID of the origin node (chronologically)
 *  @return the ID of the origin node (chronologically)
 */
int Segment::getOrigin() const { return links[ORIG_IX]; }

/**
 *  @brief Retrieve the ID of the terminal node (chronologically)
 *  @return the ID of the terminal node (chronologically)
 */
int Segment::getTerminus() const { return links[TERM_IX]; }

/**
 *  @brief Retrieve the time-of-flight along this segment
 *  @return time-of-flight along this segment, units consistent with the parent system
 */
double Segment::getTOF() const { return tof; }

/**
 *  @brief Retrieve the STM associated with this segment
 *  @return the STM associated with this segment
 */
MatrixXRd Segment::getSTM() const{
	double stmData[36];
	std::copy(stm, stm+36, stmData);
	MatrixXRd temp = Eigen::Map<MatrixXRd>(stmData, 6, 6);
	return temp;
}//====================================================

/**
 *  @brief Retrieve the STM elements in row-major order
 *  @return the STM elements in row-major order
 */
std::vector<double> Segment::getSTMElements() const{
	std::vector<double> stmVec;
	stmVec.insert(stmVec.begin(), stm, stm+36);
	return stmVec;
}//====================================================

/**
 *	@brief Retrieve a vector describing which of the velocity states
 *	should be made continuous with the node before this one in a nodeset.
 *
 *	@return a 3-element boolean vector. The first element describes
 *	whether or not the x-velocity is continuous, the second describes
 *	the y-velocity continuity, etc.
 */
std::vector<bool> Segment::getVelCon() const {
	std::vector<bool> temp;
	temp.insert(temp.end(), flags.begin(), flags.begin()+3);
	return temp;
}//====================================================

/**
 *	@brief Remove the specified constraint
 *	@param ix the index of the constraint. If the ix < 0, it will
 *	count backwards from the end of the set
 */
void Segment::removeConstraint(int ix){
	if(ix < 0)
		ix += cons.size();
	
	cons.erase(cons.begin() + ix);
}//====================================================

/**
 *	@brief Set the list of constraints for this arc segment
 *	@param constraints a vector of constraints
 */
void Segment::setConstraints(std::vector<Constraint> constraints){
	cons = constraints;
}//====================================================

/**
 *  @brief Set the ID of the node at the origin (chronologically) of this segment
 *  @param o the ID of the node at the origin (chronologically) of this segment
 */
void Segment::setOrigin(int o){ links[ORIG_IX] = o; }

/**
 *	@brief Set the STM for this step
 *	@param m a 6x6 state transition matrix (non-dim)
 *	@throw Exception if <tt>m</tt> is not 6x6
 */
void Segment::setSTM(MatrixXRd m){
	if(m.rows() != 6 || m.cols() != 6)
		throw Exception("Segment::setSTM: STM must be 6x6");
	
	std::copy(m.data(), m.data()+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements an array of stm elements in 
 *	row-major order. Note that the array must have at least 36
 *	elements, or un-initialized memory may be accessed
 */
void Segment::setSTM(const double *elements){
	std::copy(elements, elements+36, stm);
}//====================================================

/**
 *	@brief Set the STM for this step
 *	@param elements a 36-element vector of STM elements
 *	(non-dimensional) in row-major order
 *	@throw Exception if <tt>elements</tt> does not have 36 elements
 */
void Segment::setSTM(std::vector<double> elements){
	if(elements.size() != 36)
		throw Exception("Segment::setSTM: input vector must have 36 elements");
	std::copy(elements.begin(), elements.begin()+36, stm);
}//====================================================

/**
 *  @brief Set the ID of the node at the terminus (chronologically) of this segment
 *  @param t the ID of the node at the terminus (chronologically) of this segment
 */
void Segment::setTerminus(int t){ links[TERM_IX] = t; }

/**
 *  @brief Set the TOF along this segment
 *  @param t the TOF along this segment, units consistent with the parent system
 */
void Segment::setTOF(double t){ tof = t; }

/**
 *	@brief Set all velocity states to be continuous
 */
void Segment::setVel_AllCon(){
	flags[0] = true;
	flags[1] = true;
	flags[2] = true;
}//====================================================

/**
 *	@brief Set all velocity states to be discontinuous
 */
void Segment::setVel_AllDiscon(){
	flags[0] = false;
	flags[1] = false;
	flags[2] = false;
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param data a three-element boolean array. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 */
void Segment::setVelCon(const bool data[3]){
	flags[0] = data[0];
	flags[1] = data[1];
	flags[2] = data[2];
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param data a three-element boolean vector. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 *	@throw Exception if <tt>data</tt> has fewer than three elements
 */
void Segment::setVelCon(std::vector<bool> data){
	if(data.size() < 3)
		throw Exception("Segment::setVelCon: Need at least three velocity continuity booleans");

	std::copy(data.begin(), data.begin()+3, flags.begin());
}//====================================================

/**
 *	@brief Set the continuity for all three velocity states
 *	@param xCon whether or not the x-velocity component should 
 *	be continuous with the node before this one
 *	@param yCon whether or not the y-velocity component should 
 *	be continuous with the node before this one
 *	@param zCon whether or not the z-velocity component should 
 *	be continuous with the node before this one
 */
void Segment::setVelCon(bool xCon, bool yCon, bool zCon){
	flags[0] = xCon;
	flags[1] = yCon;
	flags[2] = zCon;
}//====================================================
//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

 /**
 *	@brief Copy a segment into this one
 *	@param s a segment reference
 */
void Segment::copyMe(const Segment &s){
	std::copy(s.stm, s.stm+36, stm);
	cons = s.cons;
	tof = s.tof;
	flags = s.flags;
	Linkable::copyMe(s);
}//====================================================


}// END of Astrohelion namespace