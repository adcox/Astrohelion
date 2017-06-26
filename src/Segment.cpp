/**
 *  \file Segment.cpp
 *	\brief Stores information about how nodes are linked
 *
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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
 *  \brief Default constructor
 */
Segment::Segment(){}

/**
 *  \brief Construct a new segment
 * 	\details This constructor creates an STM equal to the 6x6 Identity matrix
 *  \param originID ID associated with the origin node
 *  \param terminusID ID associated with the terminal node (if no terminal node exists, 
 *  enter a value of NAN)
 *  \param tof time-of-flight along the segment
 */
Segment::Segment(int originID, int terminusID, double tof){
	addLink(originID);
	addLink(terminusID);
	this->tof = tof;
}//====================================================

/**
 *  \brief Construct a new segment
 *
 *  \param originID ID associated with the origin node
 *  \param terminusID ID associated with the terminal node (if no terminal node exists, 
 *  enter a value of NAN)
 *  \param tof time-of-flight along the segment
 *  \param stm_data a row-major order array containing a state transition matrix
 *  associated with this segment
 *  \param len the number of elements in stm_data (should be a perfect square)
 *  \param pLaw pointer to the control law implemented along this segment
 */
Segment::Segment(int originID, int terminusID, double tof, const double* stm_data, unsigned int len,
	ControlLaw* pLaw){

	addLink(originID);
	addLink(terminusID);
	this->tof = tof;
	this->pCtrlLaw = pLaw;

	stm = Eigen::Map<const MatrixXRd>(stm_data, std::sqrt(len), std::sqrt(len));
}//====================================================

/**
 *  \brief Copy constructor
 *  \param s segment reference
 */
Segment::Segment(const Segment &s) : Linkable(s){
	copyMe(s);
}//====================================================

/**
 *  \brief Destructor
 */
// Segment::~Segment(){}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	\brief Assignment operator
 *	\param s an arc segment object
 *	\return set this object equal to s and return *this
 */
Segment& Segment::operator =(const Segment &s){
	copyMe(s);
	return *this;
}//====================================================

/**
 *	\brief Determine if two segments are identical
 *
 *	Conditions for identicalness:
 *	* Same origin, terminus, ID, and time-of-flight
 *	* (Not Active) Exact same constraint vector
 *
 *	\return whether or not two segments are identical
 */
bool operator ==(const Segment &lhs, const Segment &rhs){
	
	if(lhs.tof != rhs.tof)
		return false;

	const Linkable link_lhs(lhs);
	const Linkable link_rhs(rhs);
	return link_lhs == link_rhs;
}//====================================================

/**
 *	\brief Determine if two segments are different
 *	\return whether two segments are different
 *	@see operator==
 */
bool operator != (const Segment &lhs, const Segment &rhs){
	return !(lhs == rhs);
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	\brief Add a constraint to the current set for this segment
 *	\param c a new constraint
 */
void Segment::addConstraint(Constraint c){
	cons.push_back(c);
}//====================================================

/**
 *  \brief Append a state vector to the storage vector
 *  \details The state does not have to be the full state;
 *  it is frequently useful to call this function multiple times
 *  with different pieces of the state
 * 
 *  \param q an array of state variables
 *  \param len length of the input array
 */
void Segment::appendState(const double *q, unsigned int len){
	states.insert(states.end(), q, q+len);
}//====================================================

/**
 *  \brief Append a state vector to the storage vector
 *  \details The state does not have to be the full state;
 *  it is frequently useful to call this function multiple times
 *  with different pieces of the state
 * 
 *  \param q a vector of state variables
 */
void Segment::appendState(const std::vector<double> q){
	states.insert(states.end(), q.begin(), q.end());
}//====================================================

/**
 *  \brief Append a time value to the storage vector
 *  \param t a time value (nondimensional)
 */
void Segment::appendTime(double t){
	times.push_back(t);
}//====================================================

/**
 *	\brief Clear all constraints associated with this segment
 */
void Segment::clearConstraints(){ cons.clear(); }

/**
 *	\brief Get all constraints for this segment
 *	\return a vector containing all constraints applied to this segment
 */
std::vector<Constraint> Segment::getConstraints() const{
	return cons;
}//====================================================

/**
 *  \brief Retrieve the control law ID associated with this segment
 *  \return the control law ID associated with this segment
 */
ControlLaw* Segment::getCtrlLaw() const{ return pCtrlLaw; }

/**
 *  \brief Retrieve the number of constraints stored by this object
 *  \return the number of constraints stored by this object
 */
unsigned int Segment::getNumCons() const { return cons.size(); }

/**
 *  \brief Retrieve the number of times stored in the time vector
 *  \return the number of times stored in the time vector
 */
unsigned int Segment::getNumTimes() const { return times.size(); }

/**
 *  \brief Retrieve the ID of the origin node (chronologically)
 *  \return the ID of the origin node (chronologically)
 */
int Segment::getOrigin() const { return links[ORIG_IX]; }

/**
 *  \brief Retrieve a copy of the entire state vector
 *  \details The state vector includes, in row-major order,
 *  every state along the Segment. Depending on how the propagation
 *  is conducted, there may be as few as two states. The width of
 *  one state (i.e., one row) is available from getStateWidth()
 *  
 *  \return a copy of the entire state vector
 */
std::vector<double> Segment::getStateVector() const{ return states; }

/**
 *  \brief Retrieve the size of one state vector in the array
 *  \details States are stored in row-major order in an array;
 *  the state width is the number of values in a single state, 
 *  or row, of the array
 * 
 *  \return the number of elements in a single state
 */
unsigned int Segment::getStateWidth() const{ return stateWidth; }

/**
 *  \brief Retrieve a state from the storage vector
 *  \details State vectors are stored as rows of a matrix where
 *  the storage vector represents the matrix in row-major order.
 *  The length of each row is set/retrieved via setStateWidth() and
 *  getStateWidth()
 * 
 *  \param row Row index, begins at 0. If row < 0, it will count backwards
 *  from the end of the set of states.
 * 
 *  \return A vector containing the desired state values
 */
std::vector<double> Segment::getStateByRow(int row) const{
	if(stateWidth == 0)
		throw Exception("Segment::getStateByRow: stateWidth is zero! Set the state width before calling this function.");

	if(row < 0)
		row += states.size()/stateWidth;

	if(row < 0 || row >= static_cast<int>(states.size()/stateWidth)){
		char msg[128];
		sprintf(msg, "Segment::getStateByRow: Index %d out of bounds; expected between 0 and %d", row, static_cast<int>(states.size()/stateWidth) - 1);
		throw Exception(msg);
	}

	std::vector<double> q(states.begin()+row*stateWidth, states.begin()+(row+1)*stateWidth);
	return q;
}//====================================================

/**
 *  \brief Retrieve the ID of the terminal node (chronologically)
 *  \return the ID of the terminal node (chronologically)
 */
int Segment::getTerminus() const { return links[TERM_IX]; }

/**
 *  \brief Retrieve a copy of the entire time vector.
 *  \details Each time value corresponds to one state stored
 *  in the state vector.
 *  
 *  \return a copy of the entire time vector
 */
std::vector<double> Segment::getTimeVector() const { return times; }

/**
 *  \brief Retrieve a time from the segment time vector
 *  \details [long description]
 * 
 *  \param ix index within the time vector; if ix is negative, 
 *  it will count backward from the end of the vector
 *  \return The time value at the specified index
 *  
 *  \throws Exception if the index is out of bounds
 */
double Segment::getTimeByIx(int ix) const{
	if(ix < 0)
		ix += times.size();

	if(ix < 0 || ix >= static_cast<int>(times.size())){
		char msg[128];
		sprintf(msg, "Segment::getTimeByIx: Index %d out of bounds; expected between 0 and %zu", ix, times.size() - 1);
		throw Exception(msg);
	}

	return times[ix];
}//====================================================

/**
 *  \brief Retrieve the time-of-flight along this segment
 *  \return time-of-flight along this segment, units consistent with the parent system
 */
double Segment::getTOF() const { return tof; }

/**
 *  \brief Retrieve the STM associated with this segment
 *  \return the STM associated with this segment
 */
MatrixXRd Segment::getSTM() const{ return stm; }

/**
 *	\brief Retrieve a vector describing which of the velocity states
 *	should be made continuous with the segment before this one in a nodeset.
 *
 *	\return a 3-element boolean vector. The first element describes
 *	whether or not the x-velocity is continuous, the second describes
 *	the y-velocity continuity, etc.
 */
std::vector<bool> Segment::getVelCon() const {
	std::vector<bool> temp;
	temp.insert(temp.end(), flags.begin(), flags.begin()+3);
	return temp;
}//====================================================

/**
 *	\brief Remove the specified constraint
 *	\param ix the index of the constraint. If the ix < 0, it will
 *	count backwards from the end of the set
 */
void Segment::removeConstraint(int ix){
	if(ix < 0)
		ix += cons.size();
	
	cons.erase(cons.begin() + ix);
}//====================================================

/**
 *	\brief Set the list of constraints for this arc segment
 *	\param constraints a vector of constraints
 */
void Segment::setConstraints(std::vector<Constraint> constraints){
	cons = constraints;
}//====================================================

/**
 *  \brief Set the control law for this segment
 *  \param law a pointer to the control law applied on this segment
 */
void Segment::setCtrlLaw(ControlLaw *law){ pCtrlLaw = law; }

/**
 *  \brief Set the ID and also update the ID of any 
 *  associated constraints
 * 
 *  \param id ID that uniquely identifies the node
 */
void Segment::setID(int id){
	Linkable::setID(id);
	for(unsigned int i = 0; i < cons.size(); i++){
		cons[i].setID(id);
	}
}//====================================================

/**
 *  \brief Set the ID of the node at the origin (chronologically) of this segment
 *  \param o the ID of the node at the origin (chronologically) of this segment
 */
void Segment::setOrigin(int o){ links[ORIG_IX] = o; }

/**
 *  \brief Set the state vector associated with the propagation along this segment
 *  \param v the state vector associated with the propagation along this segment, 
 *  in row-major order.
 */
void Segment::setStateVector(std::vector<double> v){ states = v; }

/**
 *  \brief Specify the size of one state vector in the array
 *  \details States are stored in row-major order in an array;
 *  the state width is the number of values in a single state, 
 *  or row, of the array
 * 
 *  \param w the number of elements in a single state
 */
void Segment::setStateWidth(unsigned int w){ stateWidth = w; }

/**
 *	\brief Set the STM for this step
 *	\param m a state transition matrix (non-dim)
 */
void Segment::setSTM(MatrixXRd m){ stm = m; }

/**
 *	\brief Set the STM for this step
 *	\param elements an array of stm elements in 
 *	row-major order.
 *	\param len number of numbers in elements
 */
void Segment::setSTM(const double *elements, unsigned int len){
	stm = Eigen::Map<const MatrixXRd>(elements, std::sqrt(len), std::sqrt(len));
}//====================================================

/**
 *	\brief Set the STM for this step
 *	\param elements a 36-element vector of STM elements
 *	(non-dimensional) in row-major order
 */
void Segment::setSTM(std::vector<double> elements){
	unsigned int len = static_cast<unsigned int>(std::floor(std::sqrt(static_cast<double>(elements.size()))));
	stm = Eigen::Map<MatrixXRd>(&(elements[0]), len, len);
}//====================================================

/**
 *  \brief Set the ID of the node at the terminus (chronologically) of this segment
 *  \param t the ID of the node at the terminus (chronologically) of this segment
 */
void Segment::setTerminus(int t){ links[TERM_IX] = t; }

/**
 *  \brief Set the time vector associated with the propagated arc along this segment
 *  \param t he time vector associated with the propagated arc along this segment
 */
void Segment::setTimeVector(std::vector<double> t){ times = t; }

/**
 *  \brief Set the TOF along this segment; 
 *  \param t the TOF along this segment, units consistent with the parent system
 */
void Segment::setTOF(double t){ tof = t; }

/**
 *	\brief Set all velocity states to be continuous
 */
void Segment::setVel_AllCon(){
	flags[0] = true;
	flags[1] = true;
	flags[2] = true;
}//====================================================

/**
 *	\brief Set all velocity states to be discontinuous
 */
void Segment::setVel_AllDiscon(){
	flags[0] = false;
	flags[1] = false;
	flags[2] = false;
}//====================================================

/**
 *	\brief Set the continuity for all three velocity states
 *	between this segment and the previous segment
 *	\param data a three-element boolean array. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 */
void Segment::setVelCon(const bool data[3]){
	flags[0] = data[0];
	flags[1] = data[1];
	flags[2] = data[2];
}//====================================================

/**
 *	\brief Set the continuity for all three velocity states
 *	between this segment and the previous segment
 *	\param data a three-element boolean vector. Each element
 *	corresponds to one of the velocity states in the order
 *	[v_x, v_y, v_z]
 *	@throw Exception if <tt>data</tt> has fewer than three elements
 */
void Segment::setVelCon(std::vector<bool> data){
	if(data.size() < 3){
		char msg[128];
		sprintf(msg, "Segment::setVelCon: data size = %zu; need at least three elements", data.size());
		throw Exception(msg);
	}

	std::copy(data.begin(), data.begin()+3, flags.begin());
}//====================================================

/**
 *	\brief Set the continuity for all three velocity states
 *	between this segment and the previous segment
 *	\param xCon whether or not the x-velocity component should 
 *	be continuous with the segment before this one
 *	\param yCon whether or not the y-velocity component should 
 *	be continuous with the segment before this one
 *	\param zCon whether or not the z-velocity component should 
 *	be continuous with the segment before this one
 */
void Segment::setVelCon(bool xCon, bool yCon, bool zCon){
	flags[0] = xCon;
	flags[1] = yCon;
	flags[2] = zCon;
}//====================================================

/**
 *  \brief Shift all times in the time vector by the 
 *  specified amount
 * 
 *  \param deltaT amount of time to shift by
 */
void Segment::shiftAllTimes(double deltaT){
	for(unsigned int i = 0; i < times.size(); i++)
		times[i] += deltaT;
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

 /**
 *	\brief Copy a segment into this one
 *	\param s a segment reference
 */
void Segment::copyMe(const Segment &s){
	stm = s.stm;
	cons = s.cons;
	tof = s.tof;
	flags = s.flags;
	pCtrlLaw = s.pCtrlLaw;	// Copying ADDRESS
	times = s.times;
	states = s.states;
	stateWidth = s.stateWidth;
	Linkable::copyMe(s);
}//====================================================

/**
 *  \brief Compute the time-of-flight along a segment from the data stored
 *  in the times vector
 *  \details If there are fewer than 2 time values in the vector,
 *  the time of flight is set to zero
 */
void Segment::updateTOF(){
	if(times.size() > 1){
		tof = times.back() - times.front();
	}else{
		tof = 0;
	}
}//====================================================

/**
 *  \brief Print a description of the segment
 */
void Segment::print() const{
	printf("Segment | id = %d\n", ID);
	printf("\tOrigin Node ID: %d, Terminus Node ID: %d\n", getOrigin(), getTerminus());
	printf("\tTOF = %.4f\n", tof);
	printf("\tControl Law = %s\n", pCtrlLaw ? pCtrlLaw->lawTypeToString(pCtrlLaw->getLawType()).c_str() : "NONE");

	printf("\tTime Vector: %zu x 1\n", times.size());
	if(times.size() > 0){
		printf("\tState Vector: %zu x %d (remainder %f)\n", times.size(), static_cast<int>(floor(states.size()/times.size())),
			states.size() - times.size()*floor(states.size()/times.size()));
	}else{
		printf("\tState Vector: %zu x 0\n", states.size());
	}
	printf("Saved state width = %u\n", stateWidth);
}//====================================================



}// END of Astrohelion namespace