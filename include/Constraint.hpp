/**
 *  @file Constraint.hpp
 *	@brief Contains information about a constraint on a nodeset, node, or segment
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
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
 
#include <string>
#include <vector>

namespace astrohelion{
/**
 *	@brief Specify the type of constraint
 *
 *	This type tells the correction engine how to interpret and apply the 
 *	information stored in this constraint object.
 */
enum class Constraint_tp {
	NONE, 		/*!< No constraint type specified. This type may be used as a
				 * placeholder, but attempting to complete a corrections process
				 * with a constraint of this type will result in a thrown exception
				 */
	STATE,		/*!< Constrain specific states to specified numeric values.
 			 	 * The <tt>id</tt> value represents the node to constrain,
 			 	 * and the <tt>data</tt> vector contains values for each state.
 			 	 * to leave a state unconstrained, place a NAN value in the
 			 	 * place of that state. 
 			 	 */
	MATCH_ALL, 	/*!< Make one node match another in all states.
 				 * The <tt>id</tt> value represents the source node, and the
 				 * first value in <tt>data</tt> represents the node that will
 				 * be matched to the first. Note that, for non-autonomous systems,
 				 * the corrector will not attempt to match epoch time because it
 				 * doesn't make sense to make two nodes occur at the same time.
				 */
	MATCH_CUST,	/*!< Make one node match another in the specified states.
 				 *	The <tt>id</tt> value represents the source node; place the
 				 *	index of the constrained node in the place of each state you
 				 *	wish to constrain. For example, to constrain node 5 to match 
 				 *	node 2 in x and z position, the data vector would contain
 				 *	the following values: {5, NAN, 5, NAN, NAN, NAN} and the 
 				 *	<tt>id</tt> value would be set to 2.
 				 */
	DIST,		/*!< Constrain a node to be at a specific distance from a primmary.
 				 *	The <tt>id</tt> field identifies the constrained node, the 
 				 *	first value in the <tt>data</tt> field (i.e. <tt>data[0]</tt>) 
 				 *	identifies the primary number (index 0) and the second value in 
 				 *	<tt>data</tt> gives the desired distance from the primary's center
 				 *	in non-dimensional units
 				 */
	MIN_DIST, 	/*!< Constrain a node to be at or above a specific distance relative
 				 *	to a specific primary. Data is input in the same way as the Constraint_tp::DIST
 				 *	consraint described above.
 				 */
	MAX_DIST,  	/*!< Constrain a node to be at or below a specific distance relative to 
 				 *	specific primary. Data is input in the same way as the Constraint_tp::DIST constraint
 				 *	described above.
 				 */
	MAX_DELTA_V,/*!< Constrain total delta V to be less than the specified amount
 				 *	The <tt>id</tt> value is unused in this constraint; make the
 				 *	first value in the data vector be the maximum delta-V in 
 				 *	non-dimensional velocity units
 				 */
	DELTA_V,	/*!< Constrain delta V to be exactly the specified amount
 				 *	Data is input in the same format as Constraint_tp::MAX_DELTA_V
 				 */
	SP, 		/*!< Constrain the node to intersect the saddle point (BCR4BPR only)
 				 * 	Place the ID of the node you want to constrain in the <tt>id</tt>
 				 * 	variable; <tt>data</tt> is unused in this constraint
 				 */
 	SP_RANGE,	/*!< Constrain the node to be within a specified range of the saddle point
 				 *	(BCR4BPR only). Place the ID of the node you want to constrain in 
 				 * 	the <tt>id</tt> variable; The first element of <tt>data</tt> should
 				 * 	hold the maximum allowable acceleration from the saddle point (zero) in non-dim units
 				 */
 	SP_DIST, 	/*!< Constrain the node to be at a specified non-dimensional distance
 				 *	from the saddle point (BCR4BPR only). Place the ID of the node you wish
 				 *	to constrain in the <tt>id</tt> variable; The first element of <tt>data</tt>
 				 *	should contain the desired distance from the saddle point
 				 *	in non-dim units, and the 2-10 elements should contain coefficients
 				 * 	for a 2nd-order polynomial approximation of the saddle point's location
 				 * 	as a function of epoch time. To generate these coefficients, run the
 				 * 	function bcr4bpr_spLoc_polyFit and input the epoch of the node you wish 
 				 * 	to constrain. The function returns the coefficients of three second-order
 				 * 	functions in Epoch time: one for the x-position, one for y, and one for z.
 				 * 	The 2-10 elements of <tt>data</tt> should be in the following order:<br>
 				 * 	<br>2nd-order x coeff.<br>1st-order x coeff.<br>0th-order x coeff<br>
 				 * 	2nd-order y coeff.<br>etc.
 				 */
 	SP_MAX_DIST,/*!< Constrain the node to be at or below a specified non-dimensional
 				 *	distance from the saddle point. All data fields should be entered in 
 				 *	same way as for SP_DIST
 				 */
 	JC, 		/*!< Constrain the node to have a specific Jacobi constant (CR3BP only)
 				 * 	Place the ID of the node you want to cosntraint in the <tt>id</tt>
 				 *	variable; <tt>data</tt> holds the value of Jacobi
 				 */
 	TOF,		/*!< Constrain the trajectory to have a total time-of-flight
 				 *	The index of the node is unused, and <tt>data</tt> holds the
 				 *	value for the total TOF in non-dimensional units
 				 */
 	APSE,		/*!< Constrain a node to be located at an apse relative to one
 				 *	of the primaries. The <tt>data</tt> field specified the index
 				 *	of the primary (e.g. 0 for P1, 1 for P2, etc.) in the FIRST
 				 *	element of the array. For example, if I want to constrain node
 				 *	seven to be an apse relative to P1, I would set <tt>id</tt> 
 				 *	to 7 and set <tt>data</tt> to [0, NAN, NAN, NAN, NAN, NAN]
 				 */
 	CONT_PV,	/*!< Constrain a segment to be continuous with its terminal node in 
				 * the specified position and velocity states. The
				 * <tt>id</tt> value specifies the ID of the node, and the
				 * <tt>data</tt> field specifies which states must be continuous.
				 * For example, if I want the propagation along segment 4 to
				 * be continuous in all position states and x_dot, then I would set
				 * <tt>id</tt> to 4 and <tt>data</tt> to [1 1 1 1 NAN NAN]. Values
				 * of NAN tell the algorithm not to force continuity in that state,
				 * so in this example y_dot and z_dot are allowed to be discontinous.
				 * NOTE: These constraints are applied automatically by the corrections
				 * algorithm: DO NOT CREATE THESE.
				 */
	CONT_EX,	/*!< Constrain one of the extra parameters stored in a nodeset to be
				 * continuous along a segment. This may apply to epoch time, spacecraft mass,
				 * etc. Place the segment ID in <tt>id</tt> and the index of the extra
				 * parameter in the first data field, i.e. <tt>data[0]</tt>. The specified
				 * node will then be made continuous with the previous node in the nodeset
				 * provided it isn't the first one.
				 * NOTE: These constraints are applied automatically by the corrections
				 * algorithm: DO NOT CREATE THESE.
				 * 
				 * **Parameter Indices**
				 * * 0 - Epoch
				 * * 1 - Mass		[Not Implemented]
				 * * 2 - Attitude 	[Not Implemented]
				 */
	SEG_CONT_PV,/*!< Constrain a segment to be continuous with another segment at their 
				 * terminal nodes in position and/or velocity. This situation may occur 
				 * if the corrections process shoots forward in time from one node and backward
				 * in time from another node; the two propagations may meet in the middle. Place 
				 * the ID of one of the segments in the <tt>id</tt> field, and the ID of the other 
				 * segment in the <tt>data</tt> field, one entry for each state that should be
				 * continuous. For example, to constrain the terminal points of segments
				 * 1 and 2 to be continuous in position, set <tt>id</tt> to 1 and set
				 * <tt>data</tt> to <tt>[2, 2, 2, NAN, NAN, NAN]</tt>
				 */
	SEG_CONT_EX,/*!< Constrain a segement to be continuous with another segment at their
				 * terminal nodes in some extra parameter, such as epoch or mass. Place
				 * the ID of one segment in <tt>id</tt>, the ID of the other segment in 
				 * <tt>data[0]</tt> and the parameter index in the elements after [0].
				 *	@see CONT_EX for a list of parameter indices
				 */
	PSEUDOARC 	/*!< Pseudo arc-length continuation constraint. This forces a trajectory to
				 *	be in the same family as the input arc. The <tt>id</tt> attribute is 
				 *	meaningless for this constraint, but set it to the index of the last node.
				 * 	To ensure this is the final constraint, add it to the nodeset last. The <tt>data</tt>
				 *	field needs to contain the *entire* free-variable vector from a converged
				 *	family member. At the end of the data vector, append a double that represents
				 *	the continuation step size.
				 */
};

/**
 *  @brief Describes how a constraint is applied, i.e., what type of object is controlled
 *  by the constraint.
 */
enum class ConstraintApp_tp{
	APP_TO_NODE,	//!< Constraint applies to nodes
	APP_TO_SEG,		//!< Constraint applies to segments
	APP_TO_ARC		//!< Constraint applies to arcs (as a whole)
};
		
/**
 *	@brief Contains information about how a particular node, segment, or arcset
 *	should be constrained during a corrections process
 *
 * 	**Adding a New Constraint**
 * 	* Create a new enumerated type and document it fully
 * 	* Update the getTypStr() function
 * 	* Update the segAppType() function
 *	* Add the new constraint type to the list of accepted constraints for
 * 		any dynamic models you wish
 *	* Define behavior for dealing with those types of constraints in the
 *		dynamic models
 *	* Add the constraint type to the models multShoot_applyConstraint()
 *		function so it will be called by the corrector
 *	* Add the constraint type to the corrector initializer so it knows how
 *		many rows the constraint occupies
 *	* If the constraint uses slack variables, write a function to compute
 *		the default value of the slack variable (if possible) in the relevant
 *		model and update the multShoot_getSlackVarVal() function.
 *
 *	@author Andrew Cox
 *	@version August 3, 2015
 *	@copyright GNU GPL v3.0
 */
class Constraint : public Core{
	public:
		Constraint();
		Constraint(Constraint_tp);
		Constraint(Constraint_tp, int, std::vector<double>);
		Constraint(Constraint_tp, int, const double*, int);
		Constraint(const Constraint&);
		~Constraint();
		
		Constraint& operator =(const Constraint&);

		int countConstrainedStates() const;
		
		ConstraintApp_tp getAppType() const;
		static const char* getAppTypeStr(ConstraintApp_tp);
		static const char* getConTypeStr(Constraint_tp);
		std::vector<double> getData() const;
		int getID() const;
		double getFirstDataValue() const;
		double getFirstDataValue(int*) const;
		Constraint_tp getType() const;
		const char* getTypeStr() const;
		
		void print() const;

		void setType(Constraint_tp);
		void setID(int);
		void setData(std::vector<double>);
		void setData(const double*, int);
		
		
	private:
		void copyMe(const Constraint&);
		void setAppType();
		
		Constraint_tp type = Constraint_tp::NONE;		//!< The type of constraint
		ConstraintApp_tp appType = ConstraintApp_tp::APP_TO_NODE;	//!< How this constraint is applied
		int id = 0;						//!< object ID that this constraint applies to
		std::vector<double> data {};			//!< Data for this constraint
};

}// END of Astrohelion namespace