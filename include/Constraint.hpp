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
 
#include <string>
#include <vector>

namespace astrohelion{

/**
 *	@brief Specify the type of constraint
 *
 *	This type tells the correction engine how to interpret and apply the 
 *	information stored in this constraint object.
 */
enum class Constraint_tp : int{
	NONE = 0, 					/*!< No constraint type specified. This type may be used as a
								 * placeholder, but attempting to complete a corrections process
								 * with a constraint of this type will result in a thrown exception
								 */
	STATE = 1,					/*!< Constrain specific core states to specified numeric values.
				 			 	 * 	The `id` value represents the node to constrain,
				 			 	 * 	and the `data` vector contains values for each state.
				 			 	 * 	To leave a state unconstrained, place a NAN value in the
				 			 	 * 	place of that state. 
				 			 	 */
	MATCH_ALL = 2, 				/*!< Make one node match another in all core states.
				 				 * 	The `id` value represents the source node, and the
				 				 * 	first value in `data` represents the node that will
				 				 * 	be matched to the first. Note that, for non-autonomous systems,
				 				 * 	the corrector will not attempt to match epoch time because it
				 				 * 	doesn't make sense to make two nodes occur at the same time.
								 */
	MATCH_CUST = 3,				/*!< Make one node match another in the specified core states.
				 				 *	The `id` value represents the source node; place the
				 				 *	index of the constrained node in the place of each state you
				 				 *	wish to constrain. For example, to constrain node 5 to match 
				 				 *	node 2 in x and z position, the data vector would contain
				 				 *	the following values: {5, NAN, 5, NAN, NAN, NAN} and the 
				 				 *	`id` value would be set to 2.
				 				 */
	EPOCH = 4,					/*!< Cosntrain the epoch of one node to a specific value.
					 			 *	The `id` value represents the node to constrain,
					 			 *	and the `data` value represents the fixed epoch time.
					 			 */
	DIST = 5,					/*!< Constrain a node to be at a specific distance from a primmary.
				 				 *	The `id` field identifies the constrained node, the 
				 				 *	first value in the `data` field (i.e. `data[0]`) 
				 				 *	identifies the primary number (index 0) and the second value in 
				 				 *	`data` gives the desired distance from the primary's center
				 				 *	in non-dimensional units
				 				 */
	MIN_DIST = 6, 				/*!< Constrain a node to be at or above a specific distance relative
				 				 *	to a specific primary. Data is input in the same way as the Constraint_tp::DIST
				 				 *	consraint described above.
				 				 */
	MAX_DIST = 7,  				/*!< Constrain a node to be at or below a specific distance relative to 
				 				 *	specific primary. Data is input in the same way as the Constraint_tp::DIST constraint
				 				 *	described above.
				 				 */
	DELTA_V = 8,				/*!< Constrain delta V to be exactly the specified amount
				 				 *	Data is input in the same format as Constraint_tp::MAX_DELTA_V
				 				 */
	MAX_DELTA_V = 9,			/*!< Constrain total delta V to be less than the specified amount
				 				 *	The `id` value is unused in this constraint; make the
				 				 *	first value in the data vector be the maximum delta-V in 
				 				 *	non-dimensional velocity units
				 				 */
	TOF_TOTAL = 10,				/*!< Constrain the trajectory to have a total time-of-flight
				 				 *	The index of the node is unused, and `data` holds the
				 				 *	value for the total TOF in non-dimensional units
				 				 */
 	APSE = 11,					/*!< Constrain a node to be located at an apse relative to one
				 				 *	of the primaries. The `data` field specified the index
				 				 *	of the primary (e.g. 0 for P1, 1 for P2, etc.) in the FIRST
				 				 *	element of the array. For example, if I want to constrain node
				 				 *	seven to be an apse relative to P1, I would set `id` 
				 				 *	to 7 and set `data` to [0, NAN, NAN, NAN, NAN, NAN]
				 				 */
	ANGLE = 12,					/*!< Constrain a node to be located on a plane normal to the XY-plane and 
								 * 	rotated by some angle about some point. Place the ID of the node you want 
								 * 	to constrain in the `id` variable. In the `data` variable, store the
								 * 	rotation base point (x0, y0, z0) and the rotation angle (rad) measured
								 * 	relative to the positive x-axis.
								 */
	CTRL = 13,					/*!< Constrain the control state at a node to have the specified values.
								 *	The `id` value represents the node to constrain,
				 			 	 * 	and the `data` vector contains values for each control state.
				 			 	 * 	To leave a state unconstrained, place a NAN value in the
				 			 	 * 	place of that state.
								 */
	JC = 100, 					/*!< Constrain the node to have a specific Jacobi constant (CR3BP only)
				 				 * 	Place the ID of the node you want to constrain in the `id`
				 				 *	variable; `data` holds the value of Jacobi
				 				 */
	SP = 200, 					/*!< Constrain the node to intersect the saddle point (BCR4BPR only)
				 				 * 	Place the ID of the node you want to constrain in the `id`
				 				 * 	variable; `data` is unused in this constraint
				 				 */
 	SP_RANGE = 201,				/*!< Constrain the node to be within a specified range of the saddle point
				 				 *	(BCR4BPR only). Place the ID of the node you want to constrain in 
				 				 * 	the `id` variable; The first element of `data` should
				 				 * 	hold the maximum allowable acceleration from the saddle point (zero) in non-dim units
				 				 */
 	SP_DIST = 202, 				/*!< Constrain the node to be at a specified non-dimensional distance
				 				 *	from the saddle point (BCR4BPR only). Place the ID of the node you wish
				 				 *	to constrain in the `id` variable; The first element of `data`
				 				 *	should contain the desired distance from the saddle point
				 				 *	in non-dim units, and the 2-10 elements should contain coefficients
				 				 * 	for a 2nd-order polynomial approximation of the saddle point's location
				 				 * 	as a function of epoch time. To generate these coefficients, run the
				 				 * 	function bcr4bpr_spLoc_polyFit and input the epoch of the node you wish 
				 				 * 	to constrain. The function returns the coefficients of three second-order
				 				 * 	functions in Epoch time: one for the x-position, one for y, and one for z.
				 				 * 	The 2-10 elements of `data` should be in the following order:<br>
				 				 * 	<br>2nd-order x coeff.<br>1st-order x coeff.<br>0th-order x coeff<br>
				 				 * 	2nd-order y coeff.<br>etc.
				 				 */
 	SP_MAX_DIST = 203,			/*!< Constrain the node to be at or below a specified non-dimensional
				 				 *	distance from the saddle point. All data fields should be entered in 
				 				 *	same way as for SP_DIST
				 				 */
	PSEUDOARC = 300, 			/*!< Pseudo arc-length continuation constraint. This forces a trajectory to
								 *	be in the same family as the input arc. The `id` attribute is 
								 *	meaningless for this constraint, but set it to the index of the last node.
								 * 	To ensure this is the final constraint, add it to the nodeset last. The `data`
								 *	field needs to contain the *entire* free-variable vector from a converged
								 *	family member. At the end of the data vector, append a double that represents
								 *	the continuation step size.
								 */
	CONT_CTRL = 500, 			/*!< Constrain a segment to be continuous with its terminal node in
								 * 	the specified control states. The `id` value specificies 
								 * 	the ID of the segment and the `data` field specifies 
								 * 	which states must be continuous. For example, to constrain the 
								 * 	first control state but not the second, the `data` field 
								 * 	would contain `[1 NAN]`. Values of NAN tell the algorithm 
								 * 	not to force continuity in that state.
								 */
 	CONT_PV = 501,				/*!< Constrain a segment to be continuous with its terminal node in 
								 * the specified position and velocity states. The
								 * `id` value specifies the ID of the segment, and the
								 * `data` field specifies which states must be continuous.
								 * For example, if I want the propagation along segment 4 to
								 * be continuous in all position states and x_dot, then I would set
								 * `id` to 4 and `data` to `[1 1 1 1 NAN NAN]`. Values
								 * of NAN tell the algorithm not to force continuity in that state,
								 * so in this example y_dot and z_dot are allowed to be discontinous.
								 * <b>These constraints are applied automatically by the corrections
								 * algorithm: DO NOT CREATE THESE.</b>
								 */
	CONT_EX = 502,				/*!< Constrain one of the extra parameters stored in a nodeset to be
								 * continuous along a segment. This may apply to epoch time, spacecraft mass,
								 * etc. Place the segment ID in `id` and the index of the extra
								 * parameter in the first data field, i.e. `data[0]`. The specified
								 * node will then be made continuous with the previous node in the nodeset
								 * provided it isn't the first one.
								 * <b>These constraints are applied automatically by the corrections
								 * algorithm: DO NOT CREATE THESE.</b>
								 * 
								 * **Parameter Indices**
								 * * 0 - Epoch
								 */
	SEG_CONT_PV = 503,			/*!< Constrain a segment to be continuous with another segment at their 
								 * terminal nodes in position and/or velocity. This situation may occur 
								 * if the corrections process shoots forward in time from one node and backward
								 * in time from another node; the two propagations may meet in the middle. Place 
								 * the ID of one of the segments in the `id` field, and the ID of the other 
								 * segment in the `data` field, one entry for each state that should be
								 * continuous. For example, to constrain the terminal points of segments
								 * 1 and 2 to be continuous in position, set `id` to 1 and set
								 * `data` to `[2, 2, 2, NAN, NAN, NAN]`
								 */
	SEG_CONT_EX = 504,			/*!< Constrain a segement to be continuous with another segment at their
								 * terminal nodes in some extra parameter, such as epoch or mass. Place
								 * the ID of one segment in `id`, the ID of the other segment in 
								 * `data[0]` and the parameter index in the elements after [0].
								 *	@see CONT_EX for a list of parameter indices
								 */
	RM_STATE = -1,				/*!< State constraint: Remove the indicated state vector from the free variable vector
								 * 	during multiple shooting. The `id` attribute specifies which node's state
								 * 	vector to remove. The `data` field is currently unused.
					     		 */
	RM_EPOCH = -4,				/*!< Epoch constraint: Remove the specified epoch from the free variable vector
								 *	during multiple shooting. The `id` attribute specifies which node's epoch
								 *	to remove. The `data` field is currently unused.
								 */
	RM_CTRL = -12, 				/*!< Control constraint: Remove the specified control state vector from the free variable
								 *	vector during multiple shooting. The `id` attribute specifies which node's 
								 *	control state vector to remove. The `data` is currently unused.
								 */
	ENDSEG_STATE = 700,			/*!< Constrain the end of a segment (not a node) to have a certain state
								 * 	The `id` value represents the segment to constrain,
				 			 	 * 	and the `data` vector contains values for each state.
				 			 	 * 	to leave a state unconstrained, place a NAN value in the
				 			 	 * 	place of that state. This constraint performs best when the constrained
				 			 	 * 	segment is the final segment on the trajectory.
								 */
	ENDSEG_APSE = 701,			/*!< Constrain the end of a segment (not a node) to be an apse
								 * 	relative to one of the primaries. The `data` field specifies 
								 * 	the index of the primary (e.g. 0 for P1, 1 for P2, etc.) in the FIRST
				 				 *	element of the array. For example, if I want to constrain segment
				 				 *	seven to be an apse relative to P1, I would set `id` 
				 				 *	to 7 and set `data` to [0]. This constraint performs best when 
				 				 *	the constrained segment is the final segment on the trajectory.
								 */
	ENDSEG_DIST = 702,			/*!< Constrain the end of a segment (not a node) to be a specified
								 * 	distance from the center of a specified primary. The `id` 
								 * 	field identifies the constrained segment, the 
				 				 *	first value in the `data` field (i.e. `data[0]`) 
				 				 *	identifies the primary number (index 0) and the second value in 
				 				 *	`data` gives the desired distance from the primary's center
				 				 *	in non-dimensional units. This constraint performs best when 
				 				 *	the constrained segment is the final segment on the trajectory.
								 */
	ENDSEG_MAX_DIST = 703,		/*!< Constrain the end of a segment (not a node) to be, at
								 * 	most, a specified distance (constraint the maximum 
								 * 	distance from the primary). Parameter storage is consistent
								 * 	with the details for Constraint_tp::ENDSEG_DIST. Note that
								 *	this constraint performs best when the constrained segment
								 *	is the final segment on the trajectory.
								 */
	ENDSEG_MIN_DIST = 704,		/*!< Constrain the end of a segment (not a node) to be at
								 * 	lease a specified distance (constraint the minimum 
								 * 	distance from the primary). Parameter storage is consistent
								 * 	with the details for Constraint_tp::ENDSEG_DIST. Note that
								 *	this constraint performs best when the constrained segment
								 *	is the final segment on the trajectory.
								 */
	ENDSEG_JC = 705,			/*!< Constrain the end of a segment (not a node) to have a specific 
								 * 	Jacobi constant (CR3BP only) Place the ID of the segment you want 
								 * 	to constrain in the `id` variable; `data` holds the value of Jacobi.
								 * 	Note that this constraint performs best when the constrained segment
								 * 	is the final segment on the trajectory.
				 				 */
	ENDSEG_ANGLE = 706			/*!< Constrain the segment end state to be located on a plane normal
							  	 * 	to the xy-plane and rotated about a specified point by a specified
								 *	angle. Place the ID of the segment you want to constrain in the `id`
								 *	variable; `data` stores the location of the base point and the rotation
								 *	angle: {x, y, z, angle}, where the angle is measured in radians.
								 *	Note that this constraint performs best when the constrained segment is
								 *	the final segment on the trajectory
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
 * 	* Update the setAppType() function
 *	* Add the new constraint type to the list of accepted constraints for
 * 		any dynamic models you wish
 *	* Define behavior for dealing with those types of constraints in the
 *		dynamic models
 *	* Add the constraint type to the model's multShoot_applyConstraint()
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
		/**
		 *  \name *structors
		 *  \{
		 */
		Constraint();
		Constraint(Constraint_tp);
		Constraint(Constraint_tp, int, std::vector<double>);
		Constraint(Constraint_tp, int, const double*, int);
		Constraint(const Constraint&);
		~Constraint();
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		Constraint& operator =(const Constraint&);
		friend bool operator ==(const Constraint&, const Constraint&);
		friend bool operator !=(const Constraint&, const Constraint&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		ConstraintApp_tp getAppType() const;
		static const char* getAppTypeStr(ConstraintApp_tp);
		static const char* getConTypeStr(Constraint_tp);
		std::vector<double> getData() const;
		int getID() const;
		double getFirstDataValue() const;
		double getFirstDataValue(int*) const;
		Constraint_tp getType() const;
		const char* getTypeStr() const;

		void setType(Constraint_tp);
		void setID(int);
		void setData(std::vector<double>);
		void setData(const double*, int);
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		int countConstrainedStates() const;
		//\}

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void print() const;
		//\}
	private:
		Constraint_tp type = Constraint_tp::NONE;		//!< The type of constraint
		ConstraintApp_tp appType = ConstraintApp_tp::APP_TO_NODE;	//!< How this constraint is applied
		int id = 0;						//!< object ID that this constraint applies to
		std::vector<double> data {};			//!< Data for this constraint

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void copyMe(const Constraint&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		void setAppType();
		//\}
};

}// END of Astrohelion namespace