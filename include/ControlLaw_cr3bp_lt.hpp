/**
 * @file ControlLaw_cr3bp_lt.hpp
 * @brief Control Law for CR3BP-LT system header file 
 * 
 * @author Andrew Cox
 * @version March 3, 2017
 * @copyright GNU GPL v3.0
 */
 
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "ControlLaw.hpp"
#include "EigenDefs.hpp"

namespace astrohelion{

// Forward declarations
class Arcset_cr3bp_lt;
class SysData;
class SysData_cr3bp_lt;

/**
 *  \ingroup model cr3bp_lt
 *  @brief Includes control laws specific to the CR3BP-LT problem
 *  @details As specified in the MathSpec document, the CR3BP-LT equations of 
 *  motion are given by
 *  \f{align*}{
 *  	\vec{a} &= \mathbf{K_v} \vec{f} + \dpd{\Omega}{\vec{r}} + \vec{a}_{lt},
 *  	\qquad \vec{a}_lt = \frac{f}{m} \hat{u}\\
 *  	\dot{m} &= \dod{m}{\tau}\\
 *  	\dot{\vec{\gamma}} &= \dod{\vec{\gamma}}{\tau}
 *  \f}
 *  where \f$ \vec{a} = \{\ddot{x}\,\, \ddot{y}\,\, \ddot{z} \}\f$ is the 
 *  rotating acceleration vector,
 *  \f$ \vec{v} = \{\dot{x}\,\, \dot{y}\,\, \dot{z} \} \f$ is the rotating 
 *  velocity vector, and \f$ \vec{r} = \{x\,\, y\,\, z \} \f$ is the spacecraft 
 *  position vector in the rotating frame. The CR3BP pseudopotential is 
 *  represented by \f$ \Omega \f$, the nondimensional low-thrust magnitude if 
 *  \f$ f \f$, the nondimensional spacecraft mass is \f$ m \in [0,1] \f$, and 
 *  the thrust pointing unit vector is \f$ \hat{u} \f$. While most of the core 
 *  state eoms, \f$ \vec{a} \f$, are computed in the DynamicsModel_cr3bp_lt
 *  class, several of the control-specific elements are computed in the control 
 *  law object. The low-thrust acceleration, \f$ \vec{a}_{lt} \f$, is obtained 
 *  from the `getOutput()` function,
 *  
 *  	a_lt = getOutput(t, s, pSysData, partials, len);
 *  	
 *  where the `s` array contains a full state, which is generally of the form
 *  \f$ \vec{s} = \begin{Bmatrix} \vec{q} & \vec{\gamma} & \phi_1 & \dots & 
 *  \phi_n & \vec{q}_{\text{extra}} \end{Bmatrix} \f$.
 *  In this vector, the core and control state variables are
 *  \f[
 *  	\vec{q} = \begin{Bmatrix} \vec{r} \,\, \vec{v} \,\, m \end{Bmatrix}^T,
 *  	\qquad \vec{\gamma} \in \mathbb{R}^{n\times 1},
 *  \f]
 *  respectively, \f$ \phi_i \f$ are the state transition matrix elements in 
 *  row-major order, and \f$ \vec{q}_{\text{extra}} \f$ are the "extra" state 
 *  variables included in some models (not in this model, however).
 *  
 *	Because the mass time derivative is a function of thrust and Isp, both 
 *	control law parameters (or states), the second equation of motion is 
 *	computed via the control law function `get_dmdt()`. The final EOM, the 
 *	control state time derivative, is computed via the control law function 
 *	`getTimeDeriv()`. Note that, in most current implementations, the control 
 *	states are held constant over propagated arcs, thus, the time derivatives 
 *	are zero.
 *  
 *  To facilitate multiple shooting corrections, the partial derivatives that 
 *  relate the equations of motion to the core and control state variables are 
 *  required. The primary function of this object is to supply both the control 
 *  output, \f$ \vec{a}_{lt} \f$, and those partial derivatives. The first task 
 *  is achieved via the `getOutput()` function, as described above. The second 
 *  of the two tasks, the computation of the partial derivatives, is acheived in
 *  several pieces. The partial derivatives are collected in the 
 *  \f$ \mathbf{A} \f$ matrix that represents the linearized dynamics,
 *  \f[
 *  	\mathbf{A} =
 *		\begin{bmatrix}
 *			\mathbf{0}_{3\times 3} & \mathbf{I}_3 & \vec{0}_{3\times 1} & 
 *				\vec{0}_{3 \times n}\\[1em]
 *			\mathbf{D_{a,r}} & \mathbf{D_{a,v}} & \dpd{\vec{a}}{m} & 
 *				\mathbf{D_{a,\gamma}}\\[1em]
 *			\vec{0}_{1\times 3} & \vec{0}_{1\times 3} & 0 & 
 *				\dpd{\dot{m}}{\vec{\gamma}}\\[1em]
 *			\dpd{\dot{\vec{\gamma}}}{\vec{r}} & 
 *				\dpd{\dot{\vec{\gamma}}}{\vec{v}} & 
 *				\dpd{\dot{\vec{\gamma}}}{m} & 
 *				\dpd{\dot{\vec{\gamma}}}{\vec{\gamma}}
 *		\end{bmatrix}\,,
 *  \f]
 *  where
 *  \f{align*}{
 *  	\mathbf{D_{a,r}} &= \dpd{\vec{a}}{\vec{r}} = \dpd[2]{\Omega}{\vec{r}} + 
 *  		\dpd{\vec{a}_{lt}}{\vec{r}},\\
 *  	\mathbf{D_{a,v}} &= \dpd{\vec{a}}{\vec{v}} = \mathbf{K_v} + 
 *  		\dpd{\vec{a}_{lt}}{\vec{v}},\\
 *  	\mathbf{D_{a,\gamma}} &= \dpd{\vec{a}}{\vec{\gamma}} = 
 *  		\dmd{\Omega}{2}{\vec{r}}{}{\vec{\gamma}}{} +  
 *  		\dpd{\vec{a}_{lt}}{\vec{\gamma}} = \dpd{\vec{a}_{lt}}{\vec{\gamma}}\,.
 *  \f}
 *  These partials are computed via several functions. First, consider the 
 *  `getPartials_OutputWRTCoreState()` function,
 *  
 *  	getPartials_OutputWRTCoreState(t, s, pSys, partials, len)
 *  	
 *  This function fills the array, `partials`, with the elements of a 3x7 matrix,
 *  \f$ \partial \vec{a}_{lt}/\partial \vec{q}\f$, in row-major order. This 
 *  matrix includes pieces of several of the matrices listed above,
 *  \f[
 *  	\dpd{\vec{a}_{lt}}{\vec{q}} = 
 *  	\begin{bmatrix}
 *  		\dpd{\vec{a}_{lt}}{\vec{r}} & \dpd{\vec{a}_{lt}}{\vec{v}} & 
 *  		\dpd{\vec{a}_{lt}}{m}
 *  	\end{bmatrix}\,.
 *  \f]
 *  Next, consider the `getPartials_EOMsWRTCtrlState()` function,
 *  
 *  	getPartials_EOMsWRTCtrlState(t, s, pSys, partials, len)
 *  	
 *  which fills the `partials` array with the elements of a \f$ 7 \times n \f$ 
 *  matrix, \f$ \partial \vec{q}/\partial \vec{\gamma} \f$, which includes 
 *  several of the blocks from the \f$ \mathbf{A} \f$ matrix above:
 *  \f[
 *  	\dpd{\vec{q}}{\vec{\gamma}} = \begin{bmatrix}
 *  		\vec{0}_{3\times n} \\ 
 *  		\mathbf{D_{a,v}} \\[0.5em] 
 *  		\dpd{\dot{m}}{\vec{\gamma}}
 *  	\end{bmatrix}
 *  \f]
 *  Note that the top quadrant is practically always zero as the velocity EOMs 
 *  are always the simple \f$ \dot{\vec{r}} = \vec{v} \f$ set, i.e., they are 
 *  not a function of the control states.
 *  
 *  The final row of the \f$ \mathbf{A} \f$ matrix contains the partial 
 *  derivatives of the control state vector equation of motion with respect to 
 *  all core and control state variables. This entire row is computed via the 
 *  function `getPartials_TimeDerivWRTAllState()`,
 *  
 *  	getPartials_TimeDerivWRTAllState(t, s, pSys, partials, len)
 *  	
 *  where the `partials` array stores the elements in row-major order.
 */
class ControlLaw_cr3bp_lt : public ControlLaw{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	ControlLaw_cr3bp_lt(unsigned int id = NO_CTRL, 
		std::vector<double> params = {});

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	bool isVarMass() const;
	std::string getTypeString() const override;

	void setVarMass(bool);
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	double get_dmdt(double, const double*, const SysData*) const;
	double getThrustMag(double, const double*, const SysData*) const;
	void getOutput(double, const double*, const SysData*, 
		double*, unsigned int) const override;
	void getPartials_EOMsWRTCtrlState(double, const double*, 
		const SysData*, double*, unsigned int) const override;	
	void getPartials_OutputWRTCoreState(double, const double*, 
		const SysData*, double*, unsigned int) const override;
	
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	static std::string typeToString(unsigned int);
	static void convertLaws(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	static double thrust_dim2nondim(double, SysData_cr3bp_lt*);
	static double thrust_nondim2dim(double, SysData_cr3bp_lt*);
	static void pointingVecToAngles(Eigen::Vector3d, double*, double*);
	void print() const;
	//\}

	/**	
	 *	\name Control Law ID
	 *	
	 *	Low-Thrust Control Law IDs are generated via bitwise operations
	 *	
	 *	The ID is an unsigned int, which stores 16 bits, or 2 bytes. The bits
	 *	are split into four sections:
	 *	
	 *	ID = [ BASE | F_MAG | M | OPTIONS ]
	 *	
	 *	* BASE: 5-bits, represents the general law, usually a pointing strategy
	 *	* F_MAG: 3-bits, represents the thrust parameterization
	 *	* M: 2-bits, represents the mass parameterization
	 *	* OPTIONS: 6-bits, represents a series of options for the control law
	 *	
	 *	\{
	 */

	/** \brief masks 5 bits that represent base type */
	static const unsigned int BASE_MASK = 31 << 11;		// 11111 000 00 000000

	/** \brief masks 3 bits that represent the thrust parameterization */
	static const unsigned int F_MASK = 7 << 8;			// 00000 111 00 000000

	/** \brief masks 2 bits that represent the mass parameterization */
	static const unsigned int M_MASK = 3 << 6;			// 00000 000 11 000000

	/** \brief masks the 6 bits that represent the options */
	static const unsigned int OP_MASK = 63;				// 00000 000 00 111111

	/** \brief masks the first option */
	static const unsigned int OP1_MASK = 1 << 5;		// 00000 000 00 100000

	/** \brief masks the second option */
	static const unsigned int OP2_MASK = 1 << 4;		// 00000 000 00 010000

	/** \brief masks the third option */
	static const unsigned int OP3_MASK = 1 << 3;		// 00000 000 00 001000

	/** \brief masks the fourth option */
	static const unsigned int OP4_MASK = 1 << 2;		// 00000 000 00 000100

	/** \brief masks the fifth option */
	static const unsigned int OP5_MASK = 1 << 1;		// 00000 000 00 000010

	/** \brief masks the sixth option */
	static const unsigned int OP6_MASK = 1 << 0;		// 00000 000 00 000001

	/** \brief Represents a control law with arbitrary pointing in 3D space, 
	 * 	relative to the rotating frame.
	 *	
	 *	Constant parameters: `params = {}`
	 *	
	 *	Variable control states: `ctrl = {alpha, beta}`
	 *	
	 *	where `alpha` represents the in-plane pointing angle and `beta` represents 
	 *	the out-of-plane pointing angle (i.e., spherical coordinates). Alpha is 
	 *	zero when the thrust is in the XZ-plane with a positive x-component and 
	 *	increases with right-handed rotation about z. Beta is zero when the 
	 *	thrust is in the XY-plane and is positive when the z-component is 
	 *	positive.
	 *	
	 *	\see GEN_INERT
	 */
	static const unsigned int GENERAL = 1 << 11;		// 00001 000 00 000000

	/** \brief Represents a control law with arbitrary pointing in 3D space,
	 * 	relative to the inertial frame
	 * 	
	 * 	Constant parameters: `params = {theta0}`
	 * 	
	 * 	where `theta0` represents the orientation of the rotating frame w.r.t. 
	 * 	the inertial frame at time t = 0.
	 * 	
	 * 	Variable control states: `ctrl = {psi, beta}`
	 * 	
	 * 	where `psi` is the in-plane pointing angle relative to the inertial 
	 * 	X-axis and `beta` is the out of plane angle relative to the XY-plane.
	 * 	The inertial frame is assumed to be coplanar with the rotating frame 
	 * 	where inertial Z = rotating z.
	 * 
	 *  \see GENERAL
	 */
	static const unsigned int GEN_INERT = 4 << 11;		// 00100 000 00 000000

	/** \brief Represents a control law with thrust aligned with the velocity 
	 * 	vector.
	 * 
	 * 	Constant parameters: `params = {}`
	 * 	
	 * 	Variable control states: `ctrl = {}`
	 * 	
	 * 	To set the pointing to pro-velocity or anti-velocity, include the
	 * 	`VEL_PRO` or `VEL_ANTI` flags in the ID.
	 */
	static const unsigned int VEL_PT = 2 << 11;			// 00010 000 00 000000

	/** \brief Represents a control law that preserves the natural Hamiltonian, 
	 * 	i.e., the Jacobi Constant.
	 * 	
	 * 	Constant parameters: `params = {}`
	 * 	
	 * 	Variable control states: `ctrl = {}`
	 * 	
	 * 	To set the pointing to be left or right w.r.t. the velocity vector,
	 * 	include the `VEL_LEFT` or `VEL_RIGHT` flags in the ID.
	 */
	static const unsigned int CONST_C_2D = 3 << 11;		// 00011 000 00 000000
	

	/** \brief Represents a control law with constant thrust magnitude
	 * 
	 * 	The nondimensional thrust magnitude is appended to the vector of 
	 * 	constant parameters, e.g., `params = {f}`
	 */
	static const unsigned int CONST_F = 0 << 8;			// 00000 000 00 000000

	/** \brief Represents variable thrust with implicit bounds.
	 * 	
	 * 	The nondimensional maximum thrust magnitude is appended to the vector
	 * 	of constant parameters, e.g., `params = {fmax}`
	 * 	
	 * 	The nondimensional thrust magnitude is modeled via a sine function,
	 * 	\f[ f = \frac{1}{2} f_{\text{max}} \left( \sin\psi + 1 \right)\,, \f]
	 * 	where \f$\psi\f$, the control variable, is unbounded.
	 * 	Accordingly, the \f$\psi\f$ is added to the control state vector, e.g.,
	 * 	`ctrl = {alpha, beta, psi}`.
	 * 	
	 * 	The thrust magnitude is zero when 
	 * 	\f$\psi = -\frac{\pi}{2} \pm 2n\pi,~n=0,1,2,\dots\f$ and is maximized when
	 * 	\f$\psi = +\frac{\pi}{2} \pm 2n\pi,~n=0,1,2,\dots\f$.
	 * 	
	 * 	Due to the formulation, if the thrust is set to zero on an arc in a 
	 * 	multiple-shooting process, the thrust will never be updated. To avoid
	 * 	this feature, use small, nonzero magnitudes, e.g., 
	 * 	\f$\psi = -\pi/2 + \epsilon\f$
	 */
	static const unsigned int VAR_F_BND = 1 << 8;		// 00000 001 00 000000

	/** \brief Represents variable thrust with no implicit upper bound
	 * 
	 * 	The nondimensional maximum thrust magnitude is appended to the vector
	 * 	of constant parameters, e.g., `params = {fmax}`
	 * 	
	 * 	The nondimensional thrust magnitude is modeled via an exponential 
	 * 	function, \f[ f = f_{\text{max}} g^2\,, \f] where `g` is the control 
	 * 	variable and is appended to the control state vector, e.g.,
	 * 	`ctrl = {alpha, beta, g}`.
	 * 	
	 * 	In this formulation, thrust magnitude tends to zero as `g` tends to
	 * 	zero and there is no upper bound on the magnitude
	 */
	static const unsigned int VAR_F_UBND = 2 << 8;		// 00000 010 00 000000
	
	/** \brief Represents a variable mass parameterization for a CSI engine
	 * 
	 * 	Mass is modeled for a CSI engine with the differential equation,
	 * 	\f[ \dot{m} = \frac{-f}{I_{sp} g_0}\,, \f] where `f` is the 
	 * 	nondimensional thrust magnitude, `Isp` is the constant specific impulse,
	 * 	and `g0` is the average Earth gravitational constant. Append `Isp` to
	 * 	the parameter vector, e.g., `params = {f, Isp}`.
	 * 	
	 * 	No control states are required for this formulation
	 */
	static const unsigned int CSI_VAR_M = 0 << 6;		// 00000 000 00 000000

	/** \brief Represents a constant mass parameterization.
	 *
	 *	The model is quite simple, \f[ \dot{m} = 0 \f]
	 *	
	 *	No additional constant parameters or variable control states are 
	 *	required for this formulation
	 * 
	 */
	static const unsigned int CONST_M = 1 << 6;		// 00000 000 01 000000
	
	/** \brief Option for VEL_PT control to orient thrust with +v vector */
	static const unsigned int VEL_PRO = 0 << 5;			// 00000 000 00 000000

	/** \brief Option for VEL_PT control to orient thrust with -v vector */	
	static const unsigned int VEL_ANTI = 1 << 5;		// 00000 000 00 100000
	
	/** \brief Option for CONST_C_2D control to orient thrust left of v vector */
	static const unsigned int VEL_LEFT = 0 << 5;		// 00000 000 00 000000

	/** \brief Option for CONST_C_2D control to orient thrust right of v vector */
	static const unsigned int VEL_RIGHT = 1 << 5;		// 00000 000 00 100000

	//\}

	/**
	 *  @brief Identify the control law
	 *  @delails For all of these control laws, the `getLaw()` function
	 *  returns a 3-dimensional thrust direction unit vector. Similarly, the 
	 *  `getPartials_State()` function returns 21 derivative values
	 *  that relate the thrust direction to the 7 core CR3BP-LT states.
	 *  
	 *  Laws are split into categories based on the type value:
	 *  * 1-99: 		Constant thrust, variable mass
	 *  * 101 - 199: 	Constant thrust, constant mass
	 *  * 1001 - 1099: 	Variable thrust, variable mass
	 *  
	 *  Procedure for adding a new control law:
	 *  * Add the type to the enumerated type, `Law_tp`, with full documentation
	 *  * Add to `init()` and `typeToString()`
	 *  * Add to `get_dmdt()`
	 *  * Add to `getOutput()` and any functions called from that switchboard
	 *  * Add to `getPartials_OutputWRTCoreState()` and any functions called 
	 *  	from that switchboard
	 *  * Add to `getPartials_EOMsWRTCtrlState()` and any functions called from 
	 *  	that switchboard
	 */
	enum Law_tp : unsigned int{
		CONST_FC_2D_LEFT = CONST_C_2D | CONST_F | CSI_VAR_M | VEL_LEFT,
		CONST_FC_2D_RIGHT = CONST_C_2D | CONST_F | CSI_VAR_M | VEL_RIGHT,
		CONST_F_PRO_VEL = VEL_PT | CONST_F | CSI_VAR_M | VEL_PRO,
		CONST_F_ANTI_VEL = VEL_PT | CONST_F | CSI_VAR_M | VEL_ANTI,
		CONST_F_GENERAL = GENERAL | CONST_F | CSI_VAR_M,
		CONST_MF_GENERAL = GENERAL | CONST_F | CONST_M,
		VAR_F_CONST_C_2D_LEFT = CONST_C_2D | VAR_F_BND | CSI_VAR_M | VEL_LEFT,
		VAR_F_CONST_C_2D_RIGHT = CONST_C_2D | VAR_F_BND | CSI_VAR_M | VEL_RIGHT,
		VAR_F_PRO_VEL = VEL_PT | VAR_F_BND | CSI_VAR_M | VEL_PRO,
		VAR_F_ANTI_VEL = VEL_PT | VAR_F_BND | CSI_VAR_M | VEL_ANTI,
		VAR_F_GENERAL = GENERAL | VAR_F_BND | CSI_VAR_M
	};
protected:

	/**
	 *  \name *structors
	 *  \{
	 */
	void init() override;
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void getAccel_AlongVel(double, const double*, const SysData_cr3bp_lt*,
		double*, unsigned int) const;
	void getAccel_ConstC_2D(double, const double*, const SysData_cr3bp_lt*,
		double*, unsigned int) const;
	void getAccel_GeneralRot(double, const double*, const SysData_cr3bp_lt*,
		double*, unsigned int) const;
	void getAccel_GeneralInert(double, const double*, const SysData_cr3bp_lt*,
		double*, unsigned int) const;

	void getPartials_AccelWRTCore_AlongVel(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getPartials_AccelWRTCore_ConstC_2D(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getPartials_AccelWRTCore_GeneralRot(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getPartials_AccelWRTCore_GeneralInert(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getPartials_EOMsWRTCtrl_GeneralDir(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getPartials_EOMsWRTCtrl_VarF(double, const double*,
		const SysData_cr3bp_lt*, double*, unsigned int) const;

	static void convertTo_GeneralConstF(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	//\}
};

}// End of astrohelion namespace