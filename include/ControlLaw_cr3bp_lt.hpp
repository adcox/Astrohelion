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
 *  @details As specified in the MathSpec document, the CR3BP-LT equations of motion
 *  are given by
 *  \f{align*}{
 *  	\vec{a} &= \mathbf{K_v} \vec{f} + \dpd{\Omega}{\vec{r}} + \vec{a}_{lt}, \qquad \vec{a}_lt = \frac{f}{m} \hat{u}\\
 *  	\dot{m} &= \dod{m}{\tau}\\
 *  	\dot{\vec{\gamma}} &= \dod{\vec{\gamma}}{\tau}
 *  \f}
 *  where \f$ \vec{a} = \{\ddot{x}\,\, \ddot{y}\,\, \ddot{z} \}\f$ is the rotating acceleration vector,
 *  \f$ \vec{v} = \{\dot{x}\,\, \dot{y}\,\, \dot{z} \} \f$ is the rotating velocity vector, and
 *  \f$ \vec{r} = \{x\,\, y\,\, z \} \f$ is the spacecraft position vector in the rotating frame. The CR3BP
 *  pseudopotential is represented by \f$ \Omega \f$, the nondimensional low-thrust magnitude if \f$ f \f$,
 *  the nondimensional spacecraft mass is \f$ m \in [0,1] \f$, and the thrust pointing unit vector is
 *  \f$ \hat{u} \f$. While most of the core state eoms, \f$ \vec{a} \f$, are computed in the DynamicsModel_cr3bp_lt
 *  class, several of the control-specific elements are computed in the control law object. The low-thrust 
 *  acceleration, \f$ \vec{a}_{lt} \f$, is obtained from the `getLaw_output()` function,
 *  
 *  	a_lt = getLaw_Output(t, s, pSysData, partials, len);
 *  	
 *  where the `s` array contains a full state, which is generally of the form
 *  \f$ \vec{s} = \begin{Bmatrix} \vec{q} & \vec{\gamma} & \phi_1 & \dots & \phi_n & \vec{q}_{\text{extra}} \end{Bmatrix} \f$.
 *  In this vector, the core and control state variables are
 *  \f[
 *  	\vec{q} = \begin{Bmatrix} \vec{r} \,\, \vec{v} \,\, m \end{Bmatrix}^T,\qquad \vec{\gamma} \in \mathbb{R}^{n\times 1},
 *  \f]
 *  respectively, \f$ \phi_i \f$ are the state transition matrix elements in row-major order, and 
 *  \f$ \vec{q}_{\text{extra}} \f$ are the "extra" state variables included in some models (not in this model, however).
 *  
 *	Because the mass time derivative is a function of thrust and Isp, both control law parameters
 *  (or states), the second equation of motion is computed via the control law function `get_dmdt()`. The
 *  final EOM, the control state time derivative, is computed via the control law function `getLaw_StateDeriv()`.
 *  Note that, in most current implementations, the control states are held constant over propagated arcs, thus,
 *  the time derivatives are zero.
 *  
 *  To facilitate multiple shooting corrections, the partial derivatives that relate the equations
 *  of motion to the core and control state variables are required. The primary function of this object is to supply
 *  both the control output, \f$ \vec{a}_{lt} \f$, and those partial derivatives. The first task is achieved via
 *  the `getLaw_Output()` function, as described above. The second of the two tasks, the computation of the partial 
 *  derivatives, is acheived in several pieces. The partial derivatives are collected in the 
 *  \f$ \mathbf{A} \f$ matrix that represents the linearized dynamics,
 *  \f[
 *  	\mathbf{A} =
 *		\begin{bmatrix}
 *			\mathbf{0}_{3\times 3} & \mathbf{I}_3 & \vec{0}_{3\times 1} & \vec{0}_{3 \times n}\\[1em]
 *			\mathbf{D_{a,r}} & \mathbf{D_{a,v}} & \dpd{\vec{a}}{m} & \mathbf{D_{a,\gamma}}\\[1em]
 *			\vec{0}_{1\times 3} & \vec{0}_{1\times 3} & 0 & \dpd{\dot{m}}{\vec{\gamma}}\\[1em]
 *			\dpd{\dot{\vec{\gamma}}}{\vec{r}} & \dpd{\dot{\vec{\gamma}}}{\vec{v}} & \dpd{\dot{\vec{\gamma}}}{m} & \dpd{\dot{\vec{\gamma}}}{\vec{\gamma}}
 *		\end{bmatrix}\,,
 *  \f]
 *  where
 *  \f{align*}{
 *  	\mathbf{D_{a,r}} &= \dpd{\vec{a}}{\vec{r}} = \dpd[2]{\Omega}{\vec{r}} + \dpd{\vec{a}_{lt}}{\vec{r}},\\
 *  	\mathbf{D_{a,v}} &= \dpd{\vec{a}}{\vec{v}} = \mathbf{K_v} + \dpd{\vec{a}_{lt}}{\vec{v}},\\
 *  	\mathbf{D_{a,\gamma}} &= \dpd{\vec{a}}{\vec{\gamma}} = \dmd{\Omega}{2}{\vec{r}}{}{\vec{\gamma}}{} +  \dpd{\vec{a}_{lt}}{\vec{\gamma}} = \dpd{\vec{a}_{lt}}{\vec{\gamma}}\,.
 *  \f}
 *  These partials are computed via several functions. First, consider the `getLaw_OutputPartials()` function,
 *  
 *  	getLaw_OutputPartials(t, s, pSys, partials, len)
 *  	
 *  This function fills the array, `partials`, with the elements of a 3x7 matrix,
 *  \f$ \partial \vec{a}_{lt}/\partial \vec{q}\f$, in row-major order. This matrix includes pieces of 
 *  several of the matrices listed above,
 *  \f[
 *  	\dpd{\vec{a}_{lt}}{\vec{q}} = \begin{bmatrix} \dpd{\vec{a}_{lt}}{\vec{r}} & \dpd{\vec{a}_{lt}}{\vec{v}} & \dpd{\vec{a}_{lt}}{m} \end{bmatrix}\,.
 *  \f]
 *  Next, consider the `getLaw_EOMPartials()` function,
 *  
 *  	getLaw_EOMPartials(t, s, pSys, partials, len)
 *  	
 *  which fills the `partials` array with the elements of a \f$ 7 \times n \f$ matrix,
 *  \f$ \partial \vec{q}/\partial \vec{\gamma} \f$, which includes several of the blocks from the
 *  \f$ \mathbf{A} \f$ matrix above:
 *  \f[
 *  	\dpd{\vec{q}}{\vec{\gamma}} = \begin{bmatrix}
 *  		\vec{0}_{3\times n} \\ \mathbf{D_{a,v}} \\[0.5em] \dpd{\dot{m}}{\vec{\gamma}}
 *  	\end{bmatrix}
 *  \f]
 *  Note that the top quadrant is practically always zero as the velocity EOMs are always the simple
 *  \f$ \dot{\vec{r}} = \vec{v} \f$ set, i.e., they are not a function of the control states.
 *  
 *  The final row of the \f$ \mathbf{A} \f$ matrix contains the partial derivatives of the control state
 *  vector equation of motion with respect to all core and control state variables. This entire row is 
 *  computed via the function `getLaw_StateDerivPartials()`,
 *  
 *  	getLaw_StateDerivPartials(t, s, pSys, partials, len)
 *  	
 *  where the `partials` array stores the elements in row-major order.
 */
class ControlLaw_cr3bp_lt : public ControlLaw{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	ControlLaw_cr3bp_lt(unsigned int id = NO_CTRL, std::vector<double> params = {});
	// ControlLaw_cr3bp_lt(unsigned int, double, double);

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	std::string getLawTypeString() const override;
	//\}

	/**
	 *  \name Analysis Functions
	 *  \{
	 */
	void getLaw_EOMPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	void getLaw_Output(double t, const double *s, const SysData *sysData, double *law, unsigned int len) const override;
	void getLaw_OutputPartials(double t, const double *s, const SysData *pSys, double *partials, unsigned int len) const override;
	double get_dmdt(double t, const double *s, const SysData *pSys) const;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	static std::string lawTypeToString(unsigned int);
	static void convertLaws(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	static double thrust_dim2nondim(double, SysData_cr3bp_lt*);
	static double thrust_nondim2dim(double, SysData_cr3bp_lt*);
	static void pointingVecToAngles(Eigen::Vector3d, double*, double*);
	void print() const;
	//\}

	/**
	 *  @brief Identify the control law
	 *  \delails For all of these control laws, the `getLaw()` function
	 *  returns a 3-dimensional thrust direction unit vector. Similarly, the 
	 *  `getPartials_State()` function returns 21 derivative values
	 *  that relate the thrust direction to the 7 core CR3BP-LT states.
	 */
	enum Law_tp : unsigned int{
		CONST_F_C_2D_LEFT = 1,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control,
									 * thrust left w.r.t. velocity direction. Thrust magnitude is constant.
									 * - params: { thrust (nondim), Isp (seconds) }
									 * - ctrlStates: {} (None)
									 */
		CONST_F_C_2D_RIGHT = 2,		/*!< Jacobi-Preserving (constant C), two-dimensional (xy-planar) control, 
									 * thrust right w.r.t. velocity direction. Thrust magnitude is constant.
									 * - params: { thrust (nondim), Isp (seconds) }
									 * - ctrlStates: {} (None)
									 */
		CONST_F_PRO_VEL = 3,		/*!< Thrust along velocity vector (maximum energy increase). Thrust
									 * 	magnitude is constant.
									 * - params: { thrust (nondim), Isp (seconds) }
									 * - ctrlStates: {} (None)
									 */
		CONST_F_ANTI_VEL = 4,		/*!< Thrust along anti-velocity vector (maximum energy decrease). Thrust
									 * 	magnitude is constant.
									 * - params: { thrust (nondim), Isp (seconds) }
									 * - ctrlStates: {} (None)
									 */
		CONST_F_GENERAL = 5,		/*!< Thrust in an arbitrary direction. Thrust magnitude is constant.
									 * - params: { thrust (nondim), Isp (seconds) }
									 * - ctrlStates: {alpha (rad), beta (rad)}
									 */
		VAR_F_CONST_C_2D_LEFT = 101,/*!< Jacobi-preserving (constant C), two-dimensional (xy-planar) control,
									 * thrust left w.r.t. velocity direction.
									 * - params: {Isp (seconds)}
									 * - ctrlStates: {f (nondim)}
									 */
		VAR_F_CONST_C_2D_RIGHT = 102,/*!< Jacobi-preserving (constant C), two-dimensional (xy-planar) control,
									 * thrust right w.r.t. velocity direction.
									 * - params: {Isp (seconds)}
									 * - ctrlStates: {f (nondim)}
									 */
		VAR_F_PRO_VEL = 103,		/*!< Thrust along velocity vector (maximum energy increase).
									 * - params: {Isp (seconds)}
									 * - ctrlStates: {f (nondim)}
									 */
		VAR_F_ANTI_VEL = 104,		/*!< Thrust along anti-velocity vector (maximum energy decrease).
									 * - params: {Isp (seconds)}
									 * - ctrlStates: {f (nondim)}
									 */
		VAR_F_GENERAL = 105			/*!< Thrust in an arbitrary direction.
									 * - params: {Isp (seconds)}
									 * - ctrlStates: {f (nondim), alpha (rad), beta (rad)}
									 */
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
	void getAccel_AlongVel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getAccel_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getAccel_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getAccelPartials_AlongVel(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getAccelPartials_ConstC_2D(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getAccelPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;

	void getEOMPartials_GeneralDir(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	void getEOMPartials_VarF(double, const double*, const SysData_cr3bp_lt*, double*, unsigned int) const;
	
	static void convertTo_GeneralConstF(Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*);
	//\}
};

}// End of astrohelion namespace