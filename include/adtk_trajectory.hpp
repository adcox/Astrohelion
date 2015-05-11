/**
 *	Trajectory Class:
 *
 *	TODO:
 *		- Overload +, += operators for concatenating trajectories
 */
#ifndef __H_TRAJECTORY_
#define __H_TRAJECTORY_

class adtk_trajectory{
	static const int STATE_WIDTH = 10;

	private:
		/** Holds state info: [time, pos, vel, accel] */
		double state[][STATE_WIDTH];

		/** A 3D array containing the STM at each step of the integration */
		double STM[][6][6];

		/** An array of Jacobi values at each step of the integration;
		only applies for autonomous systems like the CR3BP
		*/
		double Jacobi[];

	public:

		adtk_trajectory();
		adtk_trajectory(double *state);

		double* getState(int);
		double* getAllStates();
		double getJacobi(int);
		double* getAllJacobi();
		double* getSTM(int);
		double* getAllSTM();
};

#endif