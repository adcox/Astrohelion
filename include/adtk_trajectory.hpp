/**
 *	Trajectory Class:
 *
 *	TODO:
 *		- Overload +, += operators for concatenating trajectories
 */
#ifndef __H_TRAJECTORY_
#define __H_TRAJECTORY_

#include "adtk_sys_data.hpp"
#include "adtk_matrix.hpp"

#include <vector>

class adtk_trajectory{
	static const int STATE_WIDTH = 9;

	protected:
		/** Number of points along integrated path */
		int numPoints;

		/** Holds state info: [pos, vel, accel] in 1D form, so every 9 elements
		 * 	constitutes a new "row" 
		 */
		std::vector<double> state;

		/** Holds time info */
		std::vector<double> times;

		/** An array containing the STM at each step of the integration */
		std::vector<adtk_matrix> allSTM;
	public:

		adtk_trajectory();
		adtk_trajectory(int);

		adtk_trajectory& operator= (const adtk_trajectory&);

		double getLength();
		std::vector<double> getState(int);
		std::vector<double>* getState();
		double getTime(int);
		std::vector<double>* getTime();
		adtk_matrix getSTM(int);
		std::vector<adtk_matrix>* getSTM();

		void setState(std::vector<double>);
		void setTime(std::vector<double>);
		void setSTMs(std::vector<adtk_matrix>);
};

#endif