/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#ifndef H_TPAT_SEG
#define H_TPAT_SEG

#include "tpat_linkable.hpp"

#include "tpat_constraint.hpp"
#include "tpat_eigen_defs.hpp"

#include <cmath>
#include <vector>

// Forward Declarations

/**
 *	@brief A brief description
 *
 *	@author Andrew Cox
 *	@version May 1, 2016
 *	@copyright GNU GPL v3.0
 */
class TPAT_Segment : public TPAT_Linkable{

public:
	static const int ORIG_IX = 0;	//!< Index of the origin node in the links array
	static const int TERM_IX = 1;	//!< Index if the terminus node in the links array
	
	// *structors
	TPAT_Segment();
	TPAT_Segment(int, int, double);
	TPAT_Segment(int, int, double, const double[36]);
	TPAT_Segment(const TPAT_Segment&);
	// ~TPAT_Segment();

	// Operators
	TPAT_Segment& operator =(const TPAT_Segment&);
	friend bool operator ==(const TPAT_Segment&, const TPAT_Segment&);
	friend bool operator !=(const TPAT_Segment&, const TPAT_Segment&);

	// Set and Get Functions
	void addConstraint(TPAT_Constraint);
	void clearConstraints();
	std::vector<TPAT_Constraint> getConstraints() const;
	int getNumCons() const;
	int getOrigin() const;
	MatrixXRd getSTM() const;
	std::vector<double> getSTMElements() const;
	int getTerminus() const;
	double getTOF() const;
	std::vector<bool> getVelCon() const;
	void removeConstraint(int);
	void setConstraints(std::vector<TPAT_Constraint>);
	void setOrigin(int);
	void setTerminus(int);
	void setTOF(double);
	void setSTM(MatrixXRd);
	void setSTM(const double*);
	void setSTM(std::vector<double>);
	void setVel_AllCon();
	void setVel_AllDiscon();
	void setVelCon(const bool[3]);
	void setVelCon(std::vector<bool>);
	void setVelCon(bool, bool, bool);

protected:
	virtual void copyMe(const TPAT_Segment&);

	/** Stores flags, which are currently used to indicate continuity */
	std::vector<bool> flags {true, true, true};

	double tof = 0;			//!< Time-of-flight along this segment, units consistent with the system
	double stm[36] = {	1, 0, 0, 0, 0, 0,
						0, 1, 0, 0, 0, 0,
						0, 0, 1, 0, 0, 0,
						0, 0, 0, 1, 0, 0,
						0, 0, 0, 0, 1, 0,
						0, 0, 0, 0, 0, 1 };			//!< 6x6 STM, stored in row-major order

	/** Stores constraints on this segment */
	std::vector<TPAT_Constraint> cons {};
};

#endif