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
class tpat_segment : public tpat_linkable{

public:
	static const int ORIG_IX = 0;	//!< Index of the origin node in the links array
	static const int TERM_IX = 1;	//!< Index if the terminus node in the links array
	
	// *structors
	tpat_segment();
	tpat_segment(int, int, double);
	tpat_segment(int, int, double, const double[36]);
	tpat_segment(const tpat_segment&);
	// ~tpat_segment();

	// Operators
	tpat_segment& operator =(const tpat_segment&);
	friend bool operator ==(const tpat_segment&, const tpat_segment&);
	friend bool operator !=(const tpat_segment&, const tpat_segment&);

	// Set and Get Functions
	void addConstraint(tpat_constraint);
	void clearConstraints();
	std::vector<tpat_constraint> getConstraints() const;
	int getNumCons() const;
	int getOrigin() const;
	MatrixXRd getSTM() const;
	std::vector<double> getSTMElements() const;
	int getTerminus() const;
	double getTOF() const;
	void removeConstraint(int);
	void setConstraints(std::vector<tpat_constraint>);
	void setOrigin(int);
	void setTerminus(int);
	void setTOF(double);
	void setSTM(MatrixXRd);
	void setSTM(const double*);
	void setSTM(std::vector<double>);

protected:
	virtual void copyMe(const tpat_segment&);

	double tof = 0;			//!< Time-of-flight along this segment, units consistent with the system
	double stm[36] = {	1, 0, 0, 0, 0, 0,
						0, 1, 0, 0, 0, 0,
						0, 0, 1, 0, 0, 0,
						0, 0, 0, 1, 0, 0,
						0, 0, 0, 0, 1, 0,
						0, 0, 0, 0, 0, 1 };			//!< 6x6 STM, stored in row-major order

	/** Stores constraints on this segment */
	std::vector<tpat_constraint> cons {};
};

#endif