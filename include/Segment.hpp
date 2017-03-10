/**
 *  \file Segment.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
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

#include "Linkable.hpp"

#include "Constraint.hpp"
#include "ControlLaw.hpp"
#include "EigenDefs.hpp"

#include <cmath>
#include <vector>

namespace astrohelion{

// Forward Declarations

/**
 *	\ingroup traj
 *	\brief Links nodes together to describe the flow of a trajectory
 *
 *	\author Andrew Cox
 *	\version May 1, 2016
 *	\copyright GNU GPL v3.0
 */
class Segment : public Linkable{

public:
	static const int ORIG_IX = 0;	//!< Index of the origin node in the links array
	static const int TERM_IX = 1;	//!< Index if the terminus node in the links array
	
	/**
	 *  \name *structors
	 *  \{
	 */
	Segment();
	Segment(int, int, double);
	Segment(int, int, double, const double*, unsigned int, unsigned int ctrlLawID = ControlLaw::NO_CTRL);
	Segment(const Segment&);
	// ~Segment();
	//\}

	// Operators
	Segment& operator =(const Segment&);
	friend bool operator ==(const Segment&, const Segment&);
	friend bool operator !=(const Segment&, const Segment&);

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	void addConstraint(Constraint);
	void clearConstraints();
	std::vector<Constraint> getConstraints() const;
	unsigned int getCtrlLaw() const;
	int getNumCons() const;
	int getOrigin() const;
	MatrixXRd getSTM() const;
	int getTerminus() const;
	double getTOF() const;
	std::vector<bool> getVelCon() const;
	
	void removeConstraint(int);
	void setConstraints(std::vector<Constraint>);
	void setCtrlLaw(unsigned int);
	void setID(int) override;
	void setOrigin(int);
	void setTerminus(int);
	void setTOF(double);
	void setSTM(MatrixXRd);
	void setSTM(const double*, unsigned int);
	void setSTM(std::vector<double>);
	void setVel_AllCon();
	void setVel_AllDiscon();
	void setVelCon(const bool[3]);
	void setVelCon(std::vector<bool>);
	void setVelCon(bool, bool, bool);
	//\}

	void print() const;
protected:
	virtual void copyMe(const Segment&);

	/** Stores flags, which are currently used to indicate continuity */
	std::vector<bool> flags {true, true, true};
	
	unsigned int ctrlLawID = ControlLaw::NO_CTRL;	//!< Control law applied during this segment; by default, set to NO_CTRL

	double tof = 0;			//!< Time-of-flight along this segment, units consistent with the system
	MatrixXRd stm = MatrixXRd::Identity(6,6);	// Initialize as 6x6 (most common use), but can be easily resized

	/** Stores constraints on this segment */
	std::vector<Constraint> cons {};
};

}// END of Astrohelion namespace