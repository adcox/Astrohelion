/**
 *  \file Arcset_periodic.cpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version October 5, 2017
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

#include <Eigen/Eigenvalues>

#include "Arcset_periodic.hpp"

#include "Exceptions.hpp"

namespace astrohelion{

//-----------------------------------------------------------------------------
//      Constructors and Desctructor
//-----------------------------------------------------------------------------

Arcset_periodic::Arcset_periodic(const SysData *pSys) : Arcset(pSys) {}
Arcset_periodic::Arcset_periodic(const Arcset_periodic &a) : Arcset(a) {}
Arcset_periodic::Arcset_periodic(const Arcset &a) : Arcset(a) {}
Arcset_periodic::~Arcset_periodic() {}

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  \brief Compute the eigenvalues and eigenvectors of the monodromy matrix
 *
 *	\param pVals pointer to a vector of complex doubles to store the eigenvalues in
 *	\param pVecs pointer to a matrix to store the eigenvectors in (as columns)
 *	
 *	\throws Exception if either input pointer is NULL
 */
void Arcset_periodic::getEigData(std::vector<cdouble> *pVals, MatrixXRcd *pVecs){
	if(pVals == nullptr && pVecs == nullptr)
		return;

	// Compute the eigen data
	MatrixXRd mono = getMonodromy();
	Eigen::EigenSolver<MatrixXRd> eigensolver(mono);
    if(eigensolver.info() != Eigen::Success)
        throw Exception("FamGenerator::cr3bp_pseudoArcCont: Could not compute eigenvalues of monodromy matrix");

    // Extract the eigen data to return to the user
    if(pVals){
    	pVals->clear();		// Delete any existing data
    	Eigen::VectorXcd vals = eigensolver.eigenvalues();
    	pVals->insert(pVals->begin(), vals.data(), vals.data() + vals.size());
    }

    if(pVecs){
    	*pVecs = eigensolver.eigenvectors();
    }
}//====================================================

/**
 *  \brief Retrieve the monodromy matrix for this periodic orbit
 *  \return the monodromy matrix
 */
MatrixXRd Arcset_periodic::getMonodromy(){
	setSTMs_cumulative();	// Puts arcset in chronological order if it isn't already and computes the cumulative STMs
	return getSTMByIx(-1);	// the final STM represents the transition for the full orbit
}//====================================================

/**
 *  \brief Retrieve the x-amplitude of the orbit
 *  \details the amplitude is the difference between the maximum and minimum extents
 *  \return the x-amplitude of the orbit
 */
double Arcset_periodic::getXAmp() const { return getAmp(0); }

/**
 *  \brief Retrieve the y-amplitude of the orbit
 *  \details the amplitude is the difference between the maximum and minimum extents
 *  \return the y-amplitude of the orbit
 */
double Arcset_periodic::getYAmp() const { return getAmp(1); }

/**
 *  \brief Retrieve the z-amplitude of the orbit
 *  \details the amplitude is the difference between the maximum and minimum extents
 *  \return the z-amplitude of the orbit
 */
double Arcset_periodic::getZAmp() const { return getAmp(2); }

/**
 *  \brief Retrieve the amplitude of the orbit measured across a specific state
 *  \details the amplitude is the difference between the maximum and minimum extents
 * 
 *  \param stateIx the index of the state variable (e.g., x = 0, y = 1, etc.)
 *  \return the amplitude of the oscillation in the specified variable
 */
double Arcset_periodic::getAmp(unsigned int stateIx) const{

	double minVal = 0, maxVal = 0, val = 0;
	for(unsigned int s = 0; s < segs.size(); s++){
		std::vector<double> states = segs[s].getStateVector();
		unsigned int w = segs[s].getStateWidth();

		if(stateIx >= w){
			char msg[128];
			sprintf(msg, "Arcset_periodic::getAmp: stateIx = %u is out of bounds on segment ID = %u",
				stateIx, segs[s].getID());
			throw Exception(msg);
		}

		for(unsigned int i = 0; i < states.size()/w; i++){
			val = states[w*i + stateIx];
			if(s == 0 && i == 0){
				minVal = val;
				maxVal = minVal;
			}else{
				if(val < minVal){ minVal = val; }
				if(val > maxVal){ maxVal = val; }
			}
		}
	}

	return maxVal - minVal;
}//====================================================

}// End of astrohelion namespace