/**
 * \file FamMember_cr3bp.hpp
 * 
 * \author Andrew Cox
 * \version June 27, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

namespace astrohelion{

#include "FamMember.hpp"
#include "Arcset_cr3bp.hpp"

// Forward Declarations
class SysData_cr3bp;

static const char PARAMKEY_FAM_JC = "JC";

class FamMember_cr3bp : public FamMember, public Arcset_cr3bp{
public:
	/**
	 *  \name *structors
	 *  \{
	 */
	FamMember_cr3bp(const SysData_cr3bp *pSys) : FamMember(), Arcset_cr3bp(pSys) {}
	FamMember_cr3bp(const FamMember_cr3bp &m) : FamMember(m), Arcset_cr3bp(m) {}
	~FamMember_cr3bp(){}
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	double getJacobi() const { return params.at(PARAMKEY_FAM_JC); }
	
	std::vector<double> getIC() const{
		if(nodes.size() == 0)
			throw Exception("FamMember_cr3bp::getIC(): There are no nodes... cannot retrieve IC");

		putInChronoOrder();	// No action if already in chronological order
		return getStateByIx(0);
	}//================================================

	MatrixXRd getSTM() const{
		if(segs.size() == 0)
			throw Exception("FamMember_cr3bp::getSTM(): There are no segments... cannot returieve the STM");

		putInChronoOrder(); // no action if already in chronological order
		return getSTMByIx(-1);
	}//================================================

	//\}

};

}// End of astrohelion namespace