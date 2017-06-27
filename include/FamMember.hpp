/**
 * \file FamMember.hpp
 * \brief Defines an interface to access information specific to
 * an orbit that is a member of a family.
 * 
 * \details Derived classes are children of this interface as 
 * well as one of the system-specific Arcset objects
 * 
 * \author Andrew Cox
 * \version June 27, 2017
 * \copyright GNU GPL v3.0
 */

#pragma once

#include "Core.hpp"

#include <vector>

#include "Common.hpp"
#include "EigenDefs.hpp"

namespace astrohelion{

static const char PARAMKEY_FAM_XAMP[] = "Ax";
static const char PARAMKEY_FAM_YAMP[] = "Ay";
static const char PARAMKEY_FAM_ZAMP[] = "Az";

class FamMember : public Core{

public:
	/**
	 *  \name *structors
	 *  \{
	 */
	FamMember(){}
	
	FamMember(const FamMember &m){
		copyMe(m);
		return *this;
	}//================================================

	virtual ~FamMember(){}
	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	std::vector<cdouble> getEigVals() const{ return eigVals; }
	MatrixXRcd getEigVecs() const { return eigVecs; }
	
	double getXAmplitude() const{
		return params.at(PARAMKEY_FAM_XAMP);
	}//================================================
	
	double getYAmplitude() const{
		return params.at(PARAMKEY_FAM_YAMP);
	}//================================================
	
	double getZAmplitude() const{
		return params.at(PARAMKEY_FAM_ZAMP);
	}//================================================

	virtual std::vector<double> getIC() const = 0;
	virtual MatrixXRd getSTM() const = 0;

	void setEigVals(std::vector<cdouble> vals){ eigVals = vals; }
	void setEigVecs(MatrixXRcd vecs){ eigVecs = vecs; }
	void setXAmplitude(double amp){ params[PARAMKEY_FAM_XAMP] = amp; }
	void setYAmplitude(double amp){ params[PARAMKEY_FAM_YAMP] = amp; }
	void setZAmplitude(double amp){ params[PARAMKEY_FAM_ZAMP] = amp; }

	//\}
	
protected:
	MatrixXRcd eigVecs = MatrixXRcd::Zero(6,6);	//!< Vector of eigenvectors (as columns)
	std::vector<cdouble> eigVals(6, {NAN, 0});	//!< Vector of eigenvalues
	std::map<std::string, double> params {};	//!< Map of parameters that describe the family member

	void copyMe(const FamMember &m){
		eigVecs = m.eigVecs;
		eigVals = m.eigVals;
		params = m.params;
	}//================================================
};

}//End of astrohelion namespace