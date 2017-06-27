/**
 * \file Family.hpp
 */

#pragma once

#include "Core.hpp"

namespace astrohelion{

/**
 *	\brief How to sort members of this family
 *	
 *	Families typically evolve naturally through one variable, and
 *	it makes sense to analyze trends using this parameter as the 
 *	independent variable.
 *
 *	The integer value of each type is set up to map to the index
 *	of the independent variable in the continuation process
 */
enum class FamSort_tp : int{
	SORT_X 		= 0,	//!< Sort by the x-coordinate in the IC
	SORT_Y 		= 1,	//!< Sort by the y-coordinate in the IC
	SORT_Z 		= 2,	//!< Sort by the z-coordinate in the IC
	SORT_VX 	= 3,	//!< Sort by the x velocity component in the IC
	SORT_VY  	= 4,	//!< Sort by the y velocity component in the IC
	SORT_VZ 	= 5,	//!< Sort by the z velocity component in the IC
	SORT_TOF 	= 6,	//!< Sort by Time-Of-Flight
	SORT_JC 	= 7,	//!< Sort by Jacobi Constant
	SORT_NONE 	= 999	//!< Do not adjust sorting; sortMembers() will do nothing
};

class Family : public Core{
public:
	Family();
	virtual ~Family();

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	std::string getName() const;
	int getNumMembers() const;
	FamSort_tp getSortType() const;
	const char* getSortTypeStr() const;

	void setName(std::string);
	void setSortType(FamSort_tp);
	//\}
	
	std::vector<int> findBifurcations();
	void saveToMat(const char*);
	void sortEigs();
	void sortMembers();
};

}// End of astrohelion namespace