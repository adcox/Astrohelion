/**
 *	@file adtk_bcr4bpr_nodeset.cpp
 */

#include "adtk_bcr4bpr_nodeset.hpp"

adtk_bcr4bpr_constraint adtk_bcr4bpr_nodeset::getConstraint(int i){ return constraints.at(i); }

void adtk_bcr4bpr_nodeset::addConstraint(adtk_bcr4bpr_constraint c){ constraints.push_back(c); }