/**
 *	@file adtk_cr3bp_nodeset.cpp
 */

#include "adtk_cr3bp_nodeset.hpp"

adtk_cr3bp_constraint adtk_cr3bp_nodeset::getConstraint(int i){ return constraints.at(i); }

void adtk_cr3bp_nodeset::addConstraint(adtk_cr3bp_constraint c){ constraints.push_back(c); }