/**
 *  @file tpat_model.cpp
 *
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat.hpp"

#include "tpat_model.hpp"

/**
 *	@brief Default constructor
 *	@param type the model type
 */
tpat_model::tpat_model(dynamicModel_t type){
	modelType = type;
}//===========================================

/**
 *	@brief Copy constructor
 */
tpat_model::tpat_model(const tpat_model &m){
	copyMe(m);
}//===========================================

/**
 *	@brief Destructor; virtual to allow derivative classes to fall back
 */
tpat_model::~tpat_model(){}

/**
 *	@brief Copy Operator
 */	
tpat_model& tpat_model::operator =(const tpat_model &m){
	copyMe(m);
	return *this;
}//============================================

/**
 *	@brief Copies all data from a dynamic model into this one
 *	@param m another dynamic model
 */
void tpat_model::copyMe(const tpat_model &m){
	modelType = m.modelType;
	coreStates = m.coreStates;
	stmStates = m.stmStates;
	extraStates = m.extraStates;
}//============================================

/**
 *	@brief Retrieve the number of core states
 *	@return the number of core states
 */
int tpat_model::getCoreStateSize() const { return coreStates; }

/**
 *	@brief Retrieve the number of STM elements stored in the state vector
 *	@return the number of STM elements stored in the state vector
 */
int tpat_model::getSTMStateSize() const { return stmStates; }

/**
 *	@brief Retrieve the number of extra states stored after the core states and STM elements
 *	@return the number of extra states stored after the core states and STM elements
 */
int tpat_model::getExtraStateSize() const { return extraStates; }



