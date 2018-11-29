/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <cstring>
#include <mex.h>
#include <vector>

#include "MexStreamOutput.hpp"

#include "AllIncludes.hpp"

using namespace astrohelion;

namespace astroHelper{

/* *****************************************************************************
 * 	Function Declarations
 **************************************************************************** */


SysData_tp getSysType(const mxArray*, unsigned int ix = 0);

// Converting Astrohelion objects to mex data
void arcset2MexData(const Arcset*, const char**, unsigned int, mxArray*,
	unsigned int structIx = 0);
void cons2MexData(const Arcset*, mxArray*);		// TODO
void linkTable2MexData(const Arcset*, mxArray*);
void nodeCtrlState2MexData(const Arcset*, mxArray*);
void nodeState2MexData(const Arcset*, mxArray*);
void nodeStateDeriv2MexData(const Arcset*, mxArray*);	// TODO
void nodeTime2MexData(const Arcset*, mxArray*);
void segCtrl2MexData(const Arcset*, mxArray*);
void segState2MexData(const Arcset*, mxArray*);
void segTime2MexData(const Arcset*, mxArray*);
void TOF2MexData(const Arcset*, mxArray*);
void STM2MexData(const Arcset*, mxArray*);

// Converting mex data to Astrohelion objects
void mexData2Arcset(const mxArray*, Arcset*, std::vector<ControlLaw*>&,
	unsigned int ix = 0);
void mexData2Cons(const mxArray*, Arcset*);		// TODO
void mexData2Ctrl(const mxArray*, unsigned int, unsigned int&, std::vector<double>&);
void mexData2Event(const mxArray*, std::vector<Event>&);
void mexData2LinkTable(const mxArray*, Arcset*);
void mexData2NodeState(const mxArray*, Arcset*);
void mexData2NodeStateDeriv(const mxArray*, Arcset*);	// TODO
void mexData2NodeTime(const mxArray*, Arcset*);	
void mexData2NodeCtrl(const mxArray*, Arcset*);
void mexData2SegCtrl(const mxArray*, Arcset*, std::vector<ControlLaw*>&);
void mexData2SegState(const mxArray*, Arcset*);
void mexData2SegTime(const mxArray*, Arcset*);
void mexData2TOF(const mxArray*, Arcset*);
void mexData2STM(const mxArray*, Arcset*);

// Converting mex data to SysData objects
SysData_cr3bp mexData2SysData_cr3bp(const mxArray*, unsigned int ix = 0);
SysData_cr3bp_lt mexData2SysData_cr3bp_lt(const mxArray*, unsigned int ix = 0);

const char* CTRL_FIELDS[] = {VARNAME_CTRL_TP, VARNAME_CTRL_NSTATE, VARNAME_CTRL_PARAM};

/* *****************************************************************************
 * 	Function Definitions: Utility
 **************************************************************************** */

/**
 * @brief Retrieve the system type of the arcset stored in the mxArray
 * @details [long description]
 * 
 * @param pMxArray pointer to the mxArray that represents the structure array
 * @param ix index of the arcset within the structure array (by default ix = 0)
 * 
 * @return The system type
 */
SysData_tp getSysType(const mxArray *pMxArray, unsigned int ix){
	// Super-hacked together way to figure out what kind of arcset is loaded
	mxArray *pP1 = mxGetField(pMxArray, ix, VARNAME_P1);
	mxArray *pP2 = mxGetField(pMxArray, ix, VARNAME_P2);
	mxArray *pP3 = mxGetField(pMxArray, ix, VARNAME_P3);
	mxArray *pRefMass = mxGetField(pMxArray, ix, VARNAME_M0);

	if(pP1 != nullptr && pP2 == nullptr && pP3 == nullptr){
		return SysData_tp::R2BP_SYS;
	}else if(pP1 != nullptr && pP2 != nullptr && pP3 == nullptr && pRefMass == nullptr){
		return SysData_tp::CR3BP_SYS;
	}else if(pP1 != nullptr && pP2 != nullptr && pP3 != nullptr && pRefMass == nullptr){
		return SysData_tp::BCR4BPR_SYS;
	}else if(pP1 != nullptr && pP2 != nullptr && pP3 == nullptr && pRefMass != nullptr){
		return SysData_tp::CR3BP_LT_SYS;
	}else{
		return SysData_tp::UNDEF_SYS;
	}
}//====================================================

/* *****************************************************************************
 * 	Function Definitions: Convert TO Mex Data
 **************************************************************************** */

/**
 * @brief Convert an arcset to MATLAB mex data
 * 
 * @param pArc pointer to an arcset to transform
 * @param fieldnames array of fieldnames
 * @param nNames number of names in the fieldnames array
 * @param pStructArray pointer to a struct array that the arcset data will be 
 * stored in
 * @param structIx (Optional) index of the structure within the pStructArray 
 * array. By default, this value is set to zero
 * 
 * @throws mexErrMsgTxt if the arcset cannot be converted into the correct data
 * structures or if the requested fieldnames cannot be populated.
 */
void arcset2MexData(const Arcset* pArc, const char **fieldnames, 
	unsigned int nNames, mxArray *pStructArray, unsigned int structIx){

	const SysData* pSys = pArc->getSysData();
	const unsigned int nNodes = pArc->getNumNodes();
	const unsigned int nSegs = pArc->getNumSegs();
	const unsigned int core_dim = pSys->getDynamicsModel()->getCoreStateSize();

	for(unsigned int f = 0; f < nNames; f++){
		// Create a temporary mxArray; MUST be instantiated in this loop!
		// However, different fields require different types
		mxArray *pTempArray;

		// Fill the pTempArray array with the data; each subfunction resizes and 
		// reallocates the memory within pTempArray
		if(astroHelper::strcmpi(fieldnames[f], VARNAME_LINKTABLE)){
			mwSize dims[] = {nSegs,4};
			pTempArray = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
			linkTable2MexData(pArc, pTempArray);
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_NODESTATE)){
			// Allocate the data
			pTempArray = mxCreateDoubleMatrix(nNodes, core_dim, mxREAL);
			//  Populate the data
			nodeState2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_NODETIME)){
			pTempArray = mxCreateDoubleMatrix(nNodes, 1, mxREAL);
			nodeTime2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_NODECTRL)){
			pTempArray = mxCreateCellMatrix(nNodes, 1);
			nodeCtrlState2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_SEGSTATE)){
			pTempArray = mxCreateCellMatrix(nSegs, 1);
			segState2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_SEGTIME)){
			pTempArray = mxCreateCellMatrix(nSegs, 1);
			segTime2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_TOF)){
			pTempArray = mxCreateDoubleMatrix(nSegs, 1, mxREAL);
			TOF2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_STM)){
			pTempArray = mxCreateCellMatrix(nSegs, 1);
			STM2MexData(pArc, pTempArray);
		}
		else if(astroHelper::strcmpi(fieldnames[f], VARNAME_SEGCTRL)){
			
			pTempArray = mxCreateStructMatrix(nSegs, 1, 3, CTRL_FIELDS);
			segCtrl2MexData(pArc, pTempArray);
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_P1)){
			pTempArray = mxCreateString(pSys->getPrimary(0).c_str());
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_P2)){
			if(pSys->getType() != SysData_tp::R2BP_SYS){
				pTempArray = mxCreateString(pSys->getPrimary(1).c_str());
			}else{
				char msg[64];
				snprintf(msg, 64, "%s does not have field %s\n", 
					pSys->getTypeStr().c_str(), fieldnames[f]);
				mexErrMsgTxt(msg);
			}
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_P3)){
			if(pSys->getType() == SysData_tp::BCR4BPR_SYS){
				pTempArray = mxCreateString(pSys->getPrimary(2).c_str());
			}else{
				char msg[64];
				snprintf(msg, 64, "%s does not have field %s\n", 
					pSys->getTypeStr().c_str(), fieldnames[f]);
				mexErrMsgTxt(msg);
			}
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_MU)){
			if(pSys->getType() == SysData_tp::CR3BP_SYS || pSys->getType() == SysData_tp::CR3BP_LT_SYS){
				const SysData_cr3bp *pTempSys = static_cast<const SysData_cr3bp*>(pSys);
				pTempArray = mxCreateDoubleScalar(pTempSys->getMu());	
			}else{
				char msg[64];
				snprintf(msg, 64, "%s does not have field %s\n", 
					pSys->getTypeStr().c_str(), fieldnames[f]);
				mexErrMsgTxt(msg);
			}
		}else if(astroHelper::strcmpi(fieldnames[f], VARNAME_M0)){
			if(pSys->getType() == SysData_tp::CR3BP_LT_SYS){
				const SysData_cr3bp_lt *pTempSys = static_cast<const SysData_cr3bp_lt*>(pSys);
				pTempArray = mxCreateDoubleScalar(pTempSys->getRefMass());
			}else{
				char msg[64];
				snprintf(msg, 64, "%s does not have field %s\n", 
					pSys->getTypeStr().c_str(), fieldnames[f]);
				mexErrMsgTxt(msg);
			}
		}else{
			char msg[64];
			snprintf(msg, 64, "%s field is not handled; this may crash MATLAB\n",
				fieldnames[f]);
			mexErrMsgTxt(msg);
		}
		
		// Set the output field
		mxSetField(pStructArray, structIx, fieldnames[f], pTempArray);
	}
}//====================================================

/**
 * @brief Convert an Arcset link table to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void linkTable2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);
	std::vector<int> linkTable{};
	pArc->getLinkTable(linkTable);

	// Copy data from link table to output
	int *pData = static_cast<int*>(mxGetData(pMxArray));
	std::copy(linkTable.begin(), linkTable.end(), pData);
}//====================================================

/**
 * @brief Convert an Arcset node state to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void nodeState2MexData(const Arcset* pArc, mxArray *pMxArray){
	std::vector<double> data {};
	pArc->getNodeStates(data);
	std::copy(data.begin(), data.end(), mxGetPr(pMxArray));
}//====================================================

/**
 * @brief Convert an Arcset node time to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void nodeTime2MexData(const Arcset* pArc, mxArray *pMxArray){
	std::vector<double> epochs {};
	pArc->getNodeEpochs(epochs);
	std::copy(epochs.begin(), epochs.end(), mxGetPr(pMxArray));
}//====================================================

/**
 * @brief Convert an Arcset node control state to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void nodeCtrlState2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nNodes = mxGetM(pMxArray);
	for(unsigned int n = 0; n < nNodes; n++){
		// Get control vector, if it exists
		std::vector<double> ctrl {};
		try{
			ctrl = pArc->getNodeRefByIx(n).getExtraParamVec(PARAMKEY_CTRL);
		}catch(const Exception &e){}

		// create double matrix
		mxArray *pCtrlArray;
		if(ctrl.empty()){
			pCtrlArray = mxCreateDoubleMatrix(0, 0, mxREAL);
		}else{
			pCtrlArray = mxCreateDoubleMatrix(1, ctrl.size(), mxREAL);
			double *pData = mxGetPr(pCtrlArray);
			std::copy(ctrl.begin(), ctrl.end(), pData);
		}

		mxSetCell(pMxArray, n, pCtrlArray);
	}
}//====================================================

/**
 * @brief Convert an Arcset segment state to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void segState2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);
	for(unsigned int s = 0; s < nSegs; s++){
		std::vector<double> state = pArc->getSegRefByIx(s).getStateVector();
		unsigned int w = pArc->getSegRefByIx(s).getStateWidth();
		unsigned int h = state.size()/w;

		// Copy data from row-major storage to column-major storage
		mxArray *pStateArray = mxCreateDoubleMatrix(h, w, mxREAL);
		double *pData = mxGetPr(pStateArray);
		for(unsigned int r = 0; r < h; r++){
			for(unsigned int c = 0; c < w; c++){
				pData[c*h + r] = state[r*w + c];
			}
		}

		mxSetCell(pMxArray, s, pStateArray);
	}
}//====================================================

/**
 * @brief Convert an Arcset segment time to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void segTime2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);
	for(unsigned int s = 0; s < nSegs; s++){
		std::vector<double> t = pArc->getSegRefByIx(s).getTimeVector();
		mxArray *pTimeArray = mxCreateDoubleMatrix(t.size(), 1, mxREAL);
		double *pData = mxGetPr(pTimeArray);
		std::copy(t.begin(), t.end(), pData);
		mxSetCell(pMxArray, s, pTimeArray);
	}
}//====================================================

/**
 * @brief Convert an Arcset segment control state to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void segCtrl2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);

	for(unsigned int s = 0; s < nSegs; s++){
		const ControlLaw* law = pArc->getCtrlLawByIx(s);
		
		if(law){
			// Save the law type
			mwSize dims[] = {1,1};
			mxArray *pType = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
			unsigned int* pData = static_cast<unsigned int*>(mxGetData(pType));
			pData[0] = law->getType();
			mxSetField(pMxArray, s, CTRL_FIELDS[0], pType);

			// Save the number of states
			mxArray *pNStates = mxCreateDoubleScalar(law->getNumStates());
			mxSetField(pMxArray, s, CTRL_FIELDS[1], pNStates);

			// Save the parameters
			std::vector<double> params = law->getParamsRef();
			mxArray *pParams;
			if(params.empty())
				pParams = mxCreateDoubleMatrix(0, 0, mxREAL);
			else{
				pParams = mxCreateDoubleMatrix(1, params.size(), mxREAL);
				double *pParamData = mxGetPr(pParams);
			std::copy(params.begin(), params.end(), pParamData);
			}

			mxSetField(pMxArray, s, CTRL_FIELDS[2], pParams);
		}else{
			// Fill struct with empty matrices for each field
			mxArray *pArray = mxCreateDoubleMatrix(0,0, mxREAL);
			mxSetField(pMxArray, s, CTRL_FIELDS[0], pArray);
			mxArray *pArray2 = mxDuplicateArray(pArray);
			mxSetField(pMxArray, s, CTRL_FIELDS[1], pArray2);
			mxArray *pArray3 = mxDuplicateArray(pArray);
			mxSetField(pMxArray, s, CTRL_FIELDS[2], pArray3);

		}
	}
}//====================================================

/**
 * @brief Convert an Arcset times of flight to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void TOF2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);
	double *pData = mxGetPr(pMxArray);

	for(unsigned int s = 0; s < nSegs; s++){
		pData[s] = pArc->getTOFByIx(s);
	}
}//====================================================

/**
 * @brief Convert an Arcset state transition matrices to mex data
 * 
 * @param pArc Arcset to extract the data from
 * @param pMxArray MATLAB data pointer to store the data in
 */
void STM2MexData(const Arcset* pArc, mxArray *pMxArray){
	const unsigned int nSegs = mxGetM(pMxArray);

	for(unsigned int s = 0; s < nSegs; s++){
		std::vector<double> stm = pArc->getSTMElementsByIx(s);	// row-major order
		unsigned int w = std::sqrt(stm.size());

		mxArray *pSTMArray = mxCreateDoubleMatrix(w, w, mxREAL);
		double *pData = mxGetPr(pSTMArray);

		for(unsigned int r = 0; r < w; r++){
			for(unsigned int c = 0; c < w; c++){
				pData[c*w + r] = stm[r*w + c];
			}
		}

		mxSetCell(pMxArray, s, pSTMArray);
	}
}//====================================================

/* *****************************************************************************
 * 	Function Definitions: Convert FROM Mex Data
 **************************************************************************** */

/**
 * @brief Populate an arcset with data from the MATLAB mex data
 * @details [long description]
 * 
 * @param pMxArray Pointer to struct array that includes the data necessary to 
 * construct the arcset
 * @param pArc an arcset that has been constructed with the proper
 * SystemData object
 * @param ix If the struct array includes more than one element, ix specifies 
 * the index of the struct within the array. By default ix = 0.
 * 
 * @see mexData2SysData_cr3bp, mexData2SysData_cr3bp_lt
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2Arcset(const mxArray *pMxArray, Arcset *pArc, 
	std::vector<ControlLaw*>& refLaws, unsigned int ix){

	if(!mxIsStruct(pMxArray))
		return;

	const int nFields = mxGetNumberOfFields(pMxArray);

	// First task: must load the link table
	if(const mxArray *pLinkTable = mxGetField(pMxArray, ix, VARNAME_LINKTABLE)){
		mexData2LinkTable(pLinkTable, pArc);
	}else{
		mexErrMsgTxt("mexData2Arcset: No link table, cannot initialize arc\n");
	}

	// Parse Node information
	if(const mxArray *pNodeState = mxGetField(pMxArray, ix, VARNAME_NODESTATE))
		mexData2NodeState(pNodeState, pArc);
	if(const mxArray *pNodeStateDeriv = mxGetField(pMxArray, ix, VARNAME_STATE_DERIV))
		mexData2NodeStateDeriv(pNodeStateDeriv, pArc);
	if(const mxArray *pNodeTime = mxGetField(pMxArray, ix, VARNAME_NODETIME))
		mexData2NodeTime(pNodeTime, pArc);
	if(const mxArray *pNodeCtrl = mxGetField(pMxArray, ix, VARNAME_NODECTRL))
		mexData2NodeCtrl(pNodeCtrl, pArc);

	// Parse Segment information
	if(const mxArray *pSegCtrl = mxGetField(pMxArray, ix, VARNAME_SEGCTRL))
		mexData2SegCtrl(pSegCtrl, pArc, refLaws);
	if(const mxArray *pSegState = mxGetField(pMxArray, ix, VARNAME_SEGSTATE))
		mexData2SegState(pSegState, pArc);
	if(const mxArray *pSegTime = mxGetField(pMxArray, ix, VARNAME_SEGTIME))
		mexData2SegTime(pSegTime, pArc);
	if(const mxArray *pTOF = mxGetField(pMxArray, ix, VARNAME_TOF))
		mexData2TOF(pTOF, pArc);
	if(const mxArray *pSTM = mxGetField(pMxArray, ix, VARNAME_STM))
		mexData2STM(pSTM, pArc);

	// Read/Parse constraint information
	if(const mxArray *pCons = mxGetField(pMxArray, ix, VARNAME_CONSTRAINTS))
		mexData2Cons(pCons, pArc);
}//====================================================

/**
 * @brief Load constraints from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores constraints
 * @param pArc pointer to Arcset object that will be populated with Consraint 
 * data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2Cons(const mxArray *pMxArray, Arcset *pArc){
	mexPrintf("mexData2Cons: Not implemented yet\n");
}//====================================================

/**
 * @brief Load the link table from mex data and it to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores the link table. Data is 
 * expected to be an integer matrix with one row for each segment
 * @param pArc pointer to Arcset object that will be initialized from the link
 * table. NOTE: This function MUST be called before populating the arcset with 
 * other data (e.g., before populating node states, segment states, etc.)
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2LinkTable(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsNumeric(pMxArray))
		mexErrMsgTxt("mexData2LinkTable: expects numeric data\n");

	if(mxGetN(pMxArray) != 4)
		mexErrMsgTxt("mexData2LinkTable: expects 4 columns\n");

	const unsigned int nSegs = mxGetM(pMxArray);
	const int* pData = static_cast<const int*>(mxGetData(pMxArray));
	pArc->initFromLinktable(pData, nSegs);
}//====================================================

/**
 * @brief Load node states from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores node states. Data is expected
 * to be a matrix with one row for each node state
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2NodeState(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsDouble(pMxArray))
		mexErrMsgTxt("mexData2NodeState: expects double matrix\n");

	const unsigned int nNodes = mxGetM(pMxArray);
	const unsigned int nStates = mxGetN(pMxArray);

	if(nNodes != pArc->getNumNodes()){
		pArc->print();
		mexErrMsgTxt("mexData2NodeState: node state rows and number of nodes"
			" initialized from the link table are not consistent\n");
	}
	if(nStates != pArc->getSysData()->getDynamicsModel()->getCoreStateSize()){
		mexErrMsgTxt("mexData2NodeState: state size (# cols) is inconsistent"
			" dynamics model core size\n");
	}

	const double *pData = mxGetPr(pMxArray);
	pArc->setNodeStates(pData, nNodes, nStates);
}//====================================================

/**
 * @brief Load node state derivatives from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores node state derivatives. 
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2NodeStateDeriv(const mxArray *pMxArray, Arcset *pArc){
	mexPrintf("mexData2NodeStateDeriv: Not implemented yet\n");
}//====================================================

/**
 * @brief Load node times (i.e., epochs) from mex data and add them to the 
 * Arcset object
 * 
 * @param pMxArray pointer to mex data that stores node times. Data is expected
 * to be a vector with one element for each node epoch
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2NodeTime(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsDouble(pMxArray))
		mexErrMsgTxt("mexData2NodeTime: expects double matrix\n");

	if(mxGetN(pMxArray) != 1)
		mexErrMsgTxt("mexData2NodeTime: expect column vector\n");

	const unsigned int nNodes = mxGetM(pMxArray);

	if(nNodes != pArc->getNumNodes()){
		mexErrMsgTxt("mexData2NodeTime: # node time rows is inconsistent"
			" with # nodes created from the link table\n");
	}

	const double *pData = mxGetPr(pMxArray);
	pArc->setNodeEpochs(pData, nNodes);
}//====================================================

/**
 * @brief Load node control states from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores node control states. Data is 
 * expected to be a cell array with one cell for each node control vector
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2NodeCtrl(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsCell(pMxArray))
		mexErrMsgTxt("mexData2NodeCtrl: expects cell array\n");

	if(mxGetN(pMxArray) != 1)
		mexErrMsgTxt("mexData2NodeCtrl: expect column array\n");

	const unsigned int nNodes = mxGetM(pMxArray);

	if(nNodes != pArc->getNumNodes()){
		mexErrMsgTxt("mexData2NodeCtrl: # node control rows is inconsistent with"
			" the # of nodes created from the link table\n");
	}

	for(unsigned int n = 0; n < nNodes; n++){
		mxArray *pCell = mxGetCell(pMxArray, n);

		if(!mxIsDouble(pCell))
			mexErrMsgTxt("mexData2NodeCtrl: expect double values in cell contents");

		const unsigned int r = mxGetM(pCell);
		const unsigned int c = mxGetN(pCell);
		
		if(r*c > 0){
			const double *pData = mxGetPr(pCell);
			std::vector<double> ctrl(pData, pData + r*c);
			pArc->getNodeRefByIx(n).setExtraParamVec(PARAMKEY_CTRL, ctrl);
		}
	}
}//====================================================

/**
 * @brief Load segment control laws from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores segment control. Data is expected
 * to be a structure array with one structure for each segment control law
 * @param pArc pointer to Arcset object that will be populated with data
 * @refLaws stores pointers to ControlLaws that are loaded from the data. These
 * laws are allocated on the stack and must be deleted before the program ends
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2SegCtrl(const mxArray *pMxArray, Arcset *pArc,
	std::vector<ControlLaw*> &refLaws){
	
	if(!mxIsStruct(pMxArray))
		mexErrMsgTxt("mexData2SegCtrl: expects struct matrix\n");

	if(mxGetM(pMxArray) != 1 && mxGetN(pMxArray) != 1)
		mexPrintf("segment structure array is not 1D; may yield unexpected results\n");

	const unsigned int nSegs = mxGetN(pMxArray) * mxGetM(pMxArray);

	if(nSegs != pArc->getNumSegs()){
		mexErrMsgTxt("mexData2SegCtrl: # seg control elements is inconsistent"
			" with # of segs created from the link table\n");
	}

	for(unsigned int s = 0; s < nSegs; s++){
		unsigned int tp = 0;
		std::vector<double> params {};

		mexData2Ctrl(pMxArray, s, tp, params);
		pArc->setSegCtrl(tp, params, s, refLaws);
	}
}//====================================================

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param pMxArray poitner to mex data that stores segment control. Data is 
 * expected to be a structure array with one structure for each control law
 * @param ix index of the control law structure within the larger array
 * @param tp reference to a storage variable for the parsed control law type
 * @param params reference to a storage vector for the parsed parameter vector
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2Ctrl(const mxArray *pMxArray, unsigned int ix, unsigned int &tp,
		std::vector<double> &params){

	if(~params.empty())
		params.clear();

	if(const mxArray *pType = mxGetField(pMxArray, ix, VARNAME_CTRL_TP)){
		if(!mxIsNumeric(pType))
			mexErrMsgTxt("mexData2Ctrl: expects type to be numeric\n");
		tp = *static_cast<const int*>(mxGetData(pType));
	}else{
		mexErrMsgTxt("mexData2Ctrl: seg ctrl does not specify type\n");
	}

	if(const mxArray *pParam = mxGetField(pMxArray, ix, VARNAME_CTRL_PARAM)){
		if(!mxIsDouble(pParam))
			mexErrMsgTxt("mexData2Ctrl: expect control parameters as doubles\n");
		const double *pData = mxGetPr(pParam);
		params.insert(params.begin(), pData, pData + mxGetM(pParam)*mxGetN(pParam));
	}else{
		mexErrMsgTxt("mexData2Ctrl: seg ctrl does not specify parameters\n");
	}
}//====================================================

/**
 * @brief Load segment states from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores segment states. Data is 
 * expected to be a cell array, with one cell for each segment state matrix
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2SegState(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsCell(pMxArray))
		mexErrMsgTxt("mexData2SegState: expects cell array\n");

	if(mxGetM(pMxArray) != 1 && mxGetN(pMxArray) != 1)
		mexPrintf("segment state cell array is not 1D; may yield unexpected results\n");

	const unsigned int nSegs = mxGetN(pMxArray) * mxGetM(pMxArray);

	if(nSegs != pArc->getNumSegs()){
		mexErrMsgTxt("mexData2SegState: # seg state elements is inconsistent"
			" with # of segs created from the link table\n");
	}

	for(unsigned int s = 0; s < nSegs; s++){
		mxArray *pCell = mxGetCell(pMxArray, s);
		if(!mxIsDouble(pCell))
			mexErrMsgTxt("mexData2SegState: expects double marix\n");

		const double *pData = mxGetPr(pCell);
		pArc->setSegState(s, pData, mxGetM(pCell), mxGetN(pCell));
	}
}//====================================================

/**
 * @brief Load segment times from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores segment times. Data is
 * expected to be a cell array with one cell for each vector of segment times.
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2SegTime(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsCell(pMxArray))
		mexErrMsgTxt("mexData2SegTime: expects cell array\n");

	if(mxGetM(pMxArray) != 1 && mxGetN(pMxArray) != 1)
		mexPrintf("segment time cell array is not 1D; may yield unexpected results\n");

	const unsigned int nSegs = mxGetN(pMxArray) * mxGetM(pMxArray);

	if(nSegs != pArc->getNumSegs()){
		mexErrMsgTxt("mexData2SegTime: # seg time elements is inconsistent"
			" with # of segs created from the link table\n");
	}

	for(unsigned int s = 0; s < nSegs; s++){
		mxArray *pCell = mxGetCell(pMxArray, s);
		if(!mxIsDouble(pCell))
			mexErrMsgTxt("mexData2SegTime: expects double marix\n");
		const double *pData = mxGetPr(pCell);
		std::vector<double> time(pData, pData + mxGetM(pCell)*mxGetN(pCell));
		pArc->getSegRefByIx(s).setTimeVector(time);
	}
}//====================================================

/**
 * @brief Load times-of-flight from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores times-of-flight. Data is
 * expected to be stored in a double array
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2TOF(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsDouble(pMxArray))
		mexErrMsgTxt("mexData2TOF: expects double array\n");

	if(mxGetM(pMxArray) != 1 && mxGetN(pMxArray) != 1)
		mexPrintf("TOF array is not 1D; may yield unexpected results\n");

	const unsigned int nSegs = mxGetN(pMxArray) * mxGetM(pMxArray);

	if(nSegs != pArc->getNumSegs()){
		mexErrMsgTxt("mexData2TOF: # TOF elements is inconsistent"
			" with # of segs created from the link table\n");
	}

	const double *pData = mxGetPr(pMxArray);
	for(unsigned int s = 0; s < nSegs; s++){
		pArc->getSegRefByIx(s).setTOF(pData[s]);
	}
}//====================================================

/**
 * @brief Load STMs from mex data and add them to the Arcset object
 * 
 * @param pMxArray pointer to mex data that stores STMs. Data is expected to be
 * stored in a cell array, with one cell for each STM
 * @param pArc pointer to Arcset object that will be populated with data
 * 
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2STM(const mxArray *pMxArray, Arcset *pArc){
	if(!mxIsCell(pMxArray))
		mexErrMsgTxt("mexData2STM: expects cell array\n");

	if(mxGetM(pMxArray) != 1 && mxGetN(pMxArray) != 1)
		mexPrintf("segment STM cell array is not 1D; may yield unexpected results\n");

	const unsigned int nSegs = mxGetN(pMxArray) * mxGetM(pMxArray);

	if(nSegs != pArc->getNumSegs()){
		mexErrMsgTxt("mexData2STM: # seg time elements is inconsistent"
			" with # of segs created from the link table\n");
	}

	unsigned int s = 0, r = 0, c = 0;
	for(s = 0; s < nSegs; s++){
		const mxArray *pCell = mxGetCell(pMxArray, s);
		if(!mxIsDouble(pCell))
			mexErrMsgTxt("mexData2STM: expects double marix\n");

		const double *pData = mxGetPr(pCell);
		unsigned int rows = mxGetM(pCell), cols = mxGetN(pCell);
		std::vector<double> stm(rows*cols, 0);

		for(r = 0; r < rows; r++){
			for(c = 0; c < cols; c++){
				stm[r*cols + c] = pData[c*rows + r];
			}
		}

		pArc->getSegRefByIx(s).setSTM(stm);
	}
}//====================================================

/**
 * @brief Convert MATLAB mex data to a set of events
 * 
 * @param pMxArray Pointer to a struct array that contains the events. The struct
 * array should have 1 column and n rows; each struct should include at least
 * four fields: the event type, trigger direction, stop propagation flag, and
 * parameter array. Optionally, the stop count may also be specified.
 * @param events reference to a storage vector for the events parsed from pMxArray
 * 
 * @see astrohelion::Event for explanations of the fields
 * @see Common.hpp for required field names (case insensitive)
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
void mexData2Event(const mxArray *pMxArray, std::vector<Event>& events){
	if(!mxIsStruct(pMxArray)){
		mexErrMsgTxt("mexData2Event: Expecting structure or structure array\n");
		return;
	}

	const int nFields = mxGetNumberOfFields(pMxArray);
	if(nFields < 4){
		mexErrMsgTxt("mexData2Event: Expecting at least 4 fields for each event\n");
		return;
	}
	
	const int nRows = mxGetM(pMxArray);
	const int nCols = mxGetN(pMxArray);
	
	if(nRows == 0 || nCols == 0)
		return;

	if(nCols != 1 && nRows != 1){
		mexErrMsgTxt("mexData2Event: Expecting 1 column or 1 row for event struct array\n");
		return;
	}

	for(unsigned int r = 0; r < nRows*nCols; r++){
		// Initialize emtpy data values
		Event_tp tp = Event_tp::NONE;
		int dir = 0, stopCount = -1;
		bool willStop = true;
		std::vector<double> params {};

		for(unsigned int f = 0; f < nFields; f++){
			const char* fName = mxGetFieldNameByNumber(pMxArray, f);
			const mxArray *pField = mxGetField(pMxArray, r, fName);

			if(astroHelper::strcmpi(fName, VARNAME_EVT_TP)){
				if(mxIsDouble(pField)){
					const double *pData = mxGetPr(pField);
					tp = static_cast<Event_tp>(pData[0]);
				}else{
					mexErrMsgTxt("mexData2Event: Expect numeric value for tp\n");
				}
			}else if(astroHelper::strcmpi(fName, VARNAME_EVT_DIR)){
				if(mxIsDouble(pField)){
					const double *pData = mxGetPr(pField);
					if(mxGetM(pField) == 0 || mxGetN(pField) == 0)
						mexErrMsgTxt("mexData2Event: trigger direction is empty\n");
					dir = static_cast<int>(pData[0]);
				}else{
					mexErrMsgTxt("mexData2Event: Expect numeric value for dir\n");
				}
			}else if(astroHelper::strcmpi(fName, VARNAME_EVT_STOP)){
				if(mxIsLogicalScalar(pField)){
					willStop = mxIsLogicalScalarTrue(pField);
				}else{
					mexErrMsgTxt("mexData2Event: Expect event StopProp to be a scalar boolean\n");
				}
			}else if(astroHelper::strcmpi(fName, VARNAME_EVT_PARAM)){
				if(mxIsDouble(pField)){
					const double *pData = mxGetPr(pField);
					if(mxGetM(pField) != 0 || mxGetN(pField) != 0){
						params.insert(params.begin(), pData,
							pData + mxGetM(pField)*mxGetN(pField));
					}
				}else{
					mexErrMsgTxt("mexData2Event: Expect event params to store doubles\n");
				}
			}else if(astroHelper::strcmpi(fName, VARNAME_EVT_STOPCOUNT)){
				if(mxIsDouble(pField)){
					const double *pData = mxGetPr(pField);
					if(mxGetM(pField) != 0 || mxGetN(pField) != 0)
						stopCount = static_cast<int>(pData[0]);
				}else{
					mexErrMsgTxt("mexData2Event: Expect numeric value for trigger direction\n");
				}
			}
		}

		// Create event and add it to the output array
		Event evt(tp, dir, willStop, params);
		if(stopCount > 0)
			evt.setStopCount(stopCount);
		events.push_back(evt);
	}
}//====================================================

/**
 * @brief Construct a CR3BP system data object from MATLAB mex data
 * 
 * @param pMxArray Pointer to struct array that includes the data necessary to 
 * construct the system data (P1, P2)
 * @param ix If the struct array includes more than one element, ix specifies 
 * the index of the struct within the array. By default ix = 0.
 * 
 * @return A system data object
 * @throws mexErrMsgTxt if the data cannot be parsed.
 */
SysData_cr3bp mexData2SysData_cr3bp(const mxArray *pMxArray, unsigned int ix){
	if(!mxIsStruct(pMxArray))
		mexErrMsgTxt("mexData2SysData_cr3bp: data array is not a structure");

	const int nFields = mxGetNumberOfFields(pMxArray);
	
	if(mxGetM(pMxArray) == 0 || mxGetN(pMxArray) == 0)
		mexErrMsgTxt("mexData2SysData_cr3bp: data array has zero elements");

	std::string P1 = "", P2 = "";
	bool bP1 = false, bP2 = false;
	for(unsigned int f = 0; f < nFields; f++){
		const char* fName = mxGetFieldNameByNumber(pMxArray, f);
		const mxArray *pField = mxGetField(pMxArray, ix, fName);

		if(astroHelper::strcmpi(fName, VARNAME_P1)){
			P1 = std::string(mxArrayToString(pField));
			bP1 = true;
		}else if(astroHelper::strcmpi(fName, VARNAME_P2)){
			P2 = std::string(mxArrayToString(pField));
			bP2 = true;
		}

		if(bP1 && bP2)
			break;
	}

	if(!bP1 || !bP2)
		mexErrMsgTxt("mexData2SysData_cr3bp: Could not load both P1 and P2 from mex data\n");
	
	// Create system
	return SysData_cr3bp(P1.c_str(), P2.c_str());
}//====================================================

/**
 * @brief Construct a low-thrust CR3BP system data object from MATLAB mex data
 * 
 * @param pMxArray Pointer to struct array that includes the data necessary to 
 * construct the system data (P1, P2, initial mass)
 * @param ix If the struct array includes more than one element, ix specifies 
 * the index of the struct within the array. By default ix = 0.
 * 
 * @return A system data object
 * @throws mexErrMsgTxt if the data cannot be read properly
 */
SysData_cr3bp_lt mexData2SysData_cr3bp_lt(const mxArray *pMxArray, unsigned int ix){
	if(!mxIsStruct(pMxArray))
		mexErrMsgTxt("mexData2SysData_cr3bp_lt: data array is not a structure");

	const int nFields = mxGetNumberOfFields(pMxArray);
	
	if(mxGetM(pMxArray) == 0 || mxGetN(pMxArray) == 0)
		mexErrMsgTxt("mexData2SysData_cr3bp_lt: data array has zero elements");

	std::string P1 = "", P2 = "";
	double M0 = 0;
	bool bP1 = false, bP2 = false, bM0 = false;
	for(unsigned int f = 0; f < nFields; f++){
		const char* fName = mxGetFieldNameByNumber(pMxArray, f);
		const mxArray *pField = mxGetField(pMxArray, ix, fName);

		if(astroHelper::strcmpi(fName, VARNAME_P1)){
			P1 = std::string(mxArrayToString(pField));
			bP1 = true;
		}else if(astroHelper::strcmpi(fName, VARNAME_P2)){
			P2 = std::string(mxArrayToString(pField));
			bP2 = true;
		}else if(astroHelper::strcmpi(fName, VARNAME_M0)){
			const double *pData = mxGetPr(pField);
			
			if(mxGetM(pField) == 0 || mxGetN(pField) == 0)
				mexErrMsgTxt("mexData2SysData_cr3bp_lt: reference mass field is empty\n");

			M0 = pData[0];
			bM0 = true;
		}

		if(bP1 && bP2 && bM0)
			break;
	}

	if(!bP1 || !bP2 || !bM0){
		mexErrMsgTxt("mexData2SysData_cr3bp_lt: Could not load P1, P2, "
			"and RefMass from mex data\n");
	}
	
	// Create system
	return SysData_cr3bp_lt(P1.c_str(), P2.c_str(), M0);
}

} // End of astroHelper namespace 