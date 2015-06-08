/**
 *	@file adtk_bcr4bpr_nodeset.cpp
 */

#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_utilities.hpp"

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
adtk_bcr4bpr_nodeset::adtk_bcr4bpr_nodeset(double IC[6], adtk_bcr4bpr_sys_data data, 
		double t0, double tof, int numNodes, node_distro_t type) : adtk_nodeset(6){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, type);

	// Compute epoch times for each node
	epochs.reserve(numNodes);

	double ellapsed = t0;
	epochs.push_back(t0);
	for(int n = 0; n < numNodes-1; n++){
		ellapsed += tofs.at(n);
		epochs.push_back(ellapsed);
	}
}//======================================================================

/**
 *	@brief Destructor
 */
adtk_bcr4bpr_nodeset::~adtk_bcr4bpr_nodeset(){
	epochs.clear();
}//=================================

/**
 *	@brief Copy constructor - calls base class copy constructor to handle basic copy
 *	@param n a BCR4BPR nodeset
 */
adtk_bcr4bpr_nodeset::adtk_bcr4bpr_nodeset(const adtk_bcr4bpr_nodeset &n) : adtk_nodeset(n){
	sysData = n.sysData;
	epochs = n.epochs;
}//=========================================

/**
 *	@brief Assignment operator - calls base class assignment operator to handle basic assignment
 *	@param n an input BCR4BPR nodeset
 *	@return this nodeset, made equal to the input nodeset
 */
adtk_bcr4bpr_nodeset& adtk_bcr4bpr_nodeset::operator =(const adtk_bcr4bpr_nodeset &n){
	adtk_nodeset::operator =(n);
	sysData = n.sysData;
	epochs = n.epochs;
	return *this;
}//=========================================

/**
 *	@brief Retrieve a pointer to the vector of epochs
 *	@return a pointer to the beginning of the epochs vector
 */
std::vector<double>* adtk_bcr4bpr_nodeset::getEpochs(){ return &epochs; }

/**
 *	@brief Retrieve a specifi epoch
 *	@param i epoch index (begins with zero)
 *	@return the epoch, non-dimensional units
 */
double adtk_bcr4bpr_nodeset::getEpoch(int i) const {
	return epochs.at(i);
}//=====================================

/**
 *	@brief  Retrieve a pointer to the system data object
 *	@return a pointer to the system data object for this nodeset
 */
adtk_sys_data* adtk_bcr4bpr_nodeset::getSysData() { return &sysData; }

/**
 *	@brief Add an epoch to the nodeset
 *	@param d an epoch (non-dimensional time) to add
 */
void adtk_bcr4bpr_nodeset::appendEpoch(double d){
	epochs.push_back(d);
}//=====================================

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void adtk_bcr4bpr_nodeset::print() const {
	printf("BCR4BPR Nodeset:\n  Nodes:\n");
	for(int n = 0; n < getNumNodes(); n++){
		printf("  > %02d -> [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]\n", n,
			nodes[n*nodeSize+0], nodes[n*nodeSize+1], nodes[n*nodeSize+2], 
			nodes[n*nodeSize+3], nodes[n*nodeSize+4], nodes[n*nodeSize+5]);
	}
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}//========================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void adtk_bcr4bpr_nodeset::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file @version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		saveNodes(matfp);
		saveTOFs(matfp);
		saveEpochs(matfp);
		sysData.saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the epoch values to a file
 *	@param matFile a pointer to the destination matlab file
 */
void adtk_bcr4bpr_nodeset::saveEpochs(mat_t *matFile){
	size_t dims[2] = {epochs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("Epochs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(epochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Epochs", MAT_COMPRESSION_NONE);
}//=========================================



