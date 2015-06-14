/**
 *	@file tpat_bcr4bpr_nodeset.cpp
 */
#include "tpat.hpp"
#include "tpat_bcr4bpr_nodeset.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

using namespace std;

/**
 *	@brief Construct a nodeset with no data other than the system
 *	@param data system data object describing the system the nodes exist in
 */
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(tpat_bcr4bpr_sys_data data) : tpat_nodeset(6){
	sysData = data;
}

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically set to splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(double IC[6], tpat_bcr4bpr_sys_data data, 
	double t0, double tof, int numNodes) : tpat_nodeset(6){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, tpat_nodeset::TIME);
	initEpochs(numNodes, t0);
}//======================================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically set to splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(std::vector<double> IC, tpat_bcr4bpr_sys_data data, 
	double t0, double tof, int numNodes) : tpat_nodeset(6){

	sysData = data;

	initSetFromICs(&(IC[0]), &sysData, t0, tof, numNodes, tpat_nodeset::TIME);
	initEpochs(numNodes, t0);
}//======================================================================

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
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(double IC[6], tpat_bcr4bpr_sys_data data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(6){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, type);
	initEpochs(numNodes, t0);
}//======================================================================

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
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(std::vector<double> IC, tpat_bcr4bpr_sys_data data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(6){

	sysData = data;

	initSetFromICs(&(IC[0]), &sysData, t0, tof, numNodes, type);
	initEpochs(numNodes, t0);
}//======================================================================

void tpat_bcr4bpr_nodeset::initEpochs(int numNodes, double t0){
	// Compute epoch times for each node
	epochs.reserve(numNodes);

	double ellapsed = t0;
	epochs.push_back(t0);
	for(int n = 0; n < numNodes-1; n++){
		ellapsed += tofs.at(n);
		epochs.push_back(ellapsed);
	}
}//=========================================

/**
 *	@brief Destructor
 */
tpat_bcr4bpr_nodeset::~tpat_bcr4bpr_nodeset(){
	epochs.clear();
}//=================================

/**
 *	@brief Copy constructor - calls base class copy constructor to handle basic copy
 *	@param n a BCR4BPR nodeset
 */
tpat_bcr4bpr_nodeset::tpat_bcr4bpr_nodeset(const tpat_bcr4bpr_nodeset &n) : tpat_nodeset(n){
	sysData = n.sysData;
	epochs = n.epochs;
}//=========================================

/**
 *	@brief Assignment operator - calls base class assignment operator to handle basic assignment
 *	@param n an input BCR4BPR nodeset
 *	@return this nodeset, made equal to the input nodeset
 */
tpat_bcr4bpr_nodeset& tpat_bcr4bpr_nodeset::operator =(const tpat_bcr4bpr_nodeset &n){
	tpat_nodeset::operator =(n);
	sysData = n.sysData;
	epochs = n.epochs;
	return *this;
}//=========================================

/**
 *	@brief Concatenate two nodesets
 *
 *	If the final node in <tt>lhs</tt> is the same as the first node in <tt>rhs</tt>, the
 *	concatenation will delete one occurence of the node to achieve continuity. Otherwise
 *	the nodes from <tt>rhs</tt> are concatenated to the end of <tt>lhs</tt>. The velocity
 *	continuity specifications and constraints for <tt>rhs</tt> are updated to reflect the 
 *	new indices of the nodes they describe.
 *
 *	@param lhs
 *	@param rhs
 *	@return a nodeset containing the concatenated input nodesets
 */
tpat_bcr4bpr_nodeset operator +(const tpat_bcr4bpr_nodeset &lhs, const tpat_bcr4bpr_nodeset &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot add nodesets from different systems; please transform them to be in the same system");
	}

	tpat_bcr4bpr_nodeset temp(lhs.sysData);
	tpat_bcr4bpr_nodeset::basicConcat(lhs, rhs, &temp);

	// Concatenate Epoch Times
	temp.epochs.clear();
	if(lhs.getNode(-1) == rhs.getNode(0)){
		temp.epochs.insert(temp.epochs.end(), lhs.epochs.begin(), lhs.epochs.end()-1);
	}else{
		temp.epochs.insert(temp.epochs.end(), lhs.epochs.begin(), lhs.epochs.end());
	}
	temp.epochs.insert(temp.epochs.end(), rhs.epochs.begin(), rhs.epochs.end());

	return temp;
}//=====================================================

/**
 *	@brief Retrieve a pointer to the vector of epochs
 *	@return a pointer to the beginning of the epochs vector
 */
std::vector<double>* tpat_bcr4bpr_nodeset::getEpochs(){ return &epochs; }

/**
 *	@brief Retrieve a specifi epoch
 *	@param i epoch index (begins with zero). If i is negative, the epoch will be
 *	selected from the end of the vector, i.e. -1 will give the last value, -2 will
 *	give the second to last, etc.
 *	@return the epoch, non-dimensional units
 */
double tpat_bcr4bpr_nodeset::getEpoch(int i) const {
	if(i < 0)
		i += epochs.size();
	return epochs.at(i);
}//=====================================

/**
 *	@brief  Retrieve a pointer to the system data object
 *	@return a pointer to the system data object for this nodeset
 */
tpat_sys_data* tpat_bcr4bpr_nodeset::getSysData() { return &sysData; }

/**
 *	@brief Add an epoch to the nodeset
 *	@param d an epoch (non-dimensional time) to add
 */
void tpat_bcr4bpr_nodeset::appendEpoch(double d){
	epochs.push_back(d);
}//=====================================

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void tpat_bcr4bpr_nodeset::print() const {
	printf("BCR4BPR Nodeset:\n  Nodes:\n");
	for(int n = 0; n < getNumNodes(); n++){
		printf("  > %02d -> @ %9.4f [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]", n,
			epochs[n], nodes[n*nodeSize+0], nodes[n*nodeSize+1], nodes[n*nodeSize+2], 
			nodes[n*nodeSize+3], nodes[n*nodeSize+4], nodes[n*nodeSize+5]);

		if(n < getNumNodes()-1){
			printf("  TOF = %.4f\n", tofs[n]);
		}else{
			printf("\n");
		}
	}
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}//========================================

void tpat_bcr4bpr_nodeset::reverseOrder(){
	
}

/**
 *	@brief Sort the nodes by their epoch time in chronological order
 */	
void tpat_bcr4bpr_nodeset::sortChrono(){
	// Get the indices to sort time
	vector<int> sortedIx = getSortedInd(epochs);

	vector<double> nodeCpy;
	vector<double> epochCpy;
	vector<double> tofCpy;

	for(int n = 0; n < ((int)epochs.size()); n++){
		cout << sortedIx[n] << endl;
		vector<double> node = getNode(sortedIx[n]);
		nodeCpy.insert(nodeCpy.end(), node.begin(), node.end());
		epochCpy.push_back(getEpoch(sortedIx[n]));
		
		if(n < (int)tofs.size())
			tofCpy.push_back(getTOF(sortedIx[n]));
	}

	nodes = nodeCpy;
	epochs = epochCpy;
	tofs = tofCpy;

	// Adjust nodes on constraints and velcon
	for(int n = 0; n < ((int)constraints.size()); n++){
		constraints[n].setNode(sortedIx[constraints[n].getNode()]);
	}

	// Adjust velcon nodes
	for(int n = 0; n < ((int)velConNodes.size()); n++){
		velConNodes[n] = sortedIx[velConNodes[n]];
	}
}//======================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_bcr4bpr_nodeset::saveToMat(const char* filename){
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
void tpat_bcr4bpr_nodeset::saveEpochs(mat_t *matFile){
	size_t dims[2] = {epochs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("Epochs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(epochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Epochs", MAT_COMPRESSION_NONE);
}//=========================================



