/**
 *	@file adtk_nodeset.cpp
 */

#include "adtk_nodeset.hpp"

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_nodeset::adtk_nodeset(){
	nodes.assign(1,0);
	tofs.assign(1,0);
}

adtk_nodeset::adtk_nodeset(const adtk_nodeset& n){
	nodes = n.nodes;
	tofs = n.tofs;
	velCont = n.velCont;
	constraints = n.constraints;
	nodeSize = n.nodeSize;
}

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

adtk_nodeset& adtk_nodeset::operator =(const adtk_nodeset &n){
	nodes = n.nodes;
	tofs = n.tofs;
	velCont = n.velCont;
	constraints = n.constraints;
	nodeSize = n.nodeSize;
	return *this;
}

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@param i the index of the node (begins at zero)
 *	@return the node specified by <tt>i</tt>
 */
std::vector<double> adtk_nodeset::getNode(int i){
	vector<double> node(nodes.begin()+i*nodeSize, nodes.begin()+(i+1)*nodeSize);
	return node;
}

/**
 *	@oaram i the index of the TOF; note that if there are n nodes, there will be n-1 TOFs.
 *	@return the TOF between nodes <tt>i</tt> and <tt>i+1</tt>
 */
double adtk_nodeset::getTOF(int i){ return tofs.at(i); }

/**
 *	@return the number of nodes contained in this node set
 */
double adtk_nodeset::getNumNodes(){ return nodes.size(); }

/**
 *	@return the number of states contained in one node
 */
int adtk_nodeset::getNodeSize(){ return nodeSize; }

/**
 *	Append a new node to the end of the nodes vector
 *	@param node a new node; should have <tt>nodeSize</tt> non-dimensional states
 */
void adtk_nodeset::appendNode(std::vector<double> node){
	nodes.insert(nodes.end(), node.begin(), node.end());
}

/**
 *	Append a new TOF to the end of the TOF vector
 *	@param tof a new TOF, non-dimensional units
 */
void adtk_nodeset::appendTOF(double tof){
	tofs.push_back(tof);
}

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------


