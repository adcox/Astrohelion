/**
 *	@file adtk_nodeset.cpp
 */

#include "adtk_nodeset.hpp"

#include <cstdio>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

adtk_nodeset::adtk_nodeset(const int n) : nodeSize(n){
	// Reserve space for one node and TOF
	nodes.reserve(n);
	tofs.reserve(1);
}//========================================================

adtk_nodeset::adtk_nodeset(const adtk_nodeset& n) : nodeSize(n.nodeSize){
	if(nodeSize == n.nodeSize){
		nodes = n.nodes;
		tofs = n.tofs;
	}else{
		fprintf(stderr, "Cannot create a nodeset with %d-size nodes from one with %d-size nodes\n",
			nodeSize, n.nodeSize);
		throw;
	}
}//========================================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

adtk_nodeset& adtk_nodeset::operator =(const adtk_nodeset &n){
	if(nodeSize == n.nodeSize){
		nodes = n.nodes;
		tofs = n.tofs;
		return *this;
	}else{
		fprintf(stderr, "Cannot create a nodeset with %d-size nodes from one with %d-size nodes\n",
			nodeSize, n.nodeSize);
		throw;
	}
}//======================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@param i the index of the node (begins at zero)
 *	@return the node specified by <tt>i</tt>
 */
std::vector<double> adtk_nodeset::getNode(int i) const {
	vector<double> node(nodes.begin()+i*nodeSize, nodes.begin()+(i+1)*nodeSize);
	return node;
}

/**
 *	@oaram i the index of the TOF; note that if there are n nodes, there will be n-1 TOFs.
 *	@return the TOF between nodes <tt>i</tt> and <tt>i+1</tt>
 */
double adtk_nodeset::getTOF(int i) const { return tofs.at(i); }

/**
 *	@return the number of nodes contained in this node set
 */
double adtk_nodeset::getNumNodes() const { return nodes.size(); }

/**
 *	@return the number of states contained in one node
 */
int adtk_nodeset::getNodeSize() const { return nodeSize; }

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


