#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

/**
 * @brief Free the memory containing the loaded control laws
 * 
 * @param laws reference to a vector of control law pointers that are stored
 * on the stack
 */
void freeLaws(std::vector<ControlLaw*>& laws){
	for(auto L : laws){
		delete L;
		L = nullptr;
	}
}//====================================================

int main(int argc, char** argv){
	(void) arc;
	(void) argv;

	char filename[128];

	// Load a family, propagate to y = L4, on the right-side of L4
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt fam(&sys);

	std::vector<ControlLaw*> laws {};
	fam.readFromMat(filename, laws);

	std::vector<double> yL4 {0.5*sqrt(3)};
	Event evt_yL4 = Event(Event_tp::XZ_PLANE, 0, true, yL4);

	int oID = 0;
	double tof = 0;
	std::vector<double> q0 {}, ctrl0 {};
	for(unsigned int m = 0; m < fam.getNumMembers(); m++){
		Arcset_Periodic mem = fam.getMember(m);
		mem.putInChronoOrder();

		for(unsigned int s = 0; s < mem.getNumSegs(s); s++){
			const Segment seg& = mem.getSegRefByIx_const(s);
			oID = seg.getOrigin();
			if(oID == Linkable::INVALID_ID){
				freeLaws(laws);
				throw Exception("Origin ID is INVALID_ID?!");
			}

			q0 = mem.getState(oID);
			if(seg.getCtrlLaw()){
				ctrl0 = mem.getNodeRef_const(oID).getExtraParamVec(PARAMKEY_CTRL);
			}

			
		}
	}
	return 0;
}