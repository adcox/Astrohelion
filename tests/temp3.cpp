#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

void freeLaws(std::vector<ControlLaw*>&);
int main(int, char**);

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
	if(argc != 2)
		printErr("Expected 1 input but got %d\n", argc-1);

	char *famfile = argv[1];
	// char famfile[] = "L4_Lyap_f1.0e-01_Isp3000_law6144.mat";
	char filename[128];
	sprintf(filename, "%s/%s", "../../data/families/cr3bp-lt_earth-moon", famfile);

	// Load a family, propagate to y = L4, on the right-side of L4
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt oldFam(&sys), newFam(&sys);

	std::vector<ControlLaw*> laws {};
	oldFam.readFromMat(filename, laws);

	std::vector<double> yL4 {0.5*sqrt(3)};
	Event evt_yL4 = Event(Event_tp::XZ_PLANE, 0, true, yL4);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	
	unsigned int m = 0;

	int oID = 0;
	std::vector<double> q0 {}, ctrl0 {};
	bool foundEvt = false;
	
	printf("Reconverging member %03u\n", m);
	Arcset_periodic mem = oldFam.getMember(m);
	// std::vector<cdouble> vals = oldFam.getEigVals(m);
	// printf("Eigenvalues:\n");
	// for(unsigned int i = 0; i < vals.size(); i++){
	// 	std::cout << "  " << complexToStr(vals[i]) << std::endl;
	// }
	// mem.print();

	sim.addEvent(evt_yL4);

	for(unsigned int s = 0; s < mem.getNumSegs(); s++){
		const Segment& seg = mem.getSegRefByIx_const(s);
		oID = seg.getOrigin();
		if(oID == Linkable::INVALID_ID){
			freeLaws(laws);
			throw Exception("Origin ID is INVALID_ID?!");
		}

		q0 = mem.getState(oID);
		if(seg.getCtrlLaw() && seg.getCtrlLaw()->getNumStates() > 0){
			ctrl0 = mem.getNodeRef_const(oID).getExtraParamVec(PARAMKEY_CTRL);
		}

		Arcset_cr3bp_lt temp(&sys);
		sim.runSim(q0, ctrl0, 0, seg.getTOF(), &temp, seg.getCtrlLaw());

		if(temp.getNodeRefByIx(-1).getTriggerEvent() == evt_yL4.getType()){
			// This arc hit the event!
			q0 = temp.getStateByIx(-1);
			q0[6] = 1;
			foundEvt = true;
			break;
		}
	}

	// Re-converge orbit
	Arcset_cr3bp_lt ltpo(&sys);
	if(!foundEvt){
		freeLaws(laws);
		char msg[128];
		sprintf(msg, "Could not find event on member %u", m);
		throw Exception(msg);
	}else{
		MultShootEngine shooter;
		shooter.setVerbosity(Verbosity_tp::SOME_MSG);

		Arcset_cr3bp_lt temp(&sys);
		sim.clearEvents();
		sim.runSim_manyNodes(q0, ctrl0, 0, mem.getTotalTOF(),
			mem.getNumNodes(), &temp, mem.getSegRefByIx(0).getCtrlLaw());

		std::vector<Constraint> cons = mem.getAllConstraints();
		for(auto &c : cons){
			// if(c.getID() == mem.getNodeRefByIx(-1).getID()){
			// 	c.setID(temp.getNodeRefByIx(-1).getID());
			// }

			temp.addConstraint(c);
		}

		// temp.saveToMat("tempArc.mat");
		// waitForUser();
		// temp.print();
		// waitForUser();

		try{
			shooter.multShoot(&temp, &ltpo);
		}catch(const Exception &e){
			freeLaws(laws);
			throw e;
		}
	}
	
	PseudoArcEngine pae;
	std::vector<Arcset> allArcs;
	std::vector<int> initDir {-1, 0, 0, -1, 1, 0};
	Arcset_cr3bp_lt a1(&sys), a2(&sys);

	pae.setNumOrbits(500);
	pae.pac(&ltpo, &a1, &a2, allArcs, initDir);

	// Do the reverse direction
	std::vector<Arcset> otherArcs;
	for(int &val : initDir)
		val*= -1;

	pae.pac(&ltpo, &a1, &a2, otherArcs, initDir);

	// Add from otherArcs in reverse order
	for(int i = otherArcs.size()-1; i >= 0; i--)
		newFam.addMember(static_cast<Arcset_cr3bp_lt>(otherArcs[i]));

	// Add from allArcs in forward order
	for(Arcset arc : allArcs)
		newFam.addMember(static_cast<Arcset_cr3bp_lt>(arc));

	newFam.sortEigs();

	newFam.setName(oldFam.getName());
	newFam.saveToMat(famfile);

	freeLaws(laws);
	return 0;
}