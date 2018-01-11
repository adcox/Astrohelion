/**
 * Do a simple natural parameter continuation on low-thrust
 * L4/5 periodic orbits while maintaining a constant thrust magnitude
 * and vary the Jacobi Constant
 */
#include "astrohelion/AllIncludes.hpp"

using namespace astrohelion;

int main(){
	unsigned int lawID = ControlLaw_cr3bp_lt::CONST_F_C_2D_LEFT;
	double f = 1e-2;

	// Load the L3 family at the appropriate thrust level
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp-lt_earth-moon/L3_Lyap_f%.1e_Isp3000_law%u.mat", f, lawID);
	printf("Opening %s\n", filename);
	SysData_cr3bp_lt ltSys(filename);

	Family_PO lyapFam(&ltSys);
	std::vector<ControlLaw*> loadedLaws {};
	lyapFam.readFromMat(filename, loadedLaws);

	// printf("Loaded %zu control laws:\n", loadedLaws.size());
	// for(ControlLaw* p : loadedLaws)
	// 	p->print();

	unsigned int ix0 = 48;	// for lawID = 2, f = 1e-2

	Arcset_periodic member = lyapFam.getMember(ix0);
	member.clearAllConstraints();

	double L4pos[3];
	DynamicsModel_cr3bp::getEquilibPt(&ltSys, 4, 1e-6, L4pos);

	std::vector<double> yVal{lawID == 1 ? L4pos[1] : -L4pos[1]};
	Event planeEvt(Event_tp::XZ_PLANE, 0, true, yVal);

	printf("Stop at y = %f\n", yVal[0]);
	Arcset_cr3bp_lt temp(&ltSys), newMember(&ltSys);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.addEvent(planeEvt);

	sim.runSim(member.getStateByIx(0), member.getTotalTOF(), &temp, loadedLaws[0]);

	sim.clearEvents();
	sim.runSim_manyNodes(temp.getStateByIx(-1), member.getTotalTOF(), 5, &newMember, loadedLaws[0]);

	std::vector<double> initStateData {NAN, yVal[0], 0, NAN, NAN, 0, 1};
	std::vector<double> perData {0, 0, NAN, 0, NAN, NAN};

	Constraint periodicityCon(Constraint_tp::MATCH_CUST, newMember.getNodeByIx(-1).getID(), perData);
	Constraint fixInitState(Constraint_tp::STATE, newMember.getNodeByIx(0).getID(), initStateData);

	newMember.addConstraint(fixInitState);
	newMember.addConstraint(periodicityCon);

	newMember.saveToMat("data/tempInitGuess.mat");

	MultShootEngine shooter;
	Arcset_cr3bp_lt corrected(&ltSys);

	try{
		shooter.multShoot(&newMember, &corrected);
	}catch(Exception &e){
		printErr("%s\n", e.what());
	}

	corrected.saveToMat("data/tempCorrected.mat");

	PseudoArcEngine pae;
	std::vector<Arcset> allArcs;
	std::vector<int> initDir {1, 0, 0, 1, -1, 0};
	Arcset_cr3bp_lt a1(&ltSys), a2(&ltSys);

	pae.setNumOrbits(50);
	pae.pac(&newMember, &a1, &a2, allArcs, initDir);

	Family_PO fam(&ltSys);
	for(Arcset arc : allArcs)
		fam.addMember(static_cast<Arcset_cr3bp_lt>(arc));

	sprintf(filename, "../../data/families_toBeChecked/cr3bp-lt_L4Lyap_f%.1e_Isp3000_law%u.mat", f, lawID);
	fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}