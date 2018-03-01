/**
 * Do a simple natural parameter continuation on low-thrust
 * L4/5 periodic orbits while maintaining a constant thrust magnitude
 * and vary the Jacobi Constant
 */
#include "astrohelion/AllIncludes.hpp"

using namespace astrohelion;

int main(){
	unsigned int lawID = ControlLaw_cr3bp_lt::CONST_FC_2D_LEFT;
	double f_in = 1e-2;		// Load family with this thrust magnitude
	double f_out = 2.5e-2;		// Compute family with this thrust magnitude

	// Load the L4 SPO family at the appropriate thrust level
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp-lt_earth-moon/L4_Lyap_f%.1e_Isp3000_law%u.mat", f_in, lawID);
	printf("Opening %s\n", filename);
	SysData_cr3bp_lt ltSys(filename);

	Family_PO spoFam(&ltSys);
	std::vector<ControlLaw*> loadedLaws {};
	spoFam.readFromMat(filename, loadedLaws);

	if(loadedLaws.empty())
		throw Exception("Did not load any control laws; cannot proceed");

	printf("Loaded %zu control laws:\n", loadedLaws.size());
	for(ControlLaw* p : loadedLaws)
		p->print();

	unsigned int ix0 = 1;

	Arcset_periodic member = spoFam.getMember(ix0);

	std::vector<Constraint> cons = member.getAllConstraints();
	member.clearAllConstraints();
	for(unsigned int i = 0; i < cons.size(); i++){
		if(cons[i].getType() != Constraint_tp::PSEUDOARC){
			member.addConstraint(cons[i]);
		}
	}
	member.print();
	
	// Change control law thrust value
	loadedLaws[0]->setParam(0, sqrt(f_out));

	PseudoArcEngine pae;
	std::vector<Arcset> allArcs;
	std::vector<int> initDir {1, 0, 0, 1, -1, 0};
	Arcset_cr3bp_lt a1(&ltSys), a2(&ltSys);

	pae.setNumOrbits(5);
	pae.pac(&member, &a1, &a2, allArcs, initDir);

	Family_PO fam(&ltSys);
	for(Arcset arc : allArcs)
		fam.addMember(static_cast<Arcset_cr3bp_lt>(arc));

	// fam.setSortType(FamSort_tp::SORT_TOF);
	// fam.sortMembers();
	fam.sortEigs();

	sprintf(filename, "../../data/families_toBeChecked/cr3bp-lt_L4Lyap_f%.1e_Isp3000_law%u.mat", f_out, lawID);
	fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}