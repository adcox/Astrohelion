/**
 *	Test out the trajectory data object
 */

#include "adtk_cr3bp_traj.hpp"
#include "adtk_cr3bp_sys_data.hpp"

#include <iostream>
#include <cstdio>

using namespace std;

int main(void){
	adtk_cr3bp_traj t;

	adtk_cr3bp_sys_data data("earth", "moon");
	cout << "Assigning sysData to traj " << data.getCharL() << endl;
	
	t.setSysData(data);

	adtk_cr3bp_sys_data data2 = t.getSysData();
	cout << data2.getCharL();
}