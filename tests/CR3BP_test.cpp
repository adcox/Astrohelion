/**
* Test the CR3BP numerical integrator
*/

#include <iostream>
#include <math.h>

#include "adtk_body_data.hpp"
#include "adtk_integrators.hpp"

using namespace std;

int main(void){
	adtk_body_data earth("earth");
	adtk_body_data moon("moon");

	double mu = moon.getMass()/(earth.getMass() + moon.getMass());
	cout << "Mu = " << mu << "\n";

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};
	double t[] = {0, 2.776481};
	vector<double> data = adtk_cr3bp_integrate(ic, t, mu, 2);

	double n = data.size()/7;
	cout << n << endl;

	double state[((int)n)][7];

	for(int r=0; r<n; r++){
		for (int c=0; c<7; c++){
			state[r][c] = data[r*7+c];
		}
	}

	for(int r=0; r<n; r++){
		printf("%6.12f : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", state[r][6], state[r][0],
			state[r][1], state[r][2], state[r][3], state[r][4], state[r][5]);
	}
}