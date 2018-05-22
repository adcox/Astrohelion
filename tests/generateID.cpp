#include <iostream>

int main(){

	unsigned int law_gen, law_thrust, law_mass, law_options[6] = {0};
	printf("Welcome to the low-thrust control law discombobulator\n"
		"First, select the type of pointing strategy to use:\n"
		"  [0] None\n"
		"  [1] General pointing, fixed in rotating frame\n"
		"  [2] Velocity pointing\n"
		"  [3] Jacobi-preserving, 2D\n"
		"  [4] General pointing, fixed in inertial frame\n"
		"\nType the number of the pointing strategy: ");
	std::cin >> law_gen;

	if(law_gen == 0){
		printf("Law ID = %d\n", 0);
		return EXIT_SUCCESS;
	}

	printf("\nNext, select the thrust parameterization:\n"
		"  [0] Constant thrust\n"
		"  [1] Variable thrust, bounded: f = 0.5 f_max (sin(g) + 1)\n"
		"  [2] Variable thrust, unbounded: f = g^2\n"
		"\nType the number of the thrust strategy: ");
	std::cin >> law_thrust;

	printf("\nNow select the mass parameterization:\n"
		"  [0] Variable mass, constant specific impulse\n"
		"  [1] Constant mass\n"
		"\nType the number of the mass strategy:  ");
	std::cin >> law_mass;

	if(law_gen == 2){
		printf("\nFinally, select the specific pointing:\n"
			"  [0] Pro-velocity direction\n"
			"  [1] Anti-velocity direction\n"
			"\nType the number:  ");
		std::cin >> law_options[0];
	}else if(law_gen == 3){
		printf("\nFinally, select the specific pointing:\n"
			"  [0] Left\n"
			"  [1] Right\n"
			"\nType the number:  ");
		std::cin >> law_options[0];
	}

	// Concatenate the pieces to get a law ID
	unsigned int id = (law_gen << 11) | (law_thrust << 8) | (law_mass << 6);
	for(unsigned int i = 0; i < 6; i++){
		id |= (law_options[i] << (i-1));
	}

	printf("\nLaw ID = %d\n", id);
	return EXIT_SUCCESS;
}