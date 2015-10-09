/**
 *  Integrate the grid of initial conditions from IceCube_ManTubeForwardProp.cpp
 *
 *	
 */

int main(){

	std::vector<double> gridICs = readMatrixFromMat("data/gridIC_C3.150.mat", "gridICs");
	
	return EXIT_SUCCESS;
}