/**
 *	@brief Test out the Eigen library of functions to make sure it will do what I want
 */
#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::JacobiSVD;
using Eigen::Vector3d;

MatrixXd getNullspace(MatrixXd A){
	JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullV);
	svd.setThreshold(1e-14);

	MatrixXd V = svd.matrixV();
	int rank = svd.rank();

	printf("rank = %d\n", rank);
	if(rank == V.cols())
		return MatrixXd(1,1);

	MatrixXd N((int)(V.rows()), (int)(V.cols()-rank));
	for(int r = 0; r < V.rows(); r++){
		for(int c = rank; c < V.cols(); c++){
			N(r,c-rank) = V(r,c);
		}
	}

	return N;
}

int main(){
	MatrixXd Z(2,3);
	Z<<	1, 0, 0,
		0, 0, 0;

  	std::cout << Z << std::endl;

	std::cout << "Here is the matrix Z:" << std::endl << Z<< std::endl;
	JacobiSVD<MatrixXd> svd(Z, Eigen::ComputeFullU | Eigen::ComputeFullV);
	std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
	std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
	std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;

	MatrixXd N = getNullspace(Z);
	std::cout << "Nullspace:\n" << N << std::endl;
}
