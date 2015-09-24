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
	// // Test with rows > cols, 1D nullspace
	// double Z_data[] = {1,0,0,0,0,0};
	// tpat_matrix Z(3,2,Z_data);
	// double N_data[] = {0,1};
	// tpat_matrix N_ans(2,1,N_data);
	// tpat_matrix Z_null = null_svd(Z);
	// bool test1 = Z_null == N_ans;

	// // Test with rows = cols, 2D nullspace
	// tpat_matrix W(2,2);
	// tpat_matrix N_ans2 = tpat_matrix::I(2);
	// bool test2 = null_svd(W) == N_ans2;
	
	// // Test with rows = cols, no nullspace
	// double X_data[] = {1,2,3,4};
	// tpat_matrix X(2,2,X_data);
	// tpat_matrix N = null_svd(X);
	// bool test3 = (N.getCols() == 1 && N.getRows() == 1 && N.at(0,0) == 0);

	// // Test with rows < cols, 2D nullspace
	// double Y_data[] = {1,0,0,0,0,0};
	// tpat_matrix Y(2,3, Y_data);
	// double N3_data[] = {0,0,1,0,0,1};
	// tpat_matrix N_ans3(3,2, N3_data);
	// bool test4 = null_svd(Y) == N_ans3;

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
