#include "tpat/tpat_utilities.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main(void){
	
	double data[] = {1,3,2,5,4};
	vector<double> vec(data, data+5);

	vector<int> ix = getSortedInd(vec);

	for(int i = 0; i < ((int)ix.size()); i++){
		cout << ix[i] << endl;
	}

	return 0;
}