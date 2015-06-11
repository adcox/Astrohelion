/**
 *	Test precision in C++
 */
#include <iostream>
#include <cstdlib>

using namespace std;

int main(void){
	double val = 1.0;
	double x = 1;
	double y = x/2 + val;
	int count = -1;

	while(y > val){
		x /= 2;
		y = x/2 + val;
		count--;
	}

	cout << "Precision: 2^"<<count<<endl;
	return 0;
}