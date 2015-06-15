#include <cstdio>
#include <vector>

using namespace std;

int main(void){
	
	int data[] = {1, 3, 5, 7, 9};
	vector<int> vec(data, data+5);

	printf("Vector: [%d %d %d %d %d]\n", vec[0], vec[1], vec[2], vec[3], vec[4] );

	vec.insert(vec.begin()+0, 8);
	
	printf("Vector: [%d %d %d %d %d %d]\n", vec[0], vec[1], vec[2], vec[3], vec[4], vec[5] );

	return 0;
}