#include <iostream>

using namespace std;

int main() {
    double arr[] = {1, 1.57283737283837, 2, 2.5, 3};

    double o[3];

    copy(arr+1, arr+4, o);

    cout << "The sub-array is: ";
    for(int i = 0; i < 3; i++){
    	cout << o[i] << ", ";
    }

    cout << endl;

    return 0;
}