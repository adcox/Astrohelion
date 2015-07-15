#include <cstdio>
#include <vector>

using namespace std;

enum blah_t{TYPEA = 1,
	TYPEB = 3,
	TYPEC = 77,
	TYPED = 4};

int main(void){
	
	blah_t blah1 = TYPEA;
	blah_t blah2 = TYPEB;
	blah_t blah3 = TYPEC;
	blah_t blah4 = TYPED;

	printf("Enum = %d\n", ((int)blah1));
	printf("Enum = %d\n", ((int)blah2));
	printf("Enum = %d\n", ((int)blah3));
	printf("Enum = %d\n", ((int)blah4));

	return 0;
}