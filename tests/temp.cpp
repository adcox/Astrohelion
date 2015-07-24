#include <cstdio>
#include <vector>

using namespace std;

typedef struct thing{
	int a;
	int b;
	int c;
} thing;

int main(){
	
	thing t;
	t.a = 1;
	t.b = 2;
	t.c = 3;

	switch(t.a){
		case 1:
			printf("yay!\n"); break;
		default:
			printf("boooo\n");
	}

	return 0;
}