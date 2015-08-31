#include <cstdio>
#include <vector>
#include <stdexcept>

using namespace std;

class Shape{
protected:
	virtual void printMsg() = 0;
};

class Rectangle : public Shape{
public:
	Rectangle(){ printMsg(); }
protected:
	void printMsg() { printf("Rectangle\n"); }
};

class Square : public Rectangle{
public:
	Square() { printMsg(); }
protected:
	void printMsg() { printf("Square\n"); }
};

int main(){
	
	printf("Creating Rectangle:\n");	
	Rectangle r;

	printf("\nCreating Square:\n");
	Square s;

	return 0;
}