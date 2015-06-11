#include <iostream>
using namespace std;

class Animal{
public:
	Animal(int m) : mass(m){}
	int getMass() const { return mass; }
protected:
	int mass;
	static void eat(const Animal& lhs, const Animal& rhs, Animal *result){
		result->mass = lhs.getMass() + rhs.getMass();
	}
};

class Gazelle : public Animal{
public:
	Gazelle(int m) : Animal(m){}
};

class Lion : public Animal{
public:
	Lion(int m) : Animal(m){}

	friend Lion feed(const Lion &lhs, const Gazelle &rhs){
		Lion hungry(0);
		eat(lhs, rhs, &hungry);
		return hungry;
	}
};

int main(void){
	Lion leo(5);
	Gazelle greg(1);

	Lion fullLeo = feed(leo, greg);
	cout << "Full Leo has mass " << fullLeo.getMass() << endl;
}