#include <cstdio>
#include <vector>

using namespace std;

class Rectangle{
public:
	Rectangle(int w, int h) : width(w), height(h) {}
	int getArea() { return width * height; }
	int getPerimeter() {return 2*width + 2*height; }
private:
	int width;
	int height;
};

class Square : public Rectangle{
public:
	Square(int s) : Rectangle(s, s) {}
};

int main(){
	
	std::vector< std::vector<double> > data;

	data.push_back(std::vector<double>(0));

	printf("Data has %zu vectors stored\n", data.size());
	return 0;
}