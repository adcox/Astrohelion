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
	
	Square square(4);
	Rectangle rect(2,3);

	Rectangle *cast_square = static_cast<Rectangle *>(&square);

	Rectangle *rectPtr;
	Rectangle *cast_try2 = static_cast<rectPtr>(&square);

	return 0;
}