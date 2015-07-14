#include <cstdio>
#include <vector>

using namespace std;

class Rectangle{
	public:
		Rectangle(double w, double h) : width(w), height(h) {}
		double getArea(){return width*height;}
		double getWidth(){return width;}
		double getHeight(){return height;}
	protected:
		double width = 0;
		double height = 0;
};

class Square : public Rectangle{
	public:
		Square(double s) : Rectangle(s,s) {}
};

double getPerim(Rectangle r){
	return 2*r.getWidth() + 2*r.getHeight();
}

int main(void){
	Rectangle r(3,3);
	Square s(2);

	Square s2 = static_cast<Square>(r);

	printf("Rectangle area: %f\n", r.getArea());
	printf("Square area: %f\n", s.getArea());

	printf("Rectangle perimeter: %f\n", getPerim(r));
	printf("Square perimeter: %f\n", getPerim(s));

	return 0;
}