/**
 *	I use this file to test out how dependencies work between super/derivative classes
 */

#include <iostream>
#include <string>

using namespace std;

class Shape{
	protected:
		string type;
	public:
		virtual double getArea() = 0;
		virtual string getType(){ return type; }
};

class Rectangle : public Shape{
	private:
		double width, height;
	public:
		Rectangle(double w, double h) : width(w), height(h) {
			type = "RECT";
		}
		double getArea(){ return width*height; }
		double getPerim(){ return 2*width + 2*height; }
		void setWidth(double w){ width = w; }
};

class Circle : public Shape{
	private:
		double radius;
	public:
		Circle(double r) : radius(r) {
			type = "CIRC";
		}
		double getArea(){ return 3.14159*radius*radius; }
		double getCirc(){ return 3.14159*2*radius; }
};

class DoStuff{
	private:
		Shape *shape;

	public:
		DoStuff(double d){
			if(d < 5)
				shape = new Rectangle(d, d);
			else
				shape = new Circle(d);
		}
		Shape* getShape(){return shape;}
};

int main(void){
	DoStuff doer(2);
	Shape *s = doer.getShape();
	string type = s->getType();

	cout << "Made a " << type << endl;
	cout << "Area: " << s->getArea() << endl;

	if(type.compare("RECT") == 0){
		// Convert Shape pointer to Rectangle pointer (since it really is a rectangle)
		Rectangle *rect = static_cast<Rectangle*>(s);
		// Copy the rectangle into a new object
		Rectangle rectObj = *rect;
		// Show the perimeter of the converted/copied shape
		cout << "Perimeter: " << rectObj.getPerim() << endl;
		// Change the width of the COPIED shape
		cout << ">> Set width of copied rectangle to 10" << endl;
		rectObj.setWidth(10);
		cout << "Copied Object Perimeter: " << rectObj.getPerim() << endl;
		// Note that the pointer still points to the old shape; its perimeter hasn't changed
		cout << "Original Pointer Perimeter: " << rect->getPerim() << endl;
		cout << ">> Set width of pointer rectangle to 1" << endl;
		rect->setWidth(1);
		cout << "Copied Object Perimeter: " << rectObj.getPerim() << endl;
		cout << "Original Pointer Perimeter: " << rect->getPerim() << endl;
	}else if(type.compare("CIRC") == 0){
		Circle *circ = static_cast<Circle*>(s);
		cout << "Circumference: " << circ->getCirc() << endl;
	}
	return 0;
}