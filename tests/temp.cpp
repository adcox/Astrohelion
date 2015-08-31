#include <cstdio>
#include <vector>
#include <stdexcept>

using namespace std;

class Rectangle{
public:
	Rectangle(){}
	Rectangle(int w, int h) : width(w), height(h) {}
	Rectangle(const Rectangle &r){ copyMe(r); }
	int getArea() { return width * height; }
	int getPerimeter() {return 2*width + 2*height; }
	int getWidth() const { return width; }
	int getHeight() const { return height; }
protected:
	void copyMe(const Rectangle &r){
		width = r.width;
		height = r.height;
	}

	int width = 0;
	int height = 0;
};

class Square : public Rectangle{
public:
	Square(int s) : Rectangle(s, s) {}
	Square(const Rectangle &s){ 
		if(s.getWidth() == s.getHeight())
			copyMe(s);
		else
			printf("Cannot create a square from a non-square rectangle!\n");
	}

	const Square& operator =(const Square &s){
		copyMe(s);
		return *this;
	}

	const Square& operator =(const Rectangle &s){
		if(s.getWidth() == s.getHeight()){
			copyMe(s);
			return *this;
		}else
			throw runtime_error("OOPS");
	}
};

int main(){
	
	Rectangle rect(2,3);
	Square sq = rect;
	Square sq2(4);
	sq2 = rect;

	return 0;
}