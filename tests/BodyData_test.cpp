/**
* Test the BodyData object
*/

#include <iostream>
#include <math.h>

#include "tpat_body_data.hpp"

using namespace std;

int main(){
	tpat_body_data B("earth");

	cout << "Body: " << B.getName() << "\n";
	cout << "  Mass: " << B.getMass() << " kg\n";
	cout << "  Radius: " << B.getRadius() << " km\n";
	cout << "  Parent: " << B.getParent() << "\n";
	cout << "  Gravitational Parameter: " << B.getGravParam() << " km^3/s^2\n";
	cout << "  Minimum Fly-by Altitude: " << B.getMinFlyBy() << " km\n";
}