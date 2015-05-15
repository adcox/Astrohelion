/**
*	Headers for BodyData object
*/
#ifndef __H_BODYDATA__
#define __H_BODYDATA__

#include <string>

class adtk_body_data{
	private:
		/** Mean radius of the body, km */
		double radius;

		/** Mass of the body, kg */
		double mass;

		/** Mean orbital radius (distance from parent body), km*/
		double orbitRad;

		/** Gravitational parameter associated with this body, kg^3/s^2 */
		double gravParam;

		/** Minmum acceptable fly-by altitude for this body. Altitudes lower than this will
		trigger a crash event in the numerical simulation */
		double minFlyByAlt;

		/** Name of this body */
		std::string name;

		/** Name of the parent body */
		std::string parent;

	public:

		// Constructors
		adtk_body_data();
		adtk_body_data(std::string);
		adtk_body_data(double m, double R, double r, double mu, std::string name, std::string parent);

		// Set and Get Functions
		double getRadius();
		double getMass();
		double getGravParam();
		double getOrbitRad();
		double getMinFlyBy();
		std::string getName();
		std::string getParent();

		void setRadius(double);
		void setMass(double);
		void setOrbitRad(double);
		void setGravParam(double);
		void setName(std::string);
		void setParent(std::string);
};

#endif
//END