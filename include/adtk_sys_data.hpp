/**
 * 	adtk_sys_data.hpp
 *
 * 	Header for System Data class
 */

#ifndef __H_SYSDATA_
#define __H_SYSDATA_

#include <string>

	class adtk_sys_data{

 		public:
 			enum system_t {UNDEF_SYS, CR3BP_SYS, BCR4BPR_SYS};

 			adtk_sys_data();

 			adtk_sys_data& operator= (const adtk_sys_data&);

 			virtual std::string getPrimary(int n) = 0;
 			double getCharL();
 			double getCharT();
 			double getCharM();
 			system_t getType();
 			std::string getTypeStr();

 			virtual double getMu() = 0;		// This one is virtual because ALL derivative classes must define it
 			virtual double getNu(){return 0;}	// These guys are only used for BCR4BP, so not virtual
 			virtual double getK(){return 0;}
 		protected:
 			/** Number of primaries that exist in this system */
 			int numPrimaries;

 			/** Characteristic quantities: length (km), time (s), mass (kg)*/
 			double charL, charT, charM;

 			system_t type;
 	};
#endif