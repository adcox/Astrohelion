/** 
 *	Header file for CR3BP System Data object
 */

#ifndef __H_CR3BP_SYS_DATA_
#define __H_CR3BP_SYS_DATA_

#include "adtk_sys_data.hpp"

class adtk_cr3bp_sys_data : public adtk_sys_data{
	private:
		/** Mass ratio between two primaries*/
		double mu;
		std::string P1;
		std::string P2;

	public:
		adtk_cr3bp_sys_data();
		adtk_cr3bp_sys_data(std::string P1, std::string P2);

		adtk_cr3bp_sys_data& operator=(const adtk_cr3bp_sys_data&);
		
		double getMu();
		std::string getPrimary(int n);	//We override this function, so re-declare it
};

#endif