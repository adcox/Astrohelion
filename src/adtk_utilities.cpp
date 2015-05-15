
#include "adtk_utilities.hpp"

#include <vector>

/**
 *	"Transpose" a vector from row-major order to column-major order, or vise versa.
 *	Warning: operation occurs in-place!
 *	@param v an input vector
 *	@return the transposed vector
 */
template<T> vector<T> adtk_utilities::transposeVec(vector<T> v){
 	if(v.size() <= 1){
 		return v;
 	}
	for(int n = 1;)
}