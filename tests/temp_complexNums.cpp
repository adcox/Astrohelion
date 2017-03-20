#include <cstdio>
#include <complex>
#include <string>

/**
 *  \brief Turn a complex number into a string, e.g. 1.2345 + 0.9876j
 *  \param num a complex number
 *  \return the complex number as a string
 */
std::string complexToStr(std::complex<double> num){
    char buffer[64];
    sprintf(buffer, "%.4e%s%.4ej", std::real(num), std::imag(num) < 0 ? " - " : " + ", std::abs(std::imag(num)));
    return std::string(buffer);
}

int main(){

	typedef std::complex<double> cdouble;	//!< A complex double

	cdouble a(1, 2);
	cdouble b = a*a;
	printf("%s\n", complexToStr(b).c_str());

	return EXIT_SUCCESS;
}