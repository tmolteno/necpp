#include "catch.hpp"

#include "c_evlcom.h"

#include <cstdlib>
#include <iostream>

void test_bessel(nec_float r, nec_float i, nec_float j0r, nec_float j0i) {
     nec_complex z(r,i);
     
     nec_complex j0, j0p;
     
     bessel( z, &j0, &j0p );
     
     nec_float eps = std::fabs(j0r*1.0e-4);

     REQUIRE((std::fabs(j0.real() -  j0r) < eps));
     REQUIRE((std::fabs(j0.imag() -  j0i) < eps));
}

/** Generated with octave:
 * besselj(0, 1j, 0)
 * besselj(0, 1 + 1j, 0)
 * besselj(0, 10j, 0)
 * besselj(0, 10 + 10j, 0)
 */
TEST_CASE( "Bessel Functions", "[bessel]") {
     test_bessel(0, 1, 1.2661, 0);
     test_bessel(1, 1, 0.93761, -0.49653);
     test_bessel(0, 10, 2815.71, 0);
     test_bessel(10, 10, -2314.98, 411.56);
}
