#include "catch.hpp"

#include "math_util.h"

TEST_CASE( "Addition", "[nec_3vector]") {
    nec_3vector x(1,2,3);
    nec_3vector y(2,3,4);

    REQUIRE( x.norm2() == 14.0 );
    REQUIRE( x.norm() > 3.5 );

    nec_3vector r(2,4,6);
    REQUIRE( r == x*2 );

    nec_3vector s(3,5,7);
    REQUIRE( s == x+y);  
    REQUIRE( y == x+y-x);  
}
