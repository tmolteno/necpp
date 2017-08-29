#include "catch.hpp"

#include "nec_context.h"
#include "nec_results.h"
#include "c_geometry.h"

#include <iostream>

#define REQUIRE_APPROX_EQUAL(a, b) { \
  static nec_float eps = 3e-4; \
  REQUIRE(std::abs(a - b) < eps); }

TEST_CASE( "Example 1", "[example_1]") {

    /**
        CE EXAMPLE 1. CENTER FED LINEAR ANTENNA
        GW 0 7 0. 0. -.25 0. 0. .25 .001
        GE
        EX 0 0 4 0 1.
        XQ
        ld 0 0 4 4 10. 3.000E-09 5.300E-11
        PQ
        NE 0 1 1 15 .001 0 0 0. 0. .01786
        EN
    */
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(0, 7, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001, 1.0, 1.0);
    nec.geometry_complete(0);

    nec.ex_card(EXCITATION_VOLTAGE, 0, 4, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    nec.xq_card(0);

    nec.ld_card(0, 0, 4, 4, 10.0, 3.0e-9,5.3e-11);

    nec.pq_card(0,0,0,0);
    nec.ne_card(0, 1, 1, 15, 0.001, 0.0, 0.0, 0.0, 0.0, 0.01786);

/* TAG   SEG       VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)     POWER
  NO.   NO.     REAL      IMAGINARY     REAL      IMAGINARY     REAL      IMAGINARY    REAL       IMAGINARY   (WATTS)
   0     4  1.0000E+00  0.0000E+00  9.2058E-03 -5.1547E-03  8.2698E+01  4.6306E+01  9.2058E-03 -5.1547E-03  4.6029E-03
*/

    REQUIRE_APPROX_EQUAL(nec.get_impedance_real(0), 8.2698E+01 );
    REQUIRE_APPROX_EQUAL(nec.get_impedance_imag(0), 4.6306E+01 );

    nec_antenna_input* ai = nec.get_input_parameters(0);
    REQUIRE(ai->get_segment()[0] == 4 );
/*
                      ----- ANTENNA INPUT PARAMETERS -----
  TAG   SEG       VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)     POWER
  NO.   NO.     REAL      IMAGINARY     REAL      IMAGINARY     REAL      IMAGINARY    REAL       IMAGINARY   (WATTS)
   0     4  1.0000E+00  0.0000E+00  8.9547E-03 -4.0515E-03  9.2698E+01  4.1941E+01  8.9547E-03 -4.0515E-03  4.4773E-03
*/
    REQUIRE_APPROX_EQUAL(nec.get_impedance_real(1), 9.2698E+01 );
    REQUIRE_APPROX_EQUAL(nec.get_impedance_imag(1), 4.1941E+01 );
    
    ai = nec.get_input_parameters(1);
    REQUIRE(ai->get_segment()[0] == 4 );
    REQUIRE_APPROX_EQUAL(ai->get_current()[0], nec_complex(8.9547E-03, -4.0515E-03) );
}



TEST_CASE( "Voltage Excitation", "[voltage_excitation]") {

    /**
        CM A simple structure excited by a plane wave, and at multiple frequencies.
        CE go blue ! 
        GW 0 36 0 0 0 -0.042 0.008 0.017 0.001
        GE 0
        GN -1
        LD 5 0 0 0 3.720E+07 
        FR 0 1 0 0 2400
        PT -1
        EX 1 1 1 0 0 0 0 0 0 0 0
        RP 0 1 1 0500 90 90 0 0
        EN
    */
    nec_context nec;
    nec.initialize();

    int n_freq = 10;
    
    c_geometry* geo = nec.get_geometry();
    geo->wire(0, 36, 0.0, 0.0, 0.0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
    nec.geometry_complete(0);
    nec.gn_card(-1, 0, 0, 0, 0, 0, 0, 0);
    nec.ld_card(5, 0, 0, 0, 3.72e7, 0, 0);
    nec.fr_card(0, n_freq, 2400.0, 10.0);
    nec.pt_card(-1, 0, 0, 0);
    nec.ex_card(EXCITATION_VOLTAGE, 0, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    nec.rp_card(0, 10, 10, 0,5,0,0, 0.0, 0.0, 9.0, 9.0, 0.0, 0.0);
    //nec.xq_card(0);

    // Check that we get n_freq radiation patterns
    for (int freq_index = 0; freq_index < n_freq; freq_index++) {
        nec_radiation_pattern* rp = nec.get_radiation_pattern(freq_index);
        REQUIRE(rp->get_theta_start() == 0.0 );
//         REQUIRE_APPROX_EQUAL(rp->get_rp_power_average(), 12.0);
    }

}


TEST_CASE( "Plane Wave Excitation", "[plane_wave]") {

    /**
        CE A simple structure excited by a plane wave, and at multiple frequencies.
        GW 0 36 0 0 0 -0.042 0.008 0.017 0.001
        GE 0
        GN -1
        LD 5 0 0 0 3.720E+07 
        FR 0 10 0 0 2400.0 10.0
        PT -1
        EX 1 0 1 0 0 0 0 0 0 0 0
        RP 0 1 1 0500 90 90 0 0
        EN
    */
    nec_context nec;
    nec.initialize();

    int n_freq = 10;
    
    c_geometry* geo = nec.get_geometry();
    geo->wire(0, 36, 0.0, 0.0, 0.0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
    nec.geometry_complete(0);
    nec.gn_card(-1, 0, 0, 0, 0, 0, 0, 0);
    nec.ld_card(5, 0, 0, 0, 3.72e7, 0, 0);
    nec.fr_card(0, n_freq, 2400.0, 10.0);
    nec.pt_card(-1, 0, 0, 0);
    nec.ex_card(EXCITATION_LINEAR, 0, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    nec.rp_card(0, 10, 10, 0,5,0,0, 0.0, 0.0, 9.0, 9.0, 0.0, 0.0);
    //nec.xq_card(0);

    // Check that we get n_freq radiation patterns
    for (int freq_index = 0; freq_index < n_freq; freq_index++) {
        nec_radiation_pattern* rp = nec.get_radiation_pattern(freq_index);
        REQUIRE(rp->get_theta_start() == 0.0 );
//         REQUIRE_APPROX_EQUAL(rp->get_rp_power_average(), 12.0);
    }
}


