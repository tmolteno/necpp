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


TEST_CASE( "Left-handed helix starts on +x axis", "[helix_handedness]") {
    // Regression test for #91: LH helix must start at (A1, 0, 0)
    // on the +x axis, not at (0, A1) on the +y axis.
    // Check coordinates before geometry_complete (which rescales
    // meters to wavelengths).
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    // Left-handed helix: HL < 0, tag=1, 10 segments, spacing=0.1,
    // length=-0.3, a1=0.2, b1=0.1, a2=0.2, b2=0.1, radius=0.001
    geo->helix(1, 10, 0.1, -0.3, 0.2, 0.1, 0.2, 0.1, 0.001);

    // First segment should start on +x axis: x ≈ a1, y ≈ 0
    REQUIRE(geo->x[0] == Approx(0.2).margin(1e-6));
    REQUIRE(geo->y[0] == Approx(0.0).margin(1e-6));
}

TEST_CASE( "Flat spiral does not divide by zero", "[helix_flat_spiral]") {
    // Regression test for #91: HL=0 with A2 != A1 used to
    // divide by fabs(hl), causing a divide-by-zero.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    // Flat spiral: hl=0, a1=0.1, a2=0.2 (different radii)
    geo->helix(1, 10, 0.1, 0.0, 0.1, 0.1, 0.2, 0.2, 0.001);
    REQUIRE_NOTHROW(nec.geometry_complete(0));

    // Should have 10 segments
    REQUIRE(geo->n_segments == 10);
}

TEST_CASE( "Connected wires do not trigger false intersection", "[false_intersection]") {
    // Regression test for #87: two wires sharing an endpoint
    // must not trigger a "SEGMENT MIDPOINT INTERSECTS" error.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    // Wire 1: from origin to (0, 0, 0.5)
    geo->wire(1, 11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.002, 1.0, 1.0);
    // Wire 2: from (0, 0, 0.5) to (0.36, 0.36, 0.5)
    // Shares endpoint with wire 1 — this used to cause a false
    // positive intersection error.
    geo->wire(2, 11, 0.0, 0.0, 0.5, 0.36, 0.36, 0.5, 0.002, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));
}
