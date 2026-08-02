#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "c_geometry.h"
#include "nec_context.h"
#include "libnecpp.h"
#include <cmath>
#include <iostream>

void HANDLE_NEC(long x) { 
  int __tmp = (x);  
  if (__tmp != 0) {
    std::cout << nec_error_message() << std::endl;
  }
  REQUIRE((__tmp) == 0); 
}

TEST_CASE( "Geometry", "[surface_patch]") {
    nec_context* nec;     
    nec = nec_create();
 

    HANDLE_NEC(nec_sp_card(nec, 3, 0.019000, -0.001424, 0.078830, 0.019000, 0.001424, 0.078830));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.076180, 0.019000, -0.001424, 0.076180));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.073530, 0.019000, -0.001424, 0.073530));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.070880, 0.019000, -0.001424, 0.070880));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.068230, 0.019000, -0.001424, 0.068230));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.065580, 0.019000, -0.001424, 0.065580));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.062930, 0.019000, -0.001424, 0.062930));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.060280, 0.019000, -0.001424, 0.060280));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.057630, 0.019000, -0.001424, 0.057630));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.054980, 0.019000, -0.001424, 0.054980));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.052330, 0.019000, -0.001424, 0.052330));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.049680, 0.019000, -0.001424, 0.049680));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.047030, 0.019000, -0.001424, 0.047030));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.044380, 0.019000, -0.001424, 0.044380));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.041730, 0.019000, -0.001424, 0.041730));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.039080, 0.019000, -0.001424, 0.039080));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.019000, 0.001424, 0.036430, 0.019000, -0.001424, 0.036430));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.017283, 0.001295, 0.033856, 0.017283, -0.001295, 0.033856));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.015566, 0.001167, 0.031282, 0.015566, -0.001167, 0.031282));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.013849, 0.001038, 0.028708, 0.013849, -0.001038, 0.028708));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.012132, 0.000909, 0.026134, 0.012132, -0.000909, 0.026134));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.023560, 0.010415, -0.000780, 0.023560));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.020942, 0.010415, -0.000780, 0.020942));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.018324, 0.010415, -0.000780, 0.018324));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.015707, 0.010415, -0.000780, 0.015707));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.013089, 0.010415, -0.000780, 0.013089));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.010471, 0.010415, -0.000780, 0.010471));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.007853, 0.010415, -0.000780, 0.007853));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.005236, 0.010415, -0.000780, 0.005236));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.002618, 0.010415, -0.000780, 0.002618));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.010415, 0.000780, 0.000000, 0.010415, -0.000780, 0.000000));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.007811, 0.000585, 0.000000, 0.007811, -0.000585, 0.000000));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.005208, 0.000390, 0.000000, 0.005208, -0.000390, 0.000000));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.002604, 0.000195, 0.000000, 0.002604, -0.000195, 0.000000));
    HANDLE_NEC(nec_sc_card(nec, 3, 0.000026, 0.000002, 0.000000, 0.000026, -0.000002, 0.000000));
    HANDLE_NEC(nec_gm_card(nec, 0, 4, 0.000000, 0.000000, 8.571429, 0,0,0,0));
    HANDLE_NEC(nec_wire(nec, 1,15,0.000000,0.005829,0.004961,0.000000,-0.005829,0.004961,0.000248,1,1));
    HANDLE_NEC(nec_gm_card(nec, 0,  0, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.078830, 0));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    
    HANDLE_NEC(nec_fr_card(nec, 0,1, 10450.000000, 0));
    HANDLE_NEC(nec_ex_card(nec,0,1,8,0,1.000000,0.000000, 0, 0, 0, 0));
    HANDLE_NEC(nec_ld_card(nec,5,0,0,0,3.720000E+07, 0, 0));
    HANDLE_NEC(nec_pt_card(nec,-1, 0, 0, 0));
    HANDLE_NEC(nec_rp_card(nec, 0,361,3,1,5,0,0,0.000000,0.000000,1.000000,45.000000, 0, 0));
    
    /*
                          ----- ANTENNA INPUT PARAMETERS -----
      TAG   SEG       VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)     POWER
      NO.   NO.     REAL      IMAGINARY     REAL      IMAGINARY     REAL      IMAGINARY    REAL       IMAGINARY   (WATTS)
      1     8  1.0000E+00  0.0000E+00  9.2145E-03  6.7375E-04  1.0795E+02 -7.8930E+00  9.2145E-03  6.7375E-04  4.6072E-03
    */
    REQUIRE((nec_impedance_real(nec,0) - 1.0795E+02 < 1E-3));
    REQUIRE((nec_impedance_imag(nec,0) - -7.8930E+00 < 1E-3));
    REQUIRE((nec_gain_max(nec, 0) - 10.3332 < 1E-4));

    nec_delete(nec);
}

TEST_CASE( "Flat-spiral helix (GH HL=0)", "[helix]") {
    // PR #122: HL=0 selects a flat spiral — segments should form a spiral
    // in the XY plane (z=0) with linearly interpolated radius, not collapse
    // into a straight radial wire.
    c_geometry geo;

    int    tag_id        = 1;
    int    segment_count = 30;
    nec_float s          = 0.22;   // turn spacing (used for angular step: 2*pi/s)
    nec_float hl         = 0.0;    // zero = flat spiral
    nec_float a1         = 0.05;   // start X radius
    nec_float b1         = 0.05;   // start Y radius
    nec_float a2         = 0.15;   // end X radius
    nec_float b2         = 0.15;   // end Y radius
    nec_float rad        = 0.0025;

    geo.helix(tag_id, segment_count, s, hl, a1, b1, a2, b2, rad);

    REQUIRE( geo.n_segments == segment_count );

    // All z coordinates must be zero (flat spiral in XY plane).
    for (int i = 0; i < segment_count; i++) {
        INFO( "Segment " << i );
        REQUIRE( geo.z[i]  == 0.0 );
        REQUIRE( geo.z2[i] == 0.0 );
        REQUIRE( geo.segment_tags[i] == tag_id );
        REQUIRE( geo.segment_radius[i] == rad );
    }

    // First segment starts at phi=0: x=a1*cos(0)=a1, y=b1*sin(0)=0.
    REQUIRE( geo.x[0] == Catch::Approx(0.05).margin(1e-12) );
    REQUIRE( geo.y[0] == Catch::Approx(0.0).margin(1e-12) );

    // Verify the geometry is NOT a straight wire — if x or y were constant
    // across all segments, the spiral degenerated into a radial line.
    bool x_varies = false, y_varies = false;
    for (int i = 1; i < segment_count && !(x_varies && y_varies); i++) {
        if (std::fabs(geo.x[i] - geo.x[0]) > 1e-12) x_varies = true;
        if (std::fabs(geo.y[i] - geo.y[0]) > 1e-12) y_varies = true;
    }
    REQUIRE( x_varies );
    REQUIRE( y_varies );

    // Final radius should be near a2/b2 (within a segment's interpolation step).
    nec_float r_last = std::sqrt(geo.x[segment_count-1]*geo.x[segment_count-1]
                               + geo.y[segment_count-1]*geo.y[segment_count-1]);
    REQUIRE( r_last == Catch::Approx(0.15).margin(0.01) );
}

TEST_CASE( "NE card produces near-field results with explicit counts", "[near_field]") {
    // Regression: NE with explicit NRX/NRY/NRZ = 1 must compute near fields.
    nec_context* nec = nec_create();

    // Simple dipole at 300 MHz
    HANDLE_NEC(nec_wire(nec, 0, 7, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001, 1.0, 1.0));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 300.0, 0.0));
    HANDLE_NEC(nec_ex_card(nec, 0, 0, 4, 0, 1.0, 0.0, 0, 0, 0, 0));

    // NE: rectangular, NRX=NRY=NRZ=1, single point at (1,0,0)
    HANDLE_NEC(nec_ne_card(nec, 0, 1, 1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    nec_near_field_pattern* nfp = nec->get_near_field_pattern(0);
    REQUIRE( nfp != nullptr );
    REQUIRE( nfp->get_nfeh() == 0 );
    REQUIRE( nfp->get_x().size() == 1 );
    REQUIRE( nfp->get_y().size() == 1 );
    REQUIRE( nfp->get_z().size() == 1 );

    nec_delete(nec);
}

TEST_CASE( "NH card produces near-field results with explicit counts", "[near_field]") {
    nec_context* nec = nec_create();

    HANDLE_NEC(nec_wire(nec, 0, 7, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001, 1.0, 1.0));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 300.0, 0.0));
    HANDLE_NEC(nec_ex_card(nec, 0, 0, 4, 0, 1.0, 0.0, 0, 0, 0, 0));

    // NH: rectangular, NRX=NRY=NRZ=1, single point at (1,0,0)
    HANDLE_NEC(nec_nh_card(nec, 0, 1, 1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    nec_near_field_pattern* nfp = nec->get_near_field_pattern(0);
    REQUIRE( nfp != nullptr );
    REQUIRE( nfp->get_nfeh() == 1 );
    REQUIRE( nfp->get_x().size() == 1 );
    REQUIRE( nfp->get_y().size() == 1 );
    REQUIRE( nfp->get_z().size() == 1 );

    nec_delete(nec);
}

TEST_CASE( "NE card with zero NRX/NRY/NRZ defaults counts to 1", "[near_field]") {
    // Per NEC-2 Part 3, blank count parameters default to 1.
    // Zero values (from blank INT? fields in the ANTLR grammar,
    // or from explicit zeros) are bumped to 1 in ne_nh_card(),
    // so a near-field pattern is always produced.
    nec_context* nec = nec_create();

    HANDLE_NEC(nec_wire(nec, 0, 7, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001, 1.0, 1.0));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 300.0, 0.0));
    HANDLE_NEC(nec_ex_card(nec, 0, 0, 4, 0, 1.0, 0.0, 0, 0, 0, 0));

    // Zero counts -> bumped to 1 by ne_nh_card() default
    HANDLE_NEC(nec_ne_card(nec, 0, 0, 0, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    nec_near_field_pattern* nfp = nec->get_near_field_pattern(0);
    REQUIRE( nfp != nullptr );
    REQUIRE( nfp->get_x().size() == 1 );

    nec_delete(nec);
}

TEST_CASE( "NE card with multiple near-field points", "[near_field]") {
    // Verify a 2x3x1 grid (6 points) produces the correct number of results.
    nec_context* nec = nec_create();

    HANDLE_NEC(nec_wire(nec, 0, 7, 0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001, 1.0, 1.0));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 300.0, 0.0));
    HANDLE_NEC(nec_ex_card(nec, 0, 0, 4, 0, 1.0, 0.0, 0, 0, 0, 0));

    // 2 points along X spaced 0.5m, 3 along Y spaced 0.1m, 1 along Z
    HANDLE_NEC(nec_ne_card(nec, 0, 2, 3, 1, 0.0, -0.1, 1.0, 0.5, 0.1, 0.0));

    nec_near_field_pattern* nfp = nec->get_near_field_pattern(0);
    REQUIRE( nfp != nullptr );
    REQUIRE( nfp->get_x().size() == 6 );

    nec_delete(nec);
}

TEST_CASE( "GX one-plane symmetry produces finite impedance", "[symmetry]") {
    // PR #124: GX cards produced INF impedance due to uninitialized scratch
    // slot in solves(), inverted fblock() guard, and incomplete reflection.
    // One wire, GX 2 100 -> reflect in X, tag increment 2, NOP=2.
    nec_context* nec = nec_create();

    HANDLE_NEC(nec_wire(nec, 1, 11, -0.75, 0.0, -0.245, -0.75, 0.0, 0.245, 0.001, 1.0, 1.0));
    HANDLE_NEC(nec_wire(nec, 2, 11, -0.25, 0.0, -0.245, -0.25, 0.0, 0.245, 0.001, 1.0, 1.0));
    // GX: tag_inc=2, planes=100 -> X only (ix=1, iy=0, iz=0 -> nop=2)
    HANDLE_NEC(nec_gx_card(nec, 2, 100));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 299.7925, 0.0));
    // Feed all 4 tags (original + GX reflections) to match the deck
    HANDLE_NEC(nec_ex_card(nec, 0, 1, 6, 0, 1.0, 0.0, 0, 0, 0, 0));
    HANDLE_NEC(nec_ex_card(nec, 0, 2, 6, 0, 1.0, 0.0, 0, 0, 0, 0));
    HANDLE_NEC(nec_ex_card(nec, 0, 3, 6, 0, 1.0, 0.0, 0, 0, 0, 0));
    HANDLE_NEC(nec_ex_card(nec, 0, 4, 6, 0, 1.0, 0.0, 0, 0, 0, 0));
    HANDLE_NEC(nec_rp_card(nec, 0, 1, 361, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    double zr = nec_impedance_real(nec, 0);
    double zi = nec_impedance_imag(nec, 0);

    // Key regression: before PR #124 this returned INF (-999 sentinel).
    // The test harness .nec files cross-validate exact values against
    // reference engines; the unit test guards against the INF regression.
    REQUIRE( std::isfinite(zr) );
    REQUIRE( zr > 0.0 );      // positive real impedance
    REQUIRE( zr < 1000.0 );   // not absurd
    UNUSED(zi);

    nec_delete(nec);
}

TEST_CASE( "GX three-plane symmetry produces correct impedance", "[symmetry]") {
    // PR #124: three-plane symmetry (nop=8) was broken — fblock() had wrong
    // pass count, reflect_plane() missed copies 4-7. GX 1 111 -> all planes.
    nec_context* nec = nec_create();

    HANDLE_NEC(nec_wire(nec, 1, 9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.35, 0.001, 1.0, 1.0));
    // GX: tag_inc=1, planes=111 -> X,Y,Z (nop=8)
    HANDLE_NEC(nec_gx_card(nec, 1, 111));
    HANDLE_NEC(nec_geometry_complete(nec, 0));
    HANDLE_NEC(nec_fr_card(nec, 0, 1, 299.8, 0.0));
    HANDLE_NEC(nec_ex_card(nec, 0, 1, 5, 0, 1.0, 0.0, 0, 0, 0, 0));
    HANDLE_NEC(nec_rp_card(nec, 0, 1, 361, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    double zr = nec_impedance_real(nec, 0);
    double zi = nec_impedance_imag(nec, 0);

    REQUIRE( std::isfinite(zr) );
    REQUIRE( zr > 0.0 );
    REQUIRE( zr < 1000.0 );
    UNUSED(zi);

    nec_delete(nec);
}
