#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "nec_context.h"
#include "nec_results.h"
#include "c_geometry.h"
#include "nec_wire.h"

#include <iostream>

#define REQUIRE_APPROX_EQUAL(a, b) \
  REQUIRE(std::abs((a) - (b)) < 3e-4)

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
    REQUIRE(geo->x[0] == Catch::Approx(0.2).margin(1e-6));
    REQUIRE(geo->y[0] == Catch::Approx(0.0).margin(1e-6));
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

TEST_CASE( "Optional intersection check bypass", "[intersection_check]") {
    // Regression test for #63: setting intersection check to false
    // should allow geometries that would otherwise trigger errors.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    // Disable intersection checking
    geo->set_intersection_check(false);

    // Two overlapping parallel wires that would normally be rejected
    // (same start point, same end point — definitely triggers intersection)
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    geo->wire(2, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    // Should succeed with intersection checking off
    REQUIRE_NOTHROW(nec.geometry_complete(0));
}

TEST_CASE( "Crossed wires sharing a node are accepted", "[segment_intersection]") {
    // Turnstile geometry: two dipoles crossing at the origin, each contributing
    // a segment end to the shared node. NEC-2 places no restriction on the angle
    // between connected wires, and every segment center clears the crossing
    // wire's volume by several radii.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 24, 0.0, 0.517, 0.0, 0.0, -0.517, 0.0, 0.006, 1.0, 1.0);
    geo->wire(2, 24, -0.517, 0.0, 0.0, 0.517, 0.0, 0.0, 0.006, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));
}

TEST_CASE( "Translated wire is tested at its final position", "[segment_intersection]") {
    // Tower geometry: a leg is built from the ground up, translated by its own
    // base height, and a second leg then fills the span it vacated. The two legs
    // meet end to end and share no volume in the final segment arrays.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(3, 20, 1.249, 0.0, 0.0, 1.249, 0.0, 11.24, 0.05, 1.0, 1.0);
    geo->move(0.0, 0.0, 0.0, 0.0, 0.0, 1.249, 0, 0, 0);
    geo->wire(1, 3, 1.249, 0.0, 0.0, 1.249, 0.0, 1.249, 0.05, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));
}

TEST_CASE( "Tapered wire is tested at its per-segment radius", "[segment_intersection]") {
    // Biconical geometry: a tapered cone meeting a short centre segment. The
    // cone's widest radius belongs to its far segment, and the segment actually
    // adjoining the centre is thin enough to clear it.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, 0.0, 0.0, -1.7831, 0.0, 0.0, -0.1024, 0.1514, 1.0, 0.6399);
    geo->wire(2, 1, 0.0, 0.0, -0.1024, 0.0, 0.0, 0.1024, 0.0254, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));
}

TEST_CASE( "Coincident wires are rejected", "[segment_intersection]") {
    // Two wires occupying the same space: every segment center of one lies on
    // the axis of a segment of the other.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    geo->wire(2, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    REQUIRE_THROWS(nec.geometry_complete(0));
}

TEST_CASE( "Wire penetrating another away from any node is rejected", "[segment_intersection]") {
    // Two single-segment wires crossing at their centres, sharing no end, so the
    // centre of each falls inside the volume of the other.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 1, -0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.01, 1.0, 1.0);
    geo->wire(2, 1, 0.0, -0.5, 0.0, 0.0, 0.5, 0.0, 0.01, 1.0, 1.0);
    REQUIRE_THROWS(nec.geometry_complete(0));
}

TEST_CASE( "Six wires meeting at one node are accepted", "[segment_junction]") {
    // Ground plane vertical on a mast: four radials, a radiator, and the mast
    // all contribute an end to the node at the origin. The connection record
    // links those six ends as a cycle, so the radials are not adjacent to the
    // mast in it, while every one of them is joined to it in the geometry.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 13, 0.0, 0.0, 0.0, -0.34, 0.0, -0.34, 0.0075, 1.0, 1.0);
    geo->wire(1, 13, 0.0, 0.0, 0.0, 0.34, 0.0, -0.34, 0.0075, 1.0, 1.0);
    geo->wire(1, 13, 0.0, 0.0, 0.0, 0.0, -0.34, -0.34, 0.0075, 1.0, 1.0);
    geo->wire(1, 13, 0.0, 0.0, 0.0, 0.0, 0.34, -0.34, 0.0075, 1.0, 1.0);
    geo->wire(2, 13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0075, 1.0, 1.0);
    geo->wire(3, 75, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.025, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));
    REQUIRE(geo->overlap_findings().empty());
}

TEST_CASE( "Legs meeting on the ground plane are accepted", "[segment_junction]") {
    // Two legs leaving one point on the ground plane 45 degrees apart, each
    // first center inside the other's volume and every later center clear of
    // it. build_connections() records a ground contact rather than a neighbour
    // at that end, so the junction is known from the coordinates alone.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 10, 0.0, 0.0, 0.0, 0.15307, 0.0, 0.36955, 0.02, 1.0, 1.0);
    geo->wire(2, 10, 0.0, 0.0, 0.0, -0.15307, 0.0, 0.36955, 0.02, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(1));
    REQUIRE(geo->overlap_findings().empty());
}

TEST_CASE( "Both signed ground modes reject below-plane segments", "[ground][zero_current]") {
    // NEC-2 Part 3: "When I1 is nonzero, no segment may extend below the
    // ground plane (X,Y plane) or lie in this plane." The zero-current mode
    // (GE -1) used to skip the checks that image mode (GE 1) applies.
    for (int gpflag : {1, -1}) {
      nec_context nec;
      nec.initialize();

      c_geometry* geo = nec.get_geometry();
      geo->wire(1, 5, 0.0, 0.0, -0.5, 0.0, 0.0, 0.5, 0.001, 1.0, 1.0);
      REQUIRE_THROWS(nec.geometry_complete(gpflag));
    }
}

TEST_CASE( "Both signed ground modes reject segments lying in the plane", "[ground][zero_current]") {
    // A segment whose axis lies in the ground plane is rejected in both
    // signed modes, while one that merely ends on the plane is accepted.
    for (int gpflag : {1, -1}) {
      nec_context nec;
      nec.initialize();

      c_geometry* geo = nec.get_geometry();
      geo->wire(1, 5, -0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.001, 1.0, 1.0);
      REQUIRE_THROWS(nec.geometry_complete(gpflag));
    }
}

TEST_CASE( "GE -1 accepts a segment ending on the ground plane", "[ground][zero_current]") {
    // The zero-current mode leaves a ground-touching end free (no image
    // connection), so the current goes to zero there; the segment is valid.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.001, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(-1));
    REQUIRE(geo->ground_connection() == -1);
}

TEST_CASE( "GE -1 accepts a horizontal wire close to the ground plane", "[ground][zero_current]") {
    // NEC-2 Part 3: for a horizontal wire less than 1e-3 x its segment length
    // above the ground, GE 1 connects every end to ground (or rejects the
    // wire); GE -1 must accept it with the ends left free.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, -1.0, 0.0, 1.0e-4, 1.0, 0.0, 1.0e-4, 0.001, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(-1));
}

TEST_CASE( "Ground connection without a ground model fails at simulate", "[ground][zero_current]") {
    // NEC-2 Part 3: "A positive or negative value of I1 does not cause a
    // ground to be included in the calculation... The ground parameters must
    // be specified on a program control card following the geometry cards."
    // Fail loudly instead of silently simulating free space.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.001, 1.0, 1.0);
    nec.geometry_complete(1);

    nec.ex_card(EXCITATION_VOLTAGE, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    REQUIRE_THROWS(nec.xq_card(0));
}

TEST_CASE( "Ground connection with a ground model simulates", "[ground][zero_current]") {
    // The same geometry with a GN card (perfect ground) is a valid simulation.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.001, 1.0, 1.0);
    nec.geometry_complete(1);
    nec.gn_card(1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    nec.ex_card(EXCITATION_VOLTAGE, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    REQUIRE_NOTHROW(nec.xq_card(0));
}

TEST_CASE( "Wire intruding past a shared node is rejected", "[segment_junction]") {
    // A feed leaving the hub of a radial at 20 degrees. Its first center clears
    // the radial segment it shares the hub with and comes to rest inside the
    // next one along, which it meets at no node.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->wire(2, 10, 0.0, 0.0, 0.0, 0.045, 0.0, 0.0, 0.003, 1.0, 1.0);
    geo->wire(3, 5, 0.0, 0.0, 0.0, 0.07518, 0.0, 0.02736, 0.001, 1.0, 1.0);
    REQUIRE_THROWS(nec.geometry_complete(0));
}

TEST_CASE( "Warning policy records the overlap and continues", "[segment_junction]") {
    // The geometry of the coincident-wire rejection, admitted under a policy
    // that reports rather than rejects. Every center of one wire lies on the
    // axis of the other, so both directions of every pair are recorded.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    geo->set_intersection_fatal(false);
    geo->wire(1, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    geo->wire(2, 5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.001, 1.0, 1.0);
    REQUIRE_NOTHROW(nec.geometry_complete(0));

    REQUIRE(false == geo->overlap_findings().empty());
    const nec_overlap_finding& first = geo->overlap_findings().front();
    REQUIRE(first.inside_tag == 1);
    REQUIRE(first.container_tag == 2);
    REQUIRE(first.distance <= first.container_radius);
}

TEST_CASE( "Helix rejects invalid segment_count", "[helix]") {
    // Regression test for #48: helix with segment_count < 1 must throw
    // rather than silently returning with no wires.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    REQUIRE_THROWS(geo->helix(1, 0, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.001));
    // Verify geometry has no segments from the failed helix
    REQUIRE(geo->n_segments == 0);
}

TEST_CASE( "axis_distance throws on a zero-length wire", "[wire]") {
    // A wire with coincident endpoints has a squared-length denominator of
    // zero in axis_distance(); it must throw rather than divide by zero
    // (which would produce NaN and silently poison the overlap check).
    nec_3vector p(1.0, 2.0, 3.0);
    nec_wire degenerate(p, p, 0.001, 1);

    nec_3vector q(0.0, 0.0, 0.0);
    REQUIRE_THROWS(degenerate.axis_distance(q));

    // A non-degenerate wire still yields a finite distance.
    nec_3vector a(0.0, 0.0, 0.0);
    nec_3vector b(1.0, 0.0, 0.0);
    nec_wire unit(a, b, 0.001, 1);
    nec_float d = unit.axis_distance(nec_3vector(0.5, 1.0, 0.0));
    REQUIRE( std::isfinite(d) );
    REQUIRE( d == Catch::Approx(1.0).margin(1e-12) );
}

TEST_CASE( "EX type-4 source on a patch centroid zeroes that patch", "[etmns][current_source]") {
    // Regression test for #128: an elementary current source (EX type 4)
    // placed exactly on a patch centroid used to skip assigning that patch's
    // right-hand-side rows, leaving stale/garbage values in the persisting
    // current vector (the coincident patch printed 9.88E-324 and corrupted
    // the currents coupled to it). An element centered on the source point
    // receives no field from it, so its slots must read exactly zero.
    nec_context nec;
    nec.initialize();

    c_geometry* geo = nec.get_geometry();
    // A wire (3 segments) keeps n > 0 so the patch tangential slots are based
    // at a nonzero index, exercising the offset layout.
    geo->wire(1, 3, 0.5, -0.05, 0.0, 0.5, 0.05, 0.0, 0.001, 1.0, 1.0);
    // Three patches tile a plate; the middle sits exactly at the source point
    // (the origin). SP ns=0: arbitrary shape at (x,y,z), elevation 90, az 0.
    nec.sp_card(0,  0.1, 0.0, 0.0, 90.0, 0.0, 0.01);
    nec.sp_card(0,  0.0, 0.0, 0.0, 90.0, 0.0, 0.01);
    nec.sp_card(0, -0.1, 0.0, 0.0, 90.0, 0.0, 0.01);
    nec.geometry_complete(0);

    nec.fr_card(0, 1, 300.0, 0.0);
    // z-directed elementary current source at the origin (ETAI = 90).
    nec.ex_card(EXCITATION_CURRENT, 0, 0, 0, 0.0, 0.0, 0.0, 90.0, 0.0, 1.0);
    nec.rp_card(0, 19, 37, 1000, 0, 0, 0, 0.0, 0.0, 10.0, 10.0, 0.0, 0.0);

    nec_structure_currents* sc = nec.get_structure_currents(0);
    REQUIRE(sc != nullptr);

    // Patches: the coincident (middle, index 1) one must read zero in every
    // component. Before the fix its slots carried stale values.
    REQUIRE(sc->get_patch_e_x().size() == 3);
    REQUIRE(std::abs(sc->get_patch_e_x()[1]) < 1e-9);
    REQUIRE(std::abs(sc->get_patch_e_y()[1]) < 1e-9);
    REQUIRE(std::abs(sc->get_patch_e_z()[1]) < 1e-9);

    // A z-directed source imposes antisymmetric wire currents about the
    // source plane: segments 1 and 3 are equal and opposite, the center nulls.
    REQUIRE(sc->get_current().size() == 3);
    const nec_complex c1 = sc->get_current()[0];
    const nec_complex c2 = sc->get_current()[1];
    const nec_complex c3 = sc->get_current()[2];
    REQUIRE(std::abs(c2) < 1e-8);
    REQUIRE(std::real(c1) == Catch::Approx(-std::real(c3)).margin(1e-12));
    REQUIRE(std::imag(c1) == Catch::Approx(-std::imag(c3)).margin(1e-12));
    // Absolute values agree with the reference (cross-validated against
    // nec2c / xnec2c / NEC-2 Fortran).
    REQUIRE(std::real(c1) == Catch::Approx(6.9547e-05).margin(1e-7));
    REQUIRE(std::imag(c1) == Catch::Approx(-1.5977e-05).margin(1e-7));
}
