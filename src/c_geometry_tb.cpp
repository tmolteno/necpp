#include "catch.hpp"

#include "libnecpp.h"
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
