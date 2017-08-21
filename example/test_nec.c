#include <libnecpp.h>
#include <stdio.h>

/*
  C Example using nec2++ as a library
  
  Enter the following file into test_nec.c, and compile with
  gcc -o test_nec test_nec.c -L . -lnecpp -lm
*/

#define NEC_ERROR_HANDLE(_x) { if (_x != 0) printf("Error: %s\n", nec_error_message()); }

void seven_wire_antenna() {
  nec_context* nec;	
  nec = nec_create();
  
  NEC_ERROR_HANDLE(nec_wire(nec, 1, 9, 0.0, 0.0, 0.0, -0.0166, 0.0045, 0.0714, 0.001, 1.0, 1.0));
  NEC_ERROR_HANDLE(nec_wire(nec, 2, 7, -0.0166, 0.0045, 0.0714, -0.0318, -0.0166, 0.017, 0.001, 1.0, 1.0));
  NEC_ERROR_HANDLE(nec_wire(nec, 3, 7, -0.0318, -0.0166, 0.017, -0.0318, -0.0287, 0.0775, 0.001, 1.0, 1.0));
  NEC_ERROR_HANDLE(nec_wire(nec, 4, 11, -0.0318, -0.0287, 0.0775, -0.0318, 0.0439, 0.014, 0.001, 1.0, 1.0));
/*  NEC_ERROR_HANDLE(nec_wire(nec, 5, 7, -0.0318, 0.0439, 0.014, -0.0318, 0.0045, 0.0624, 0.001, 1.0, 1.0));*/
  NEC_ERROR_HANDLE(nec_wire(nec, 6, 5, -0.0318, 0.0045, 0.0624, -0.0106, 0.0378, 0.0866, 0.001, 1.0, 1.0));
  NEC_ERROR_HANDLE(nec_wire(nec, 7, 7, -0.0106, 0.0378, 0.0866, -0.0106, 0.0257, 0.023, 0.001, 1.0, 1.0));
  NEC_ERROR_HANDLE(nec_geometry_complete(nec, 1, 0));
  NEC_ERROR_HANDLE(nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0));
  NEC_ERROR_HANDLE(nec_fr_card(nec, 0, 1, 1600.0, 0.0));
  NEC_ERROR_HANDLE(nec_ex_card(nec, 0, 1, 1,  0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0));
  NEC_ERROR_HANDLE(nec_rp_card(nec, 0, 17, 45, 0,5,0,0, 0, 0, 5, 8, 0, 0));
  
  printf("Impedance: %f, %f\n",nec_impedance_real(nec,0), nec_impedance_imag(nec,0));
  printf("Gain: %f, %f +/- %f dB\n",nec_gain_max(nec,0), nec_gain_mean(nec,0), nec_gain_sd(nec,0));
  printf("RHCP Gain: %f, %f +/- %f dB\n",nec_gain_rhcp_max(nec,0), nec_gain_rhcp_mean(nec,0), nec_gain_rhcp_sd(nec,0));
  printf("LHCP Gain: %f, %f +/- %f dB\n",nec_gain_lhcp_max(nec,0), nec_gain_lhcp_mean(nec,0), nec_gain_lhcp_sd(nec,0));
  
  nec_delete(nec);
}

void simple_example() {
  /*  GW 0 9 0. 0. 2. 0. 0. 7 .1
      GE 1
      FR 0 1 0 30.
      EX 0 5 0 1.
      GN 1
      RP 0 90 1 0000 0 90 1 0 
  */
  
  nec_context* nec;	
  nec = nec_create();
  
  nec_wire(nec, 0, 9, 0, 0, 2, 0, 0, 7, 0.1, 1, 1);
  nec_geometry_complete(nec, 1, 0);
  nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
  nec_fr_card(nec, 0, 1, 30, 0);
  nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0, 0, 0, 0, 0);
  nec_rp_card(nec, 0, 90, 1, 0,5,0,0, 0, 90, 1, 0, 0, 0);
  
  printf("Gain: %f, %f +/- %f dB\n",nec_gain_max(nec,0), nec_gain_mean(nec,0), nec_gain_sd(nec,0));
  printf("RHCP Gain: %f, %f +/- %f dB\n",nec_gain_rhcp_max(nec,0), nec_gain_rhcp_mean(nec,0), nec_gain_rhcp_sd(nec,0));
  printf("LHCP Gain: %f, %f +/- %f dB\n",nec_gain_lhcp_max(nec,0), nec_gain_lhcp_mean(nec,0), nec_gain_lhcp_sd(nec,0));
  
  nec_delete(nec);
}

void example3() {
  /*
    CMEXAMPLE 3. VERTICAL HALF WAVELENGTH ANTENNA OVER GROUND 
    CM           EXTENDED THIN WIRE KERNEL USED 
    CM           1. PERFECT GROUND 
    CM           2. IMPERFECT GROUND INCLUDING GROUND WAVE AND RECEIVING 
    CE              PATTERN CALCULATIONS 
    GW 0 9 0. 0. 2. 0. 0. 7. .03 
    GE 1 
    EK 
    FR 0 1 0 0 30. 
    EX 0 0 5 0 1. 
    GN 1
    RP 0 10 2 1301 0. 0. 10. 90. 
    GN 0 0 0 0 6. 1.000E-03  
    RP 0 10 2 1301 0. 0. 10. 90. 
    RP 1 10 1 0 1. 0. 2. 0. 1.000E+05 
    EX 1 10 1 0 0. 0. 0. 10.
    PT 2 0 5 5 
    XQ 
    EN
*/
  nec_context* nec;     
  nec = nec_create();
  
  nec_wire(nec, 0, 9, 0., 0.0, 2.0, 0.0, 0.0, 7.0, 0.03, 1.0, 1.0);
  nec_geometry_complete(nec, 1, 0);
  nec_ek_card(nec, 0);
  nec_fr_card(nec, 0, 1, 30., 0 );
  nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
  nec_rp_card(nec, 0,10,2,1,3,0,1,0.0,0.0,10.0,90.0, 0, 0);

  printf("Impedance %f\n", nec_impedance_real(nec,0)); // 83.7552291016712)
   
  nec_gn_card(nec, 0, 0, 6.0, 1.000E-03, 0.0, 0.0, 0.0, 0.0);
  nec_rp_card(nec, 0,10,2,1,3,0,1, 0.0,0.0,10.0,90.0, 0, 0);

  printf("Impedance 2 %f\n", nec_impedance_real(nec,0)); // 83.7552291016712)
}
int main() {
  simple_example();
  example3();
  seven_wire_antenna();
  return 0;
}
