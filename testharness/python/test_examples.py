from necpp import *


import unittest

class TestDipoleGain(unittest.TestCase):

  def handle_nec(self, result):
    if (result != 0):
      print nec_error_message()
    self.assertEqual(result,0)

  def test_example1(self):
    nec = nec_create()

    '''
    CE EXAMPLE 1. CENTER FED LINEAR ANTENNA
    GW 0 7 0. 0. -.25 0. 0. .25 .001
    GE
    EX 0 0 4 0 1.
    XQ
    LD 0 0 4 4 10. 3.000E-09 5.300E-11
    PQ
    NE 0 1 1 15 .001 0 0 0. 0. .01786
    EN    
    '''
    self.handle_nec(nec_wire(nec, 0, 7, 0., 0., .75, 0., 0., 1.25, .001, 1.0, 1.0))
    self.handle_nec(nec_geometry_complete(nec, 1, 0))
    self.handle_nec(nec_ex_card(nec, 0, 0, 4,0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    self.handle_nec(nec_xq_card(nec, 0))
    self.handle_nec(nec_ld_card(nec, 0, 0, 4, 4, 10., 3.000E-09, 5.300E-11))
    self.handle_nec(nec_pq_card(nec, 0, 0, 0, 0))
    self.handle_nec(nec_ne_card(nec, 0, 1, 1, 15, .001, 0, 0, 0., 0., .01786))
    
    self.assertAlmostEqual(nec_impedance_real(nec,0),82.69792906662622)
    self.assertAlmostEqual(nec_impedance_imag(nec,0),46.30603888063429)
    
    nec_delete(nec)

  def test_example3(self):
    '''
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
    '''
    nec = nec_create()
    self.handle_nec(nec_wire(nec, 0, 9, 0., 0.0, 2.0, 0.0, 0.0, 7.0, 0.03, 1.0, 1.0))
    self.handle_nec(nec_geometry_complete(nec, 1, 0))
    self.handle_nec(nec_ek_card(nec, 0))
    self.handle_nec(nec_fr_card(nec, 0, 1, 30., 0 ))
    self.handle_nec(nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    self.handle_nec(nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0))
    self.handle_nec(nec_rp_card(nec, 0,10,2,1,3,0,1,0.0,0.0,10.0,90.0, 0, 0))

    self.assertAlmostEqual(nec_impedance_real(nec,0),83.7552291016712)
    self.assertAlmostEqual(nec_impedance_imag(nec,0),45.32205265591289)
    self.assertAlmostEqual(nec_gain_max(nec,0),8.393875976328134)
   
    self.handle_nec(nec_gn_card(nec, 0, 0, 6.0, 1.000E-03, 0, 0, 0, 0))
    self.handle_nec(nec_rp_card(nec, 0,10,2,1,3,0,1, 0.0,0.0,10.0,90.0, 0, 0))
    
    self.assertAlmostEqual(nec_impedance_real(nec,1),86.415,3)
    self.assertAlmostEqual(nec_impedance_imag(nec,1),47.822,3)
    self.assertAlmostEqual(nec_gain_max(nec,1),1.44837,3)

    self.handle_nec(nec_rp_card(nec, 1,10,1,0,0,0,0, 1.0,0.0,2.0,0.0, 1.000E+05, 0))
    # Not sure what to check here.
    
    self.handle_nec(nec_ex_card(nec, 1, 10, 1, 0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0))
    self.handle_nec(nec_pt_card(nec, 2, 0, 5, 5))
    # Not sure what to check here.

    nec_delete(nec)


  def test_example4(self):
    '''
    CEEXAMPLE 4. T ANTENNA ON A BOX OVER PERFECT GROUND
    SP 0 0 .1 .05 .05 0. 0. .01
    SP 0 0 .05 .1 .05 0. 90. .01
    GX 0 110
    SP 0 0 0. 0. .1 90. 0. .04
    GW 1 4 0. 0. .1 0. 0. .3 .001
    GW 2 2 0. 0. .3 .15 0. .3 .001
    GW 3 2 0. 0. .3 -.15 0. .3 .001
    GE 1
    GN 1
    EX 0 1 1 0 1.
    RP 0 10 4 1001 0. 0. 10. 30.
    EN
    '''
    nec = nec_create()
    self.handle_nec(nec_sp_card(nec, 0, 0.1, 0.05, 0.05, 0.0, 0.0, 0.01))
    self.handle_nec(nec_sp_card(nec, 0, .05, .1, .05, 0.0, 90.0, 0.01))
    self.handle_nec(nec_gx_card(nec, 0, 110))
    self.handle_nec(nec_sp_card(nec, 0, 0, 0.0, 0.0, 0.1, 90.0, 0.0, 0.04))
    
    self.handle_nec(nec_wire(nec, 1, 4, 0., 0.0, 0.1, 0.0,  0.0, 0.3, .001, 1.0, 1.0))
    self.handle_nec(nec_wire(nec, 2, 2, 0., 0.0, 0.3, 0.15, 0.0, 0.3, .001, 1.0, 1.0))
    self.handle_nec(nec_wire(nec, 3, 2, 0., 0.0, 0.3, -.15, 0.0, 0.3, .001, 1.0, 1.0))

    self.handle_nec(nec_geometry_complete(nec, 1, 0))
    self.handle_nec(nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0))
    
    self.handle_nec(nec_ex_card(nec, 0, 1, 1, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    self.handle_nec(nec_rp_card(nec, 0,10,4,1,0,0,1,0.0,0.0,10.0,30.0, 0, 0))
    nec_delete(nec)
    

if __name__ == '__main__':
  unittest.main()