from necpp import *


import unittest

class TestDipoleGain(unittest.TestCase):

  def handle_nec(self, result):
    if (result != 0):
      print nec_error_message()
    self.assertEqual(result,0)

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
    self.handle_nec(nec_sp_card(nec, 0, 0.0, 0.0, 0.1, 90.0, 0.0, 0.04))
    
    self.handle_nec(nec_wire(nec, 1, 4, 0., 0.0, 0.1, 0.0,  0.0, 0.3, .001, 1.0, 1.0))
    self.handle_nec(nec_wire(nec, 2, 2, 0., 0.0, 0.3, 0.15, 0.0, 0.3, .001, 1.0, 1.0))
    self.handle_nec(nec_wire(nec, 3, 2, 0., 0.0, 0.3, -.15, 0.0, 0.3, .001, 1.0, 1.0))

    self.handle_nec(nec_geometry_complete(nec, 1))
    self.handle_nec(nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0))
    
    self.handle_nec(nec_ex_card(nec, 0, 1, 1, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    self.handle_nec(nec_rp_card(nec, 0,10,4,1,0,0,1,0.0,0.0,10.0,30.0, 0, 0))
    
    self.assertAlmostEqual(nec_gain_max(nec,0),5.076,3)
    
    gmax = -999.0
    
    for theta_index in range(0,10):
      for phi_index in range(0,4):
        g = nec_gain(nec,0,theta_index, phi_index)
        gmax = max(g, gmax)
        
    self.assertAlmostEqual(gmax, nec_gain_max(nec,0), 5 )

    nec_delete(nec)
    

if __name__ == '__main__':
  unittest.main()