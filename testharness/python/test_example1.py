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

if __name__ == '__main__':
  unittest.main()