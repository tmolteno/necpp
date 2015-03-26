import necpp

#
# \brief Using the python interface to the C-style nec2++ API
# 
# This is built in the 'python' directory of the source code distribution
#

def handle_nec(result):
  if (result != 0):
    print nec_error_message()

def frequency_response():
  # Scan through frequencies from 1 to 30 MHz
  for f in range(1,30):
    nec = necpp.nec_create()
    handle_nec(necpp.nec_wire(nec, 1, 17, 0, 0, 2, 0, 0, 11, 0.1, 1, 1))
    handle_nec(necpp.nec_geometry_complete(nec, 1, 0))
    handle_nec(necpp.nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0))
    handle_nec(necpp.nec_fr_card(nec, 0, 1, f, 0))
    handle_nec(necpp.nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0, 0, 0, 0, 0))
    handle_nec(necpp.nec_rp_card(nec, 0, 90, 1, 0,5,0,0, 0, 90, 1, 0, 0, 0))
    result_index = 0
    z = complex(necpp.nec_impedance_real(nec,result_index), necpp.nec_impedance_imag(nec,result_index))
    print "f=%0.2fMHz \t(%6.1f,%+6.1fI) Ohms" % (f, z.real, z.imag)
    necpp.nec_delete(nec)

frequency_response()

