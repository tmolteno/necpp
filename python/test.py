import necpp as Necpp


def frequency_response():
  # Scan through frequencies from 1 to 30 MHz
  for f in range(1,30):
    nec = Necpp.nec_create()
    Necpp.nec_wire(nec, 1, 17, 0, 0, 2, 0, 0, 11, 0.1, 1, 1);
    Necpp.nec_geometry_complete(nec, 1, 0);
    Necpp.nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
    Necpp.nec_fr_card(nec, 0, 1, f, 0);
    Necpp.nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0, 0, 0, 0, 0);
    Necpp.nec_rp_card(nec, 0, 90, 1, 0,5,0,0, 0, 90, 1, 0, 0, 0);
    result_index = 0
    z = complex(Necpp.nec_impedance_real(nec,result_index), Necpp.nec_impedance_imag(nec,result_index))
    print "f=%0.2fMHz \t(%6.1f,%+6.1fI) Ohms" % (f, z.real, z.imag)
    Necpp.nec_delete(nec)

frequency_response()