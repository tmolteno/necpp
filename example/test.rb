# Test Code for the nec2++ ruby module
# Author: Tim Molteno, tim@physics.otago.ac.nz
# GPL v3.
require 'necpp'
require 'complex'

=begin
  Collect some statistics about an antenna. If a frequency scan is done, then each individual
  result can be accessed using the index parameter. Here these are all zero.
=end
def antenna_stats(nec)
  result_index = 0
  z = Complex(Necpp.nec_impedance_real(nec,result_index), Necpp.nec_impedance_imag(nec,result_index))
  print "Impedance:#{z} Ohms\n"
  print "Gain:#{Necpp.nec_gain_max(nec,result_index)} dB Mean:#{Necpp.nec_gain_mean(nec,result_index)} dB, sd #{Necpp.nec_gain_sd(nec,result_index)} dB\n"
  print "RHCP:#{Necpp.nec_gain_rhcp_max(nec,result_index)} dB #{Necpp.nec_gain_rhcp_min(nec,result_index)} dB Mean:#{Necpp.nec_gain_rhcp_mean(nec,result_index)} dB, sd #{Necpp.nec_gain_rhcp_sd(nec,result_index)} dB\n"
  print "LHCP:#{Necpp.nec_gain_lhcp_max(nec,result_index)} dB Mean:#{Necpp.nec_gain_lhcp_mean(nec,result_index)} dB, sd #{Necpp.nec_gain_lhcp_sd(nec,result_index)} dB\n"
end

def frequency_response
  # Scan through frequencies from 1 to 30 MHz
  for f in (1..30) do
    nec = Necpp.nec_create
    Necpp.nec_wire(nec, 1, 17, 0, 0, 2, 0, 0, 11, 0.1, 1, 1);
    Necpp.nec_geometry_complete(nec, 1, 0);
    Necpp.nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
    Necpp.nec_fr_card(nec, 0, 1, f, 0);
    Necpp.nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0, 0, 0, 0, 0);
    Necpp.nec_rp_card(nec, 0, 90, 1, 0,5,0,0, 0, 90, 1, 0, 0, 0);
    result_index = 0
    z = Complex(Necpp.nec_impedance_real(nec,result_index), Necpp.nec_impedance_imag(nec,result_index))
    print "f=#{f}MHz \t#{z} Ohms\n"
    Necpp.nec_delete(nec)
  end
end


1.times do
  frequency_response
end

