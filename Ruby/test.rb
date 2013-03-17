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

def benchmark
 	nec = Necpp.nec_create
	Necpp.nec_wire(nec, 1, 9, 0.0, 0.0, 0.0, -0.0166, 0.0045, 0.0714, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 2, 7, -0.0166, 0.0045, 0.0714, -0.0318, -0.0166, 0.017, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 3, 7, -0.0318, -0.0166, 0.017, -0.0318, -0.0287, 0.0775, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 4, 11, -0.0318, -0.0287, 0.0775, -0.0318, 0.0439, 0.014, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 6, 5, -0.0318, 0.0045, 0.0624, -0.0106, 0.0378, 0.0866, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 7, 7, -0.0106, 0.0378, 0.0866, -0.0106, 0.0257, 0.023, 0.001, 1.0, 1.0)
	Necpp.nec_geometry_complete(nec, 1, 0);
	Necpp.nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
	Necpp.nec_fr_card(nec, 0, 1, 1600.0, 0.0)
	Necpp.nec_ex_card(nec, 0, 1, 1,  0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0)
#	Necpp.nec_rp_card(nec, 0, 90, 1, 0,0,0,0, 0, 90, 1, 0, 0, 0);
	Necpp.nec_rp_card(nec, 0, 17, 45, 0,0,0,0, 0, 0, 5, 8, 0, 0)
	100.times do antenna_stats(nec) end
	Necpp.nec_delete(nec)

  	nec = Necpp.nec_create
	Necpp.nec_wire(nec, 1, 17, 0, 0, 2, 0, 0, 11, 0.1, 1, 1);
	Necpp.nec_geometry_complete(nec, 1, 0);
	Necpp.nec_gn_card(nec, 1, 0, 0, 0, 0, 0, 0, 0);
	Necpp.nec_fr_card(nec, 0, 1, 30, 0);
	Necpp.nec_ex_card(nec, 0, 0, 5, 0, 1.0, 0, 0, 0, 0, 0);
	Necpp.nec_rp_card(nec, 0, 90, 1, 0,5,0,0, 0, 90, 1, 0, 0, 0);
	100.times do antenna_stats(nec) end
	Necpp.nec_delete(nec)

  nec = Necpp.nec_create
	Necpp.nec_wire(nec, 1, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 2, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 3, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 4, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 5, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 6, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 7, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 8, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 9, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0)
	Necpp.nec_wire(nec, 10, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0)
	Necpp.nec_geometry_complete(nec, 0, 0)
	
	Necpp.nec_gn_card(nec, -1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0)
	Necpp.nec_ld_card(nec, 5,0,0,0,3.72e7,0.0,0.0)
	Necpp.nec_pt_card(nec, -1, 0, 0, 0)
	Necpp.nec_ex_card(nec, 1, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	Necpp.nec_fr_card(nec, 0, 2, 2400.0, 100.0)
	Necpp.nec_rp_card(nec, 0, 1, 1, 0,5,0,0, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0)
	100.times do antenna_stats(nec) end
	Necpp.nec_delete(nec)
end


#print "Benchmark: #{Necpp.nec_benchmark}\n"

1.times do
	benchmark
end

