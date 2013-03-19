require 'Vector3.rb'
require 'NECParser.rb'
require 'complex'
require 'NecWire.rb'
require 'necpp'
require 'AntennaStatistics'
require 'Parameters'

class WireArray < Array
	def to_s
		tag = 1
		ret = "\nCE\n"
		self.each{ |w|
			ret +=  "#{w.card(tag,1.0)}\n"
			tag += 1
		}
		ret += "GE\n"
		return ret
	end
end

class Antenna
	attr_accessor :wires, :freq, :function_count
	
	def initialize(comment)
		@comment = comment
		@wires = WireArray.new
		@length = 0
		@function_count=0
	end

	def add_comment(msg)
		@comment += msg.to_s
	end

	def cleanup
#		print "Cleanup #{self}\n"
		new_wires = WireArray.new
		@wires.each_index{ |i|
			a = @wires[i]
			if (a.length < 1e-3) then
				print "Deleting empty wire \n#{a}\n"
			else
				new_wires.push(a)
			end
		}
		@wires = new_wires

		deletions = Array.new
		@wires.each_index{ |i|
			a = @wires[i]
			@wires.each_index{ |j|
				if (i <  j) then
					b = @wires[j]					
					if a.similar(b) then
						print "Similar Wires \n#{a} \t#{b}\n"
						deletions.push(j)
					end
				end
			}
		}
		deletions.uniq!
		deletions.each{|k| @wires.delete_at(k) }

		@wires = intersections(@wires)
	end

	def intersections(old_wires)
#		print "Intersections\n #{old_wires}"
		new_wires = WireArray.new(old_wires)

		old_wires.each_index { |i|
			old_wires.each_index { |j|
				if (i > j) then
		      wires = old_wires[i].intersect(old_wires[j])
					if (wires != nil)
						if (wires.length == 2) then
							old_wires[i] = wires[0]
							old_wires[j] = wires[1]
  					else
							new_wires.delete_at(i)
							new_wires.delete_at(j)
#							print "Splitting [#{old_wires[i]}] & [#{old_wires[j]}]\n"
							wires.each{ |w| 
#								print "adding #{w}\n"
								new_wires.push(w)
		                    }
							return intersections(new_wires)
						end
					end
				end
			}
		}
		return new_wires
	end

	def add_wire(x0,x1,radius=0.0005)
		@wires.push(NecWire.new(x0,x1,radius))
		@length = @length + (x0 - x1).norm
  end

#	LD 0 <tag> 0 0 <R> <L> <C>
=begin
LDTAG (I2) Tag number; identifies the wire section(s) to be loaded by its (their) tag numbers. The next two parameters can be used to further specify certain segment(s) on the wire section(s). Blank or zero here implies that absolute segment numbers are being used in the next two parameters to identify segments. If the next two paramters ara blank or zero, all segments with tag LDTAG are loaded.

LDTAGF (13) - Equal to m specifies the mth segment of the set of segments whose tag numbers equal the tag number specified in the previous parameter. If the previous parameter (LDTAG) is zero, LDTAGF then specifies an absolute segment number. If both LDTAG and LDTAGF are zero, all segments will be loaded.

LDTAGT (14) - Equal to n specifies the nth segment of the set of segments whose tag numbers equal the tag number specified in the parameter LDTAG. This parameter must be greater than or equal to the previous param- eter. The loading specified is applied to each of the mth through nth segments of the set of segments having tags equal to LDTAG. Again if LDTAG is zero, these parameters refer to absolute segment numbers. If LDTAGT is left blank, it is set equal to the previous parameter (LDTAGF). 
=end
	def add_impedance(pos, r, l, c)
		@impedances.push(NecImpedance.new(pos,r,l,c))
	end

	def get_nec_geometry(wavelength)
		self.cleanup
		tag = 1
		ret = ""
		@comment.split('(').each{ |l| ret += "CM #{l}\n" }
		ret += "CE \n"
		@wires.each{ |w|
			ret +=  "#{w.card(tag,wavelength)}\n"
			tag += 1
		}
		return ret
	end
	
	def to_s
		ret = '' # ret = get_nec_file(300)
		@wires.each{ |w|
			ret +=  "add_wire(#{w})\n"
		}
		return ret
	end

	def rp_parameters
		theta_resolution = 5.0
		phi_resolution = 15.0
		theta_min = THETA_MIN
		theta_max = THETA_MAX
		phi_min = 0.0
		phi_max = 360.0
		n_theta = ((theta_max - theta_min) / theta_resolution).round + 1
		n_phi = ((phi_max - phi_min) / phi_resolution).round + 1
		d_theta = sf((theta_max - theta_min) / (n_theta-1.0), 0.001)
		d_phi = sf((phi_max - phi_min) / (n_phi-1.0), 0.001)
		# 		print "RP 0 #{n_theta} #{n_phi} 0000 0 0 #{d_theta} #{d_phi} 0 0\n"
		return theta_min, phi_min, n_theta, n_phi, d_theta, d_phi
	end
	
	# Simulate the antenna and return an AntennaStatistics object.
	def get_statistics(freq)
		self.cleanup
		begin
			nec = Necpp.nec_create
			wavelength = 3.0e8 / (freq * 1.0e6)
			tag = 1
			@wires.each{ |w|
				w.nec_wire(nec,tag,wavelength)
				tag += 1
			}
			ex_tag = 1
			ex_seg = 1

			ret = AntennaStatistics.new
			ret.wire_count = tag - 1
			return ret if tag == 1
		
			Necpp.nec_geometry_complete(nec, 1, 0);
			
			n_radial_wires = 17
			ground_screen_radius = 0.1
			Necpp.nec_gn_card(nec, 0, n_radial_wires, GROUND_DIELECTRIC_CONSTANT, GROUND_CONDUCTIVITY, ground_screen_radius, 0.001, 0.0, 0.0); 
			
			Necpp.nec_ld_card(nec,5, 0, 0, 0, WIRE_CONDUCTIVITY,0.0, 0.0)
			
			Necpp.nec_fr_card(nec, 0, 1, freq, 0.0)
			Necpp.nec_ex_card(nec, 0, ex_tag, ex_seg,  0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0)
			theta_min, phi_min, n_theta, n_phi, d_theta, d_phi = rp_parameters
			Necpp.nec_rp_card(nec, 0, n_theta, n_phi, 0,5,0,0, theta_min, phi_min, d_theta, d_phi, 1000.0, 0.0)
# # 			Necpp.nec_geometry_complete(nec, 0, 0);
# # 			Necpp.nec_gn_card(nec, 0, 4, 0.0, 0.0, 2.0, 0.005, 0.0, 0.0); # GN 2 0 0 0 2.0000 0.0050
# # 			Necpp.nec_fr_card(nec, 0, 1, @freq, 0.0)
# # 			Necpp.nec_ex_card(nec, 0, ex_tag, ex_seg,  0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0)
# # 			Necpp.nec_rp_card(nec, 4, 17, 45, 0,5,0,0, 0.0, 0.0, 5.0, 8.0, 0.0, 0.0)

			ret.impedance = Complex(Necpp.nec_impedance_real(nec,0),Necpp.nec_impedance_imag(nec,0))
			ret.gain_mean = Necpp.nec_gain_mean(nec,0)
			ret.gain_max = Necpp.nec_gain_max(nec,0)
			ret.gain_min = Necpp.nec_gain_min(nec,0)
			ret.gain_sd = Necpp.nec_gain_sd(nec,0)
			ret.gain_rhcp_mean = Necpp.nec_gain_rhcp_mean(nec,0)
			ret.gain_rhcp_max = Necpp.nec_gain_rhcp_max(nec,0)
			ret.gain_rhcp_min = Necpp.nec_gain_rhcp_min(nec,0)
			ret.gain_rhcp_sd = Necpp.nec_gain_rhcp_sd(nec,0)
			ret.gain_lhcp_mean = Necpp.nec_gain_lhcp_mean(nec,0)
			ret.gain_lhcp_max = Necpp.nec_gain_lhcp_max(nec,0)
			ret.gain_lhcp_min = Necpp.nec_gain_lhcp_min(nec,0)
			ret.gain_lhcp_sd = Necpp.nec_gain_lhcp_sd(nec,0)
		rescue CPPError => e
#			print "Exception: #{e}\n"
#			print "."
#			STDOUT.flush
			ret = nil
		ensure
			Necpp.nec_delete(nec)
		end

		return ret
	end

	def get_nec_file(freq)
		wavelength = 3.0e8 / (freq * 1.0e6)
		ex_tag = 1
		ex_seg = 1
		geo = get_nec_geometry(wavelength)
		geo += "GE 1\n"
		
		
		geo += "GN 0 0 0 0 #{GROUND_DIELECTRIC_CONSTANT} #{GROUND_CONDUCTIVITY} 0 0 0 0\n"
		geo += "LD 5 0 0 0 #{WIRE_CONDUCTIVITY} 0.0 0.0\n"
		geo += "FR 0 1 0 0 #{freq} 0\n"
		geo += "EX 0  #{ex_tag}  #{ex_seg}  0  1  0  0  0  0  0\n"
		theta_min, phi_min, n_theta, n_phi, d_theta, d_phi = rp_parameters
		geo += "RP 0 #{n_theta} #{n_phi} 0500 #{theta_min} #{phi_min} #{d_theta} #{d_phi} 0 0\n"
		geo += "EN\n"
		return geo
	end

	def old_stats(freq)
		self.cleanup
		np = NECParser.new(get_nec_file(freq))
		stats = np.get_stats
		stats.wire_count = @wires.length
		np.finalize
		return stats
	end
end
