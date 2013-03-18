require 'AntennaStatistics'
require 'tempfile'

def deg2rad(theta)
	return theta*(Math::PI / 180.0)
end

def to_db(x)
	return AntennaStatistics::GAIN_MIN if x == 0
	return 10.0*Math.log10(x)
end

def from_db(db)
	return 0.0 if db < -100
	return 10.0**(db / 10.0)
end


class NECParser
	def initialize(geo)
		fbase = "test#{rand(10000)}"

		@outFile = Tempfile.new("#{fbase}.out")
		@necFile = Tempfile.new("#{fbase}.nec")
		@necFile.print(geo)
		@necFile.close
		sys = IO.popen("nec2++ -i #{@necFile.path} -o #{@outFile.path} 2>&1")
		b = sys.readlines
		b.join
		sys.close

		@lines = Array.new
		@outFile.each_line do |line|
			@lines.push(line.clone)
		end
		@necFile.close
		@impedance = nil
		@rp = nil
		@gains = nil
		@thetas = Array.new
		self.parse
	end


	def finalize
		@necFile.unlink
		@outFile.unlink
	end

	def impedance
		return @impedance
	end

	def scales
		return @scales unless @scales == nil
		@scales = Array.new
		n_theta = @thetas.uniq.length.to_f
		@thetas.each{ |t| @scales.push((Math.sin(t) * (n_theta - 1) / n_theta) + 1.0/n_theta) }
		return @scales
	end

	def max(pattern)
		return nil if pattern == nil
 		return pattern.max
	end

	def mean(pattern)
		return nil if pattern == nil
		f = scales
		sum = 0
 		pattern.each_index { |i| sum += pattern[i] * f[i] }
 		return sum / (pattern.length * 2.0 / Math::PI)
	end

	def sd(pattern, mean)
		return nil if pattern == nil
		f = scales
		sum = 0
		pattern.each_index { |i| sum += (pattern[i] - mean)**2 * f[i] }
		return Math.sqrt(sum / (pattern.length * 2.0 / Math::PI))
	end

	def gain_sd
		return @gain_sd unless @gain_sd == nil
		@gain_sd = sd(@gains, self.gain_mean)
	end

	def gain_mean
		return @gain_mean unless @gain_mean == nil
		@gain_mean = mean(@gains)
	end

	def gain_max
		return @gain_max unless @gain_max == nil
		@gain_max = max(@gains)
	end

######################### RHCP #############################
	def gain_rhcp_max
		return @gains_rhcp_max unless @gains_rhcp_max == nil
		@gains_rhcp_max = max(@gains_rhcp)
	end

	def gain_rhcp_mean
		return @gain_rhcp_mean unless @gain_rhcp_mean == nil
		@gain_rhcp_mean = mean(@gains_rhcp)
	end

	def gain_rhcp_sd
		return @gain_rhcp_sd unless @gain_rhcp_sd == nil
		@gain_rhcp_sd = sd(@gains_rhcp, self.gain_rhcp_mean)
	end

######################### LHCP #############################
	def gain_lhcp_max
		return @gains_lhcp_max unless @gains_lhcp_max == nil
		@gains_lhcp_max = max(@gains_lhcp)
	end

	def gain_lhcp_mean
		return @gain_lhcp_mean unless @gain_lhcp_mean == nil
		@gain_lhcp_mean = mean(@gains_lhcp)
	end

	def gain_lhcp_sd
		return @gain_lhcp_sd unless @gain_lhcp_sd == nil
		@gain_lhcp_sd = sd(@gains_lhcp, self.gain_lhcp_mean)
	end

	def get_stats
		stats = AntennaStatistics.new
		stats.impedance = impedance

		stats.gain_max = gain_max
		stats.gain_mean = gain_mean
		stats.gain_sd = gain_sd

		stats.gain_rhcp_max = gain_rhcp_max
		stats.gain_rhcp_mean = gain_rhcp_mean
		stats.gain_rhcp_sd = gain_rhcp_sd

		stats.gain_lhcp_max = gain_lhcp_max
		stats.gain_lhcp_mean = gain_lhcp_mean
		stats.gain_lhcp_sd = gain_lhcp_sd
		return stats
	end

	def parse
		nl = -1
		@lines.each{ |l|
			if (@impedance == nil) then
				if (nl == 0) then
					# 1     1  1.0000E+00  0.0000E+00  1.2616E-03 -2.3212E-03  1.8075E+02  3.3257E+02  1.2616E-03 -2.3212E-03  6.3078E-04 
					ret_r = l.scan(/[0-9.eE\+-]+/)[6].to_f
					ret_i = l.scan(/[0-9.eE\+-]+/)[7].to_f
					@impedance = Complex(ret_r, ret_i)
					nl = -1
				end
				if (l.index("PARAMETERS") != nil) then
					nl = 3
				end
				if (nl >= 0) then
					nl = nl - 1
				end
			elsif (@rp == nil) then
				#                              ---------- RADIATION PATTERNS -----------
				# 
				#  ---- ANGLES -----     --- DIRECTIVE GAINS ---      ---- POLARIZATION ----   ---- E(THETA) ----    ----- E(PHI) ------
				#   THETA      PHI       VERTC   HORIZ   TOTAL       AXIAL      TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE
				#  DEGREES   DEGREES        DB       DB       DB       RATIO   DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES
				if (nl == 0) then
					@rp = 1
					@gains = Array.new
					@axial_ratios = Array.new
					@gains_rhcp = Array.new
					@gains_lhcp = Array.new
				end
				if (l.index("RADIATION PATTERN") != nil) then
					nl = 4
				end
				if (nl >= 0) then
					nl = nl - 1
				end
			elsif (@rp == 1) then
				#                              ---------- RADIATION PATTERNS -----------
				# 
				#  ---- ANGLES -----     --- DIRECTIVE GAINS ---      ---- POLARIZATION ----   ---- E(THETA) ----    ----- E(PHI) ------
				#   THETA      PHI       VERTC   HORIZ   TOTAL       AXIAL      TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE
				#  DEGREES   DEGREES        DB       DB       DB       RATIO   DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES
				floats = l.scan(/[0-9.eE\+-]+/)
				theta = floats[0].to_f
				phi = floats[1].to_f
				db = floats[4].to_f
				a = floats[5].to_f
				tilt = floats[6].to_f
				sense = l.scan(/[A-Z]+/)[0]
				if (floats.length == 0) then
					break
				end
				if (sense == "RIGHT")
					a = -a
				end
				gain = from_db(db)
				lhcp_f =(1+2*a+a**2)/(2*(1+a**2))
				rhcp_f =(1-2*a+a**2)/(2*(1+a**2))
				rh_gain = gain * rhcp_f
				lh_gain = gain * lhcp_f
#				print "RP(#{theta},#{phi}) #{db} dB #{sense} #{a} rh=#{rhcp_f} lh=#{lhcp_f} db=(#{to_db(rh_gain)}, #{to_db(lh_gain)})\n"
				@gains_rhcp.push(to_db(rh_gain))
				@gains_lhcp.push(to_db(lh_gain))

				@thetas.push(deg2rad(theta))
				@gains.push(to_db(gain))
				@axial_ratios.push(a)
			end
		}
	end

end

if __FILE__ == $0
require "Antenna.rb"
	ant = Antenna.new("Linden")

=begin
CM Linden's Seven wire genetic antenna
CE
-GW 1 7 0.0 0.0 0.0 -0.0166 0.0045 0.0714 0.001
-GW 1 7 -0.0166 0.0045 0.0714 -0.0318 -0.0166 0.0170 0.001
-GW 1 7 -0.0318 -0.0166 0.0170 -0.0318 -0.0287 0.0775 0.001
-GW 1 7 -0.0318 -0.0287 0.0775 -0.0318 0.0439 0.0140 0.001
GW 1 7 -0.0318 0.0439 0.0140 -0.0318 0.0045 0.0624 0.001
GW 1 7 -0.0318 0.0045 0.0624 -0.0106 0.0378 0.0866 0.001
GW 1 7 -0.0106 0.0378 0.0866 -0.0106 0.0257 0.0230 0.001
GE 1 
GN 1 0 0 0 0 0 0 0 0 0 
FR 0 1 0 0 1600 0
EX 0  1  1  0  1  0  0  0  0  0  
RP 0 40 90 0000 0 0 2 4 0 0 
=end
	ant.add_wire(Vector3.new(0.0, 0.0, 0.0), Vector3.new(-0.0166, 0.0045, 0.0714))
	ant.add_wire(Vector3.new(-0.0166, 0.0045, 0.0714), Vector3.new(-0.0318, -0.0166, 0.0170))
	ant.add_wire(Vector3.new(-0.0318, -0.0166, 0.0170), Vector3.new(-0.0318, -0.0287, 0.0775))
	ant.add_wire(Vector3.new(-0.0318, -0.0287, 0.0775), Vector3.new(-0.0318, 0.0439, 0.0140))
	ant.add_wire(Vector3.new(-0.0318, 0.0439, 0.0140), Vector3.new(-0.0318, 0.0045, 0.0624))
	ant.add_wire(Vector3.new(-0.0318, 0.0045, 0.0624), Vector3.new(-0.0106, 0.0378, 0.0866))
	ant.add_wire(Vector3.new(-0.0106, 0.0378, 0.0866), Vector3.new(-0.0106, 0.0257, 0.0230))

	ant.cleanup
	print "Antenna #{ant}\n"

end