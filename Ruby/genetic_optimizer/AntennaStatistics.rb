require 'complex'

# A little class that does the fitness calculations for the Simulation Clients
class AntennaStatistics
	attr_accessor :wire_count, :impedance, :gain_mean, :gain_max, :gain_min, :gain_sd, :gain_rhcp_mean, :gain_rhcp_max, :gain_rhcp_min, :gain_rhcp_sd, :gain_lhcp_mean, :gain_lhcp_max, :gain_lhcp_min, :gain_lhcp_sd

	GAIN_MIN = -999.0
	GAIN_SD_DEFAULT = 9.0e9

	def initialize
		@wire_count = 0
		@impedance = nil
		@gain_mean = nil
		@gain_max = nil
		@gain_min = nil
		@gain_sd = nil
		@gain_rhcp_mean = nil
		@gain_rhcp_max = nil
		@gain_rhcp_min = nil
		@gain_rhcp_sd = nil
		@gain_lhcp_mean = nil
		@gain_lhcp_max = nil
		@gain_lhcp_min = nil
		@gain_lhcp_sd = nil
	end

	# Return the return loss assuming the Antenna
	# is fed with impedance z0_in
	def reflection_coeff(z0_in)
		return nil if @impedance == nil
		z0 = Complex(z0_in)
		return (@impedance - z0) / (@impedance + z0)
	end

	# Return the Voltage Standing Wave Ratio (VSWR) assuming the Antenna
	# is fed with impedance z0_in
	def vswr(z0_in)
		return nil if @impedance == nil
		rl = reflection_coeff(z0_in)
		return (1 + rl.abs) / (1 - rl.abs)
	end

	def to_s
		begin
			ret = "Stats: rc=#{sf(reflection_coeff(50.0).abs,0.01)}, vswr=#{sf(vswr(50.0),0.01)} "
			ret += "Gain: [#{sf(@gain_min,0.01)}, #{sf(@gain_max,0.01)}] #{sf(@gain_mean,0.01)} +- #{sf(@gain_sd,0.01)} "
			ret += "RHCP: [#{sf(@gain_rhcp_min,0.01)}, #{sf(@gain_rhcp_max,0.01)}] #{sf(@gain_rhcp_mean,0.01)} +- #{sf(@gain_rhcp_sd,0.01)} "
			ret += "LHCP: [#{sf(@gain_lhcp_min,0.01)}, #{sf(@gain_lhcp_max,0.01)}] #{sf(@gain_lhcp_mean,0.01)} +- #{sf(@gain_lhcp_sd,0.01)}"
		rescue Exception => e
			return ret
		end
	end
end
