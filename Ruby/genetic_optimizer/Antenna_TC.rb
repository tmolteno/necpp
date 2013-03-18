require 'Antenna'
require 'Individual'

require 'test/unit'
class Antenna_TC < Test::Unit::TestCase

	def testa
		a0 = Vector3.new(0,0,0)
		a1 = Vector3.new(1,1,0)
	
		b0 = Vector3.new(1,0,0)
		b1 = Vector3.new(0,1,0)
	
		c0 = Vector3.new(0.3,0,0)
		c1 = Vector3.new(0.3,1,0)
	
		ant = Antenna.new("Test A")
		ant.add_wire(a0,a1)
		ant.add_wire(b0,b1)
		ant.add_wire(c0,c1)
		return ant
	end

	def testb
		ant = Antenna.new("Test B")
		ant.add_wire(Vector3.new(-1.085,-0.135,-2.725), Vector3.new(-2.38,0.155,-2.86))
		ant.add_wire(Vector3.new(-2.38,0.155,-2.86), Vector3.new(-1.275,0.095,-2.66))
		ant.add_wire(Vector3.new(-1.875,0.035,-2.645), Vector3.new(-2.76,0.665,-2.585))
		return ant
	end
	
	def testc
		ant = Antenna.new("Test C")
	
		ant.add_wire(Vector3.new(0.65987, 2.9107, -0.39789), Vector3.new( 0.68, 2.13, -0.215))
		ant.add_wire(Vector3.new(0.68, 2.13, -0.215), Vector3.new(-0.21371, 1.9233, 1.15837))
		ant.add_wire(Vector3.new(0.415, 2.33, -0.41), Vector3.new( -0.81001, 3.01843, -0.73221))
		return ant
	end
	
	def testd
		ant = Antenna.new("Test D")
	
		ant.add_wire(Vector3.new(0.25708, -0.28144, 0.72579), Vector3.new( 0.23, -0.305, 0.765))
		ant.add_wire(Vector3.new(0.23, -0.305, 0.765), Vector3.new( -0.27268, -0.6793, 1.40084))
		ant.add_wire(Vector3.new(0.0, 0.0, 0.0), Vector3.new( 0.23, -0.305, 0.765))
		ant.add_wire(Vector3.new(0.23, -0.305, 0.765), Vector3.new( 0.25708, -0.28144, 0.72579))
		return ant
	end
	
	def teste
		ant = Antenna.new("Test E")
	
		ant.add_wire(Vector3.new(0.0,0.5,0.0), Vector3.new(1.0,0.5,0.0))
		ant.add_wire(Vector3.new(0.5,0.0,0.0), Vector3.new(0.5,0.4995,0.0))
	
		return ant
	end

	def linden
		ant = Antenna.new("Lindex 7 wire")
		ant.add_wire(Vector3.new(0.0, 0.0, 0.0), Vector3.new(-0.0166, 0.0045, 0.0714), 0.001)
		ant.add_wire(Vector3.new(-0.0166, 0.0045, 0.0714), Vector3.new(-0.0318, -0.0166, 0.0170), 0.001)
		ant.add_wire(Vector3.new(-0.0318, -0.0166, 0.0170), Vector3.new(-0.0318, -0.0287, 0.0775), 0.001)
 		ant.add_wire(Vector3.new(-0.0318, -0.0287, 0.0775), Vector3.new(-0.0318, 0.0439, 0.0140), 0.001)
 		ant.add_wire(Vector3.new(-0.0318, 0.0439, 0.0140), Vector3.new(-0.0318, 0.0045, 0.0624), 0.001)
 		ant.add_wire(Vector3.new(-0.0318, 0.0045, 0.0624), Vector3.new(-0.0106, 0.0378, 0.0866), 0.001)
 		ant.add_wire(Vector3.new(-0.0106, 0.0378, 0.0866), Vector3.new(-0.0106, 0.0257, 0.0230), 0.001)
		ant.freq = 1600.0
		return ant
	end

	def setup
		@a = testa
		@b = testb
		@c = testc
		@d = testd
		@e = teste

		@l = linden
	end
	
	# def teardown
	# end
	
 
	def test_simple
		assert_equal(3, @a.wires.length )
		@a.cleanup
		assert_equal(9, @a.wires.length )
		@a.cleanup 	# Test that another cleanup does nothing
		assert_equal(9, @a.wires.length )

		@e.cleanup
		assert_equal(3, @e.wires.length )
		@e.cleanup	# Test that another cleanup does nothing
		assert_equal(3, @e.wires.length )

		@l.cleanup
		assert_equal(7, @l.wires.length )
		@l.cleanup	# Test that another cleanup does nothing
		assert_equal(7, @l.wires.length )

		np = @l.old_stats(300)
		print "Impedance : #{np.impedance.abs()}\n"
		gm = np.gain_mean
		gsd = np.gain_sd
		print "Gain : #{gm} +- #{gsd}\n"
		rgm = np.gain_rhcp_mean
		rgsd = np.gain_rhcp_sd
		print "RHCP Gain : #{rgm} +- #{rgsd}\n"
		lgm = np.gain_lhcp_mean
		lgsd = np.gain_lhcp_sd
		print "LHCP Gain : #{lgm} +- #{lgsd}\n"
		assert_equal(2.0, sf(np.gain_rhcp_mean,0.5), "Linden 7-wire RHCP gain average" )
		assert_equal(-5.0, sf(np.gain_lhcp_mean,0.5), "Linden 7-wire LHCP gain average" )
		assert_equal(412.0, np.impedance.abs().round )
		
		ns = @l.get_statistics(300)
		print "LHCP Gain 2 : #{ns.gain_lhcp_max} - #{ns.gain_lhcp_min}\n"
		print "RHCP Gain 2 : #{ns.gain_rhcp_max} - #{ns.gain_rhcp_min}\n"
		assert_equal(sf(ns.impedance.abs(),0.01), sf(np.impedance.abs(),0.01), "impedance" )

		assert_equal(sf(ns.gain_max,0.01), sf(np.gain_max,0.01), "gain maximum" )
		assert_equal(sf(ns.gain_rhcp_max,0.01), sf(np.gain_rhcp_max,0.01), "gain rhcp_maximum" )
		assert_equal(sf(ns.gain_lhcp_max,0.01), sf(np.gain_lhcp_max,0.01), "gain lhcp_maximum" )

		assert_equal(sf(ns.gain_mean,0.05), sf(np.gain_mean,0.05), "gain mean" )
		assert_equal(sf(ns.gain_rhcp_mean,0.05), sf(np.gain_rhcp_mean,0.05), "gain rhcp_mean" )
		assert_equal(sf(ns.gain_lhcp_mean,0.05), sf(np.gain_lhcp_mean,0.05), "gain lhcp_mean" )

		assert_equal(sf(ns.gain_sd,0.005), sf(np.gain_sd,0.005), "gain sd" )
		assert_equal(sf(ns.gain_rhcp_sd,0.01), sf(np.gain_rhcp_sd,0.01), "gain rhcp_sd" )
		assert_equal(sf(ns.gain_lhcp_sd,0.01), sf(np.gain_lhcp_sd,0.01), "gain lhcp_sd" )

	end
	
	def test_typecheck
#		assert_raise( RuntimeError ) { SimpleNumber.new('a') }
	end
	
	def test_failure
#		assert_equal(3, SimpleNumber.new(2).add(2), "Adding doesn't work" )
	end
end