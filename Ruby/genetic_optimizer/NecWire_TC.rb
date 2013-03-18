require 'NecWire'
require 'test/unit'
class NecWire_TC < Test::Unit::TestCase
	def setup
		@a = NecWire.new(Vector3.new(0,0,0), Vector3.new(1,1,0))
		@b = NecWire.new(Vector3.new(0.3,0,0), Vector3.new(0.3,1,0))
		@c = NecWire.new(Vector3.new(1,1,0), Vector3.new(0.3,1,0))
		@d = NecWire.new(Vector3.new(-0.2178,     0.0585,    -0.2387), Vector3.new(-0.2000,     0.0950,    -0.2650))
		@e = NecWire.new(Vector3.new(-0.2000,     0.0950,    -0.2650), Vector3.new(-0.2178,     0.0585,    -0.2387))

		@l1 = NecWire.new(Vector3.new(0.0, 0.0, 0.0), Vector3.new(-0.0166, 0.0045, 0.0714),0.001)
		@l2 = NecWire.new(Vector3.new(-0.0318, -0.0166, 0.0170), Vector3.new(-0.0318, -0.0287, 0.0775),0.001)

		@ap = NecWire.new(Vector3.new(1,0,2), Vector3.new(0,1,2))
	end
	
	# def teardown
	# end
	
 
	def test_simple
		assert_nil(@l1.intersect(@l2) )

		assert_equal(false, @a.similar(@b) )
		assert_equal(false, @b.similar(@a) )
		assert_equal(false, @b.similar(@c) )
		assert_equal(true, @d.similar(@e) )
		assert_not_nil(@a.intersect(@b) )
		assert_equal(nil, @d.intersect(@e) )

		d2, sa, sb = NecWire.int_solve(@a.x0,@a.x1,@ap.x0,@ap.x1)
		assert_equal(2.0, d2)
		assert_equal(0.5, sa)
		assert_equal(0.5, sb)
		
	end
	
	def test_typecheck
#		assert_raise( RuntimeError ) { SimpleNumber.new('a') }
	end
	
	def test_failure
#		assert_equal(3, SimpleNumber.new(2).add(2), "Adding doesn't work" )
	end
	
end

