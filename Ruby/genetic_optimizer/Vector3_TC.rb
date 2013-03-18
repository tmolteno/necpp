require 'Vector3'
require 'test/unit'
class Vector3_TC < Test::Unit::TestCase
	def setup
		@a = Vector3.new(0,0,0)
		@b = Vector3.new(1,1,0)
		@x = Vector3.new(1,0,0)
		@y = Vector3.new(0,1,0)
		@z = Vector3.new(0,0,1)
		@z2 = Vector3.new(0,0,1)
	end
	
	# def teardown
	# end
	
 
	def test_simple
		assert_equal(0.0, @a.norm )
		assert_equal(1.0, @x.norm )
		assert_equal(0.0, @x.dot(@y) )
		assert_equal(false, @x == @y )
		assert_equal(true, @z == @z2 )
	end
	
	def test_typecheck
#		assert_raise( RuntimeError ) { SimpleNumber.new('a') }
	end
	
	def test_failure
#		assert_equal(3, SimpleNumber.new(2).add(2), "Adding doesn't work" )
	end
	
end

