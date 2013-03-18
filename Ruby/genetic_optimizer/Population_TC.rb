require 'Population'
require 'Gp.rb'
require 'test/unit'
class Population_TC < Test::Unit::TestCase
	def setup
		@p_even = Population.new(30, EvalQueue.new)
		@p_even.randomize
		@p_odd = Population.new(31, EvalQueue.new)
		@p_odd.randomize
	end
	
	# def teardown
	# end
	
 
	def test_simple
		assert_equal(30, @p_even.individuals.length )
		assert_equal(31, @p_odd.individuals.length )

		best_even = @p_even.evolve(0.1, 0.1)
		assert_not_nil(best_even.fitness )
		best_odd = @p_odd.evolve(0.1, 0.1)
		assert_not_nil(best_odd.fitness )

		assert_equal(30, @p_even.individuals.length )
		assert_equal(31, @p_odd.individuals.length )

	end
	
	def test_typecheck
#		assert_raise( RuntimeError ) { SimpleNumber.new('a') }
	end
	
	def test_failure
#		assert_equal(3, SimpleNumber.new(2).add(2), "Adding doesn't work" )
	end
	
end

