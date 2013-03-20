require 'Vector3.rb'
require 'Antenna.rb'
require 'FunctionNode.rb'
require 'drb/drb'

#
# The base class for a member of the GP population. Individual objects
# are aggregated into a Population object.
#
# An indivitual is essentially a Tree of FunctionNode objects. 
#
class Individual

	include DRb::DRbUndumped

	attr_accessor :function_tree, :nodelist, :my_fitness
	

	# Create a random individual
	def initialize
		@function_tree = FunctionNode.new(self)
		self.init_nodelist
		fitness_refresh
	end

	def child(n)
		@nodelist[n]
	end

	def add_node(n)
		@nodelist.push(n)
	end

	def del_node(n)
		@nodelist.delete(n)
	end

	def length
		return @nodelist.length
	end

	def init_nodelist
		@nodelist = Array.new
	end

	def duplicate
		ret = self.clone
		ret.init_nodelist
		ret.function_tree = @function_tree.deep_copy(ret)
		ret.my_fitness = @my_fitness
		#ret.my_stats = nil #@my_stats
		#ret.my_antenna = nil #@my_antenna
		return ret
	end

	def randomize
		@function_tree.add(@function_tree.random_child)
	end

	def to_s
		"[" + @function_tree.to_s + "]\n"
	end
	
	def <=>(a)
		fa = self.fitness
		fb = a.fitness
		fa <=> fb
	end

	def get_antenna
		return @my_antenna unless @my_antenna == nil
		ant = Antenna.new( @function_tree.to_s )
		pos = Vector3.new(0,0,ANT_HEIGHT) # starting position
		@function_tree.evaluate_children(pos, 0.0, ant)
		@my_antenna = ant
		@my_antenna.cleanup
		return @my_antenna
	end

	def set_stats(s)
		@my_stats = s
	end
	
	def fitness_refresh
		@my_fitness = nil
		@my_stats = nil
		@my_antenna = nil
	end

	def fitness
		return @my_fitness unless @my_fitness == nil
		
		self.fitness_calc
		return @my_fitness
	end

	def fitness_display
		fitness_refresh
		self.fitness_calc
		ant = self.get_antenna
		ret = fitness_from_stats(true)
		ant.add_comment(ret)
		return ret
	end

	def Individual::get_stats(ant)
		ret = Array.new
		ret.push(ant.get_statistics(FREQUENCY))
		ret.push(ant.get_statistics(FREQUENCY - BANDWIDTH/2))
		ret.push(ant.get_statistics(FREQUENCY + BANDWIDTH/2))
		return ret
	end
	
	def fitness_calc
		if @my_stats == nil
			ant = self.get_antenna
			@my_stats = Individual::get_stats(ant)
		end

		@my_fitness = fitness_from_stats
	end

	#
	#	Return the fitness from the antenna statistics
	#
	def fitness_from_stats(display=false)
		ret = 9999; wire_count = 0
		return -9999 if @my_stats == nil
		ret_stats = nil
		@my_stats.each { |s| 
			return -99.0 if s == nil
			return -99.0 if s.wire_count == 0
			return -99.0 if s.wire_count > 15

			imp =  s.impedance
			vswr = s.vswr(50.0)
			gain_min = s.gain_min
			glmax = s.gain_max
			gsd = s.gain_sd

			return -99.0 if s.gain_mean > GAIN_MEAN_LIMIT # 
			return -99.0 if vswr == nil
			return -99.0 if gain_min == nil
			return -99.0 if gain_min < -999.0
			return -99.0 if glmax == nil
			return -99.0 if glmax < -999.0
			return -99.0 if imp.abs() < 1.0

			transmission_coeff = 1-s.reflection_coeff(50.0).abs() # 50 Ohm impedance

		  	test_fit = to_db(transmission_coeff) # Maximize transmission coefficient
		        test_fit += gain_min # Maximize RHCP gain
		        test_fit -= gsd # Minimize LHCP gain

		      
			ant = get_antenna
			penalty = 0.0
			ant.wires.each do |w|
		    penalty += 100 * [0.0, w.x0.norm - ANT_MAX_X].max
		    penalty += 100 * [0.0, w.x1.norm - ANT_MAX_Y].max
		  end
			wire_count = s.wire_count
			penalty += wire_count / 15.0;
			
		  test_fit -= penalty*0.2
			
		  if (test_fit < ret)
		     ret = test_fit
		     ret_stats = s
		  end
	  }
		fc = @my_antenna.function_count
		ret -= fc / 30
		
		return "Fit=#{sf(ret,0.001)} #{ret_stats}" if display
		return ret
	end

	def random_node
		node_num = rand(self.length)
		return self.child(node_num)
	end
	
	def mutate
		node = self.random_node
		node.mutate
		self.fitness_refresh
	end

	def Individual.breed(a,b)
		achild = a.random_node
		bchild = b.random_node

		aparent = achild.parent
		bparent = bchild.parent

		aparent.remove!(achild)
		bparent.remove!(bchild)

		aparent.add(bchild)
		bparent.add(achild)
		
		a.fitness_refresh
		b.fitness_refresh
	
		return a,b
	end
end