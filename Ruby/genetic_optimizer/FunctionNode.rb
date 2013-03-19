require 'Vector3.rb'
require 'Antenna.rb'

#
# Base class for construction of programmes (these are trees). So this is 
# the base class for a tree node.
#
class FunctionNode
	attr_accessor :parent, :children, :root
	
	def initialize(root)
		@children = Array.new
		@parent = nil
		@root = root
	end

	def add(c)
		c.parent = self
		@root.add_node(c)
		@children.push(c)
	end

	def remove!(c)
		c.parent = nil
		@root.del_node(c)
		@children.delete(c)
		return c
	end

	def deep_copy(root)
		ret = self.clone
		ret.root = root
		ret.children = Array.new
		@children.each{|c| ret.add(c.deep_copy(root)) }
		return ret
	end

	def random_child
		if (@root.length > 15) then 
			return End.new(@root)
		end

		r = rand(9)
		if (r <= 3)
			child = Wire.new(@root)
			child.set_dr(Vector3.randomize)
			child.mutate
			child.add(child.random_child)
		elsif r == 4
			child = Fork.new(@root)
			child.fork(self.random_child, self.random_child)
		elsif r == 5
			child = Space.new(@root)
			child.set_dr(Vector3.randomize)
			child.mutate
			child.add(child.random_child)
		elsif r == 6
			child = Rotate.new(@root)
			child.rotate(Vector3.random(-3.14, 3.14))
			child.add(child.random_child)
			# 		elsif r == 6
# 			child = CapSMD.new(@root)
# 			child.randomize
# 			child.add(child.random_child)
		else
			child = End.new(@root)
		end
		return child
	end

	def to_s
		ret = "("
		@children.each { |c| ret += "#{c}," }
		ret += ") "
		return ret
	end
	
	def evaluate_children(pos, angle, ant)
		ant.function_count += 1
		@children.each { |c| c.evaluate(pos, angle, ant) }
	end
end

# An End Node
class End < FunctionNode

	def to_s
		"End" + super.to_s
	end
	
	def mutate
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		#		print "end\n"
	end
end

# An offset wire
class Wire < FunctionNode

	def set_dr(dr)
		@dr = dr
	#	@dr.z = 0.00 # do not allow changes in height
	end

	def to_s
		"Wire #{@dr} #{super}"
	end
	
	def mutate
		set_dr(@dr + Vector3.randomize/2)
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		newpos = pos.clone;
		newpos.plus(@dr)
		ant.add_wire(pos,newpos)
#		print "gw #{pos}, #{newpos}\n"
		self.evaluate_children(newpos, angle, ant)
	end
end

# An offset space
class Space < FunctionNode
	def set_dr(dr)
		@dr = dr
		@dr.z = 0.00 # do not allow changes in height
	end

	def to_s
		"Space #{@dr} #{super}"
	end

	def mutate
		set_dr(@dr + Vector3.randomize/2)
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		newpos = pos.clone;
		newpos.plus(@dr)
#		print "move #{pos}, #{newpos}\n"
# 		ant.add_wire(pos,newpos)
		self.evaluate_children(newpos, angle, ant)
	end
end

# A join (fork in the antenna tree)
class Fork < FunctionNode

	def fork(a,b)
		add(a)
		add(b)
	end
	
	def to_s
		"Fork #{super}"
	end

	def mutate
		@children[0].mutate
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		self.evaluate_children(pos, angle, ant)
	end
end

# A rotation about an angle theta
class Rotate < FunctionNode

	def rotate(theta)
		@theta = theta
	end
	
	def to_s
		"Rotate( theta=" + @theta.to_s + "," +  super.to_s
	end

	def mutate
		@theta = @theta + Vector3.random(0,2*Math::PI)
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		self.evaluate_children(pos, angle + @theta, ant)
	end
end

# An 0806 SMD capacitor. Allowed values are
class CapSMD < FunctionNode

	def set_cap(c)
		@c = c
	end

	def to_s
		"CapSMD( c=#{@c}, #{super}"
	end

	def randomize
		@c = 10**Vector3.random(-9.0,-5.0) # 1 pF -> 220 uF
	end

	def mutate
		@c = @c * Vector3.random(0.0,2.0)
	end

	def evaluate(pos, angle, ant)
		ant.function_count += 1
		#		ant.add_impedance(pos, 0.0, 0.0, @c)
		self.evaluate_children(pos, angle, ant)
	end
end
