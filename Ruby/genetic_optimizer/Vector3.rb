
def sf(x,del)
	return x if del == 0
	return (x / del).round * del
end

# A 3-vector used for managing geometry on Antennas
class Vector3
	attr_accessor :x, :y, :z
	
	def initialize(x,y,z)
		@x = x
		@y = y
		@z = z
	end
	
	# Construct a random vector
	def Vector3.randomize
		ret = Vector3.new(random(-RANDOM_RANGE,RANDOM_RANGE),random(-RANDOM_RANGE,RANDOM_RANGE),random(-RANDOM_RANGE,RANDOM_RANGE))
		ret.round(0.0001)
		return ret
	end
	
	# In place addition TODO use +=?
	def plus(v)
		@x = @x + v.x
		@y = @y + v.y
		@z = @z + v.z
	end
	
	# Distance between two vectors (points) |a - b|
	def distance(v)
		dx = (@x - v.x)
		dy = (@y - v.y)
		dz = (@z - v.z)
		return Math.sqrt(dx*dx + dy*dy + dz*dz)
	end
	
	# Equality test
	def ==(v)
		return ((@x == v.x) && (@y == v.y) && (@z == v.z))
	end
	
	# Addition of another vector
	def +(v)
		return Vector3.new(@x + v.x,@y + v.y, @z + v.z)
	end

	# Subtraction of another vector
	def -(v)
		return Vector3.new(@x - v.x,@y - v.y, @z - v.z)
	end

	# Devision by a scalar
	def /(x)
		return Vector3.new(@x/x,@y/x, @z/x)
	end

	# Multiplication by a scalar
	def *(x)
		return Vector3.new(@x*x,@y*x, @z*x)
	end

	# Dot product of self with another vector. Returns a scalar.
	def dot(v)
		return @x*v.x + @y*v.y + @z*v.z
	end

	# Euclidean norm
	def norm
		return Math.sqrt(@x*@x + @y*@y + @z*@z)
	end

	# Round to the level of precision del (e.g. 0.01)
	def round(del)
		@x = sf(@x,del)
		@y = sf(@y,del)
		@z = sf(@z,del)
	end
	
	def to_s
		"[#{@x}  #{@y}  #{@z}]"
	end

	# Return a random vector that lies on the line between two vectors x0 and x1.
	def Vector3.random(x0,x1)
		return x0 + (x1 - x0)* (rand(1000000)/1000000.0)
	end
	
end
