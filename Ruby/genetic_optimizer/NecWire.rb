require 'Vector3.rb'
require 'necpp'

require 'Parameters'

#
# A wire segment (this roughly corresponds to an nec2++ GW card
#
class NecWire
	attr_accessor :x0, :x1, :radius

 	def initialize(x0, x1, radius = 0.0005)
		@x0 = x0
		@x1 = x1
		@radius = radius
	end

	#
	# Truncate a number to four significant figures
	#
	def sf4(x)
		return sf(x,@radius / 10)
	end

	def length
		return @x0.distance(@x1)
	end

	def card(tag, wavelength)
		d = length
		oddseg = (((SEGMENTS_PER_WAVELENGTH / 2.0)*d/wavelength).round)*2 + 1
		segments = [3,oddseg].max
		return "GW #{tag} #{segments} #{sf4(@x0.x)} #{sf4(@x0.y)} #{sf4(@x0.z)} #{sf4(@x1.x)} #{sf4(@x1.y)} #{sf4(@x1.z)} #{@radius}"
	end

	def nec_wire(nec,tag, wavelength)
		d = length
		oddseg = (((SEGMENTS_PER_WAVELENGTH / 2.0)*d/wavelength).round)*2 + 1
		segments = [3,oddseg].max
		Necpp.nec_wire(nec, tag, segments, sf4(@x0.x), sf4(@x0.y), sf4(@x0.z), sf4(@x1.x), sf4(@x1.y), sf4(@x1.z), @radius, 1.0, 1.0)
#		print "Necpp.nec_wire(nec, #{tag}, #{segments}, #{@x0.x}, #{@x0.y}, #{@x0.z}, #{@x1.x}, #{@x1.y}, #{@x1.z}, #{@radius}, 1.0, 1.0)\n"
	end

	def to_s
		"Vector3.new(#{@x0.x},#{@x0.y},#{@x0.z}), Vector3.new(#{@x1.x},#{@x1.y},#{@x1.z})"
	end

	# Parametrize the wire in terms of arc length s
	def parametrize(s)
		ret = @x0
		ret += (@x1 - @x0)*s
		return ret
	end

	# Calculate whether self intersects with the wire b
	def intersect(b)
		d2,sa,sb = NecWire.int_solve(@x0,@x1, b.x0, b.x1)
#		print "d2,sa,sb #{d2}, #{sa} #{sb}\n"
		epsa = b.radius / length # Set by the radius of the other wires
		epsb = @radius / b.length # Set by the radius of the other wires

		return nil if (d2 > (@radius + b.radius))
	
		if ((sa >= -epsa) && (sa <= 1.0 + epsa) && (sb >= -epsb) && (sb <= 1.0 + epsb)) then
			a_pt = self.parametrize(sa)
			b_pt = b.parametrize(sb)

			ret = Array.new
			
			ret.push(NecWire.new(@x0,a_pt,@radius)) if (sa > epsa)
			ret.push(NecWire.new(a_pt,@x1,@radius)) if (sa < 1.0 - epsa)

			ret.push(NecWire.new(b.x0,b_pt,b.radius)) if (sb > epsb)
			ret.push(NecWire.new(b_pt,b.x1,b.radius)) if (sb < 1.0 - epsb)
			return ret
		else
			return nil
		end
	end

	def similar(b)
		# Check if the wires share two endpoints
		d1 = @x0.distance(b.x0)
		d2 = @x0.distance(b.x1)
		da = [d1,d2].min

		d3 = @x1.distance(b.x1)
		d4 = @x1.distance(b.x0)
		db = [d3,d4].min

 #		print "similar [#{self.to_s}], [#{b}] = #{da} #{db}\n"
		if ((da.abs < 1e-3) && (db.abs < 1e-3))
			return true
		else
			return false
		end
	end

=begin
	
	We use the following Mathematica expression to get the solution for the intersection
	of two cylinders. We set up the distance between the two center lines as d2
	
	The equations are derived from sympy in the file intersections.py
	
	d2 =    (a0x - b0x + sa*(a1x - a0x) - sb*(b1x - b0x))^2 + 
		(a0z - b0z + sa*(a1z - a0z) - sb*(b1z - b0z))^2 + 
		(a0y - b0y + sa*(a1y - a0y) - sb*(b1y - b0y))^2;

	soln = Solve[{D[d2,sa] == 0, D[d2,sb] == 0}, {sa,sb}];
	Simplify[soln]
	
	sa -> (2*a0x*(a0x - a1x) + 2*a0y*(a0y - a1y) + 2*a0z*(a0z - a1z) + 2*(-a0x \
	+ a1x)*b0x + 2*(-a0y + a1y)*b0y + 2*(-a0z + a1z)*b0z + ((-4*(a0x^2 + a0y^2 + \
	a1x*b0x - a0x*(a1x + b0x) + a1y*b0y - a0y*(a1y + b0y) + (a0z - a1z)*(a0z - \
	b0z))*((a0x - a1x)*(b0x - b1x) + (a0y - a1y)*(b0y - b1y) + (a0z - a1z)*(b0z \
	- b1z)) + 4*((a0x - a1x)^2 + (a0y - a1y)^2 + (a0z - a1z)^2)*((a0x - \
	b0x)*(b0x - b1x) + (a0y - b0y)*(b0y - b1y) + (a0z - b0z)*(b0z - \
	b1z)))*(-2*(a0x - a1x)*(b0x - b1x) - 2*(a0y - a1y)*(b0y - b1y) - 2*(a0z - \
	a1z)*(b0z - b1z)))/(-4*((a0x - a1x)*(b0x - b1x) + (a0y - a1y)*(b0y - b1y) + \
	(a0z - a1z)*(b0z - b1z))^2 + 4*((a0x - a1x)^2 + (a0y - a1y)^2 + (a0z - \
	a1z)^2)*((b0x - b1x)^2 + (b0y - b1y)^2 + (b0z - b1z)^2)))/(2*((a0x - a1x)^2 \
	+ (a0y - a1y)^2 + (a0z - a1z)^2))
	
	sb -> -((-4*(a0x^2 + a0y^2 + a1x*b0x - \
	a0x*(a1x + b0x) + a1y*b0y - a0y*(a1y + b0y) + (a0z - a1z)*(a0z - b0z))*((a0x \
	- a1x)*(b0x - b1x) + (a0y - a1y)*(b0y - b1y) + (a0z - a1z)*(b0z - b1z)) + \
	4*((a0x - a1x)^2 + (a0y - a1y)^2 + (a0z - a1z)^2)*((a0x - b0x)*(b0x - b1x) + \
	(a0y - b0y)*(b0y - b1y) + (a0z - b0z)*(b0z - b1z)))/(-4*((a0x - a1x)*(b0x - \
	b1x) + (a0y - a1y)*(b0y - b1y) + (a0z - a1z)*(b0z - b1z))^2 + 4*((a0x - \
	a1x)^2 + (a0y - a1y)^2 + (a0z - a1z)^2)*((b0x - b1x)^2 + (b0y - b1y)^2 + \
	(b0z - b1z)^2)))
	
=end
	def NecWire.int_solve(a0,a1,b0,b1)
		a0x = a0.x; a0y = a0.y; a0z = a0.z
		b0x = b0.x; b0y = b0.y; b0z = b0.z
		a1x = a1.x; a1y = a1.y; a1z = a1.z
		b1x = b1.x; b1y = b1.y; b1z = b1.z

		a01x = (a0x - a1x);
		a01y = (a0y - a1y);
		a01z = (a0z - a1z);

		b01x = (b0x - b1x);
		b01y = (b0y - b1y);
		b01z = (b0z - b1z);

		moda = (a01x*a01x + a01y*a01y + a01z*a01z);
		modb = (b01x*b01x + b01y*b01y + b01z*b01z);

		tmp2 = (a01x*b01x + a01y*b01y + a01z*b01z);

		den = (-4.0*tmp2*tmp2 + 4.0*moda*modb);

		return 9.0e9, 2, 2 if den == 0
		
		tmp3 = (-4.0*(a0x*a0x + a0y*a0y + a1x*b0x - a0x*(a1x + b0x) + a1y*b0y - a0y*(a1y + b0y) + a01z*(a0z - b0z))*tmp2 + 4.0*moda*((a0x - b0x)*b01x + (a0y - b0y)*b01y + (a0z - b0z)*b01z));

		sa = (a0x*a01x + a0y*a01y + a0z*a01z - a01x*b0x - a01y*b0y - a01z*b0z - tmp3*tmp2/den)/moda;

		sb = -(tmp3/den);

		d2 = (a0x - b0x + sa*(a1x - a0x) - sb*(b1x - b0x))**2 + (a0z - b0z + sa*(a1z - a0z) - sb*(b1z - b0z))**2 + (a0y - b0y + sa*(a1y - a0y) - sb*(b1y - b0y))**2

		return Math.sqrt(d2),sa,sb
	end
end
