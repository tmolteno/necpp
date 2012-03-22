/*
	Copyright (C) 2008-2011  Timothy C.A. Molteno
	tim@physics.otago.ac.nz
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __nec_wire__
#define __nec_wire__

#include "math_util.h"

/*!\brief A class to handle properties of wires.
*/
class nec_wire
{
public:
	nec_wire(const nec_3vector& a, const nec_3vector& b, nec_float in_radius, int id)
		: x0(a), x1(b), radius(in_radius), _tag_id(id)
	{
	}

	nec_3vector parametrize(nec_float s) const
	{
		nec_3vector ret = x0;
		ret += (x1 - x0)*s;
		return ret;
	}

	nec_float distance(const nec_3vector& a, const nec_3vector& b) const
	{
		return (a - b).norm();
	}

	nec_float length() const
	{
		return distance(x0, x1);
	}
	
	int tag_id() const
	{
		return _tag_id;
	}

/*!\brief Calculate whether two wires intersect.
\return A list of the wires that should be created.
*/
std::vector<nec_wire> intersect(nec_wire& b)
{
	nec_float d2,sa,sb;
	int_solve(x0,x1, b.x0, b.x1, d2,sa,sb);
	
	nec_float epsa = b.radius / length(); // Set by the radius of the other wires
	nec_float epsb = radius / b.length(); // Set by the radius of the other wires
	
	std::vector<nec_wire> ret;
	
	if (d2 > (radius + b.radius)) return ret;
	
	if ((sa >= -epsa) && (sa <= 1.0 + epsa) && (sb >= -epsb) && (sb <= 1.0 + epsb))
	{
		nec_3vector a_pt = parametrize(sa);
		nec_3vector b_pt = b.parametrize(sb);
		
		if (sa > epsa) ret.push_back(nec_wire(x0, a_pt, radius,_tag_id));
		if (sa < 1.0 - epsa) ret.push_back(nec_wire(a_pt, x1, radius,_tag_id));
		
		if (sb > epsb) ret.push_back(nec_wire(b.x0, b_pt, b.radius,b._tag_id));
		if (sb < 1.0 - epsb) ret.push_back(nec_wire(b_pt, b.x1, b.radius, b._tag_id));
	}
	
	return ret;
}


/*!\brief Calculate whether the point is inside the wire
\return true if the point x is inside the wire
*/
bool intersect(nec_3vector& b0)
{
	nec_float a0x = x0.x(); nec_float a0y = x0.y(); nec_float a0z = x0.z();
	nec_float a1x = x1.x(); nec_float a1y = x1.y(); nec_float a1z = x1.z();
	
	nec_float b0x = b0.x(); nec_float b0y = b0.y(); nec_float b0z = b0.z();
	/*
#
#	Do the expression for intersection of a cylinder and a point
#
#	aptitude install python-sympy
#
#
from sympy import *
from sympy.matrices import Matrix

a0x = Symbol('a0x')
a1x = Symbol('a1x')
a0y = Symbol('a0y')
a1y = Symbol('a1y')
a0z = Symbol('a0z')
a1z = Symbol('a1z')

b0x = Symbol('b0x')
b0y = Symbol('b0y')
b0z = Symbol('b0z')

a0 = Matrix([a0x,a0y,a0z])
b0 = Matrix([b0x,b0y,b0z])

a1 = Matrix([a1x,a1y,a1z])

# Each wire is parametrized by sa and sb
sa = Symbol('sa')

a = a0 + sa*(a1 - a0)
b = b0

# Distance between the two wires is d2
delta = (b-a)
d2 = delta.dot(delta)

print "Distance ="
print_python(d2)

# Closest point is the set of parameters (sa,sb) that minimize d2
# subject to the constraint that sa and sb are in the range [0,1]

eqn1 = diff(d2, sa)

equations = [Eq(eqn1,0)]
print "diff(d2,sa) = "
print eqn1


solution = solve(equations,[sa])

print solution
*/
	nec_float sa = (a1x*b0x + a1y*b0y + a1z*b0z - a0x*a1x - a0x*b0x - a0y*a1y - a0y*b0y - a0z*a1z - a0z*b0z + a0x*a0x + a0y*a0y + a0z*a0z)/
	(-2.0*a0x*a1x - 2*a0y*a1y - 2.0*a0z*a1z + a0x*a0x + a0y*a0y + a0z*a0z + a1x*a1x + a1y*a1y + a1z*a1z);
  if (sa < 0) sa = 0.0;
  if (sa > 1.0) sa = 1.0;
  
  nec_3vector a_pt = parametrize(sa);
	
	if (distance(a_pt, b0) > radius) return false;
	return true;
	
}

	bool similar(nec_wire& b)
	{
		// Check if the wires share two endpoints
		nec_float d1 = distance(x0, b.x0);
		nec_float d2 = distance(x0, b.x1);
		nec_float da = std::min(d1,d2);

		nec_float d3 = distance(x1, b.x1);
		nec_float d4 = distance(x1, b.x0);
		nec_float db = std::min(d3,d4);

		if ((std::abs(da) < radius) && (std::abs(db) < radius))
			return true;
		else
			return false;
	}

/*!
	
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
	
*/
	static void int_solve(nec_3vector& a0, nec_3vector& a1,
		nec_3vector& b0, nec_3vector& b1,
		nec_float& distance, nec_float& sa, nec_float& sb)
	{
		nec_float a0x = a0.x(); nec_float a0y = a0.y(); nec_float a0z = a0.z();
		nec_float b0x = b0.x(); nec_float b0y = b0.y(); nec_float b0z = b0.z();
		nec_float a1x = a1.x(); nec_float a1y = a1.y(); nec_float a1z = a1.z();
		nec_float b1x = b1.x(); nec_float b1y = b1.y(); nec_float b1z = b1.z();

		nec_float a01x = (a0x - a1x);
		nec_float a01y = (a0y - a1y);
		nec_float a01z = (a0z - a1z);

		nec_float b01x = (b0x - b1x);
		nec_float b01y = (b0y - b1y);
		nec_float b01z = (b0z - b1z);

		nec_float moda = (a01x*a01x + a01y*a01y + a01z*a01z);
		nec_float modb = (b01x*b01x + b01y*b01y + b01z*b01z);

		nec_float tmp2 = (a01x*b01x + a01y*b01y + a01z*b01z);

		nec_float den = (-4.0*tmp2*tmp2 + 4.0*moda*modb);

		distance = 9.0e9; sa = 2.0; sb = 2.0;
		if (0 == den) return;

		nec_float tmp3 = (-4.0*(a0x*a0x + a0y*a0y + a1x*b0x - a0x*(a1x + b0x) + a1y*b0y - a0y*(a1y + b0y) + a01z*(a0z - b0z))*tmp2 + 4.0*moda*((a0x - b0x)*b01x + (a0y - b0y)*b01y + (a0z - b0z)*b01z));

		sa = (a0x*a01x + a0y*a01y + a0z*a01z - a01x*b0x - a01y*b0y - a01z*b0z - tmp3*tmp2/den)/moda;

		sb = -(tmp3/den);

		nec_float d2 = 	pow((a0x - b0x + sa*(a1x - a0x) - sb*(b1x - b0x)),2) +
				pow((a0z - b0z + sa*(a1z - a0z) - sb*(b1z - b0z)),2) +
				pow((a0y - b0y + sa*(a1y - a0y) - sb*(b1y - b0y)),2);

		distance = std::sqrt(d2);
	}

private:
	nec_3vector x0, x1;
	nec_float radius;
	int _tag_id;
};

#endif /*__nec_wire__*/

