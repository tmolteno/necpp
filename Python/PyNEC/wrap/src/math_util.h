#ifndef __math_util__
#define __math_util__

/*
	Various Useful Math Utilities for nec2++
	
	Copyright (C) 2004-2005  Timothy C.A. Molteno
	tim@molteno.net 
	
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
 
#include <complex>

#include "safe_array.h"

/*!
	This typedef allows us to use nec2++ with
	a different complex number precision. For example
	float, or long double.
*/
typedef double nec_float;
typedef std::complex<nec_float> nec_complex;

typedef safe_array<int> int_array;
typedef safe_array<nec_float>  real_array;
typedef safe_array<nec_complex>  complex_array;


inline nec_complex cplx_00()
{
	static nec_complex _cplx00(0.0,0.0); return _cplx00;
}

inline nec_complex cplx_01()
{
	static nec_complex _cplx01(0.0,1.0); return _cplx01;
}

inline nec_complex cplx_10()
{
	static nec_complex _cplx10(1.0,0.0); return _cplx10;
}

inline nec_complex cplx_11()
{
	static nec_complex _cplx11(1.0,1.0); return _cplx11;
}


inline nec_complex cplx_exp(nec_float x)
{
	return nec_complex(cos(x),sin(x));
}


inline nec_float pi()
{
	static nec_float _pi = 3.1415926536; return _pi;
}

inline nec_float two_pi()
{
	static nec_float _tmp = 2.0 * pi(); return _tmp;
}

inline nec_float pi_two()
{
	static nec_float _tmp = pi() / 2.0; return _tmp;
}

inline nec_complex two_pi_j()
{
	static nec_complex _tmp(0.0,two_pi()); return _tmp;
}

inline nec_float rad_to_degrees(nec_float in_radians)
{
	static nec_float _rad_to_deg = 360.0 / (2 * pi()); // 57.29577951
	
	return in_radians * _rad_to_deg;
}

inline nec_float degrees_to_rad(nec_float in_degrees)
{
	static nec_float _deg_to_rad = (2 * pi()) / 360.0;
	
	return in_degrees * _deg_to_rad;
}

/**
	Create a complex number from a magnitude and
	an angle in degrees.
*/
inline nec_complex deg_polar(nec_float r, nec_float theta)
{
	return std::polar(r, degrees_to_rad(theta));
}


/**
	Get the angle of a complex number in degrees.
*/
inline nec_float arg_degrees(nec_complex z)
{
	return rad_to_degrees(std::arg(z));
}


/**
	atgn2 is the arctangent function modified to return 0 when x=y=0.
*/
inline nec_float atgn2( nec_float x, nec_float y)
{
	if ((0.0 == y) && (0.0 == x))
		return 0.0;
		
	return( std::atan2(y, x) );
}


/** This function returns db for magnitude (field) */
inline nec_float db10( nec_float x )
{
	if ( x < 1.0e-20 )
		return( -999.99 );
	
	return( 10.0 * log10(x) );
}


/** This function returns db for mag**2 (power) */
inline nec_float db20( nec_float x )
{
	if ( x < 1.0e-20 )
		return( -999.99 );
	
	return( 20.0 * log10(x) );
}


inline nec_float norm(const nec_float x, const nec_float y)
{
	return std::sqrt(x*x + y*y);
}

inline nec_float norm(const nec_float x, const nec_float y, const nec_float z)
{
	return std::sqrt(x*x + y*y + z*z);
}

/**!\brief The L1-distance (often called the Manhattan norm) */
inline nec_float normL1(const nec_float x, const nec_float y, const nec_float z)
{
	return std::fabs(x) + std::fabs(y) + std::fabs(z);
}

inline nec_float norm2(const nec_float x, const nec_float y, const nec_float z)
{
	return (x*x + y*y + z*z);
}

/** \brief A Class for handling 3 dimensional vectors */
class nec_3vector
{
public:
	 
	nec_3vector(const nec_float& in_x, const nec_float& in_y, const nec_float& in_z)
		: m_x(in_x), m_y(in_y), m_z(in_z)
	{ };
	
	nec_3vector(const nec_3vector& copy)
		: m_x(copy.m_x), m_y(copy.m_y), m_z(copy.m_z)
	{ };
	
	/**!\brief The Euclidian norm */
	inline nec_float norm() const
	{
		return ::norm(m_x, m_y, m_z);
	}
	
	inline nec_float norm2() const
	{
		return ::norm2(m_x, m_y, m_z);
	}
	
	/**!\brief The L1-distance (often called the Manhattan norm) */
	inline nec_float normL1() const
	{
		return ::normL1(m_x, m_y, m_z);
	}
	
	inline nec_3vector& operator=(const nec_3vector& copy)
	{
		m_x = copy.m_x;
		m_y = copy.m_y;
		m_z = copy.m_z;
		return *this;
	}
	
	inline int operator==(const nec_3vector& copy) const
	{
		if (m_x != copy.m_x)
			return 0;
		if (m_y != copy.m_y)
			return 0;
		if (m_z != copy.m_z)
			return 0;
		
		return 1;
	}
	
	inline nec_3vector& operator+=(const nec_3vector& a)
	{
		m_x += a.m_x;
		m_y += a.m_y;
		m_z += a.m_z;
		return *this;
	}
	inline nec_3vector operator+(nec_float a) const
	{
		return nec_3vector(m_x + a, m_y + a, m_z + a);
	}
	inline nec_3vector operator+(const nec_3vector& a) const
	{
		return nec_3vector(m_x + a.m_x, m_y + a.m_y, m_z + a.m_z);
	}
	
	inline nec_3vector& operator-=(const nec_3vector& a)
	{
		m_x -= a.m_x;
		m_y -= a.m_y;
		m_z -= a.m_z;
		return *this;
	}
	inline nec_3vector operator-(const nec_3vector& a) const
	{
		return nec_3vector(m_x - a.m_x, m_y - a.m_y, m_z - a.m_z);
	}
	
	
	inline nec_3vector& operator/=(const nec_float& a)
	{
		m_x /= a;
		m_y /= a;
		m_z /= a;
		return *this;
	}
	inline nec_3vector operator/(nec_float a) const
	{
		return nec_3vector(m_x / a, m_y / a, m_z / a);
	}
	
	inline nec_float dot(const nec_3vector& a) const
	{
		return (m_x * a.m_x) + (m_y * a.m_y) + (m_z * a.m_z);
	}
	
	inline nec_3vector& operator*=(const nec_float& a)
	{
		m_x *= a;
		m_y *= a;
		m_z *= a;
		return *this;
	}
	inline nec_3vector operator*(nec_float a) const
	{
		return nec_3vector(m_x * a, m_y * a, m_z * a);
	}
	
	/**\brief Cross-product */
	nec_3vector operator*(const nec_3vector& a) const
	{
		return nec_3vector(
				m_y*a.m_z - m_z*a.m_y,
				m_z*a.m_x - m_x*a.m_z,
				m_x*a.m_y - m_y*a.m_x
			);
	}
	
	inline nec_float x() const	{ return m_x; }
	inline nec_float y() const	{ return m_y; }
	inline nec_float z() const	{ return m_z; }
	
	
private:
	nec_float m_x, m_y, m_z;
};

/**!\brief The Euclidian norm */
inline nec_float norm(const nec_3vector& v)
{
	return v.norm();
}

/**!\brief The L1-distance (often called the Manhattan norm) */
inline nec_float normL1(const nec_3vector& v)
{
	return v.normL1();
}

#endif /* __math_util__ */
