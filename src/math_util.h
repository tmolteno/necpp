#ifndef __math_util__
#define __math_util__

/*
  Various Useful Math Utilities for nec2++
  
  Copyright (C) 2004-2015  Timothy C.A. Molteno
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
 
#include "common.h"


/*! \brief      Use the Eigen package for arrays 
    \todo Work through how this should be done.
*/
//#define USING_EIGEN_3VECT 1

#if USING_EIGEN_ARRAY
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wconversion"
    #include <Eigen/Dense>
    #pragma GCC diagnostic pop

  typedef Eigen::Matrix<int32_t, Eigen::Dynamic, 1>  int_array;
  typedef Eigen::Matrix<nec_float, Eigen::Dynamic, 1>  real_array;
  typedef Eigen::Matrix<nec_complex, Eigen::Dynamic, 1>  complex_array;
#else
  // Use our own types rather than Eigen
  #include "safe_array.h"
  typedef safe_array<int32_t>  int_array;
  typedef safe_array<nec_float>  real_array;
  typedef safe_array<nec_complex>  complex_array;

  typedef safe_matrix<int32_t>  int_matrix;
  typedef safe_matrix<nec_float>  real_matrix;
  typedef safe_matrix<nec_complex>  complex_matrix;
#endif


inline void vector_fill(complex_array& x, int64_t start, int64_t N, const nec_complex& y) {
  x.fill(start, N, y);
}


inline nec_complex cplx_00() {
  static nec_complex _cplx00(0.0,0.0); return _cplx00;
}

inline nec_complex cplx_01() {
  static nec_complex _cplx01(0.0,1.0); return _cplx01;
}

inline nec_complex cplx_10() {
  static nec_complex _cplx10(1.0,0.0); return _cplx10;
}

inline nec_complex cplx_11() {
  static nec_complex _cplx11(1.0,1.0); return _cplx11;
}


inline nec_complex cplx_exp(const nec_float& x) {
  return nec_complex(cos(x),sin(x));
}


inline nec_float pi() {
  static nec_float _pi = 3.1415926536; return _pi;
}

inline nec_float two_pi() {
  static nec_float _tmp = 2.0 * pi(); return _tmp;
}

inline nec_float four_pi() {
  static nec_float _tmp = 4.0 * pi(); return _tmp;
}

inline nec_float pi_two() {
  static nec_float _tmp = pi() / 2.0; return _tmp;
}

inline nec_float sqrt_pi() {  // was SP from common.h
  static nec_float _tmp = sqrt(pi()); return _tmp;
}


inline nec_complex two_pi_j() {
  static nec_complex _tmp(0.0,two_pi()); return _tmp;
}



inline nec_float rad_to_degrees(nec_float in_radians) {
  static nec_float _rad_to_deg = 360.0 / (2 * pi()); // 57.29577951
  return in_radians * _rad_to_deg;
}

inline nec_float degrees_to_rad(nec_float in_degrees) {
  static nec_float _deg_to_rad = (2 * pi()) / 360.0;
  return in_degrees * _deg_to_rad;
}

/*! \brief Create a complex number from a magnitude and	an angle in degrees.
*/
inline nec_complex deg_polar(nec_float r, nec_float theta) {
  return std::polar(r, degrees_to_rad(theta));
}


/*! \brief Get the angle of a complex number in degrees.
*/
inline nec_float arg_degrees(nec_complex z) {
  return rad_to_degrees(std::arg(z));
}


/*!\brief arctangent modified to return 0 when x=y=0.
*/
inline nec_float atgn2( nec_float x, nec_float y) {
  if ((0.0 == y) && (0.0 == x))
    return 0.0;
          
  return( std::atan2(y, x) );
}


/*! \brief Decibel dB for magnitude (field) */
inline nec_float db10( nec_float x ) {
  if ( x < 1.0e-20 )
    return( -999.99 );
  
  return( 10.0 * log10(x) );
}

/*! \brief magnitude from Decibel dB */
inline nec_float from_db10( nec_float x ) {
  if ( x < -99.9 )
    return( 0.0 );
  
  return( std::pow(10,x / 10.0));
}


/*! \brief Decibel dB for mag**2 (power) */
inline nec_float db20( nec_float x ) {
  if ( x < 1.0e-20 )
    return( -999.99 );
  
  return( 20.0 * log10(x) );
}


inline nec_float norm(const nec_float x, const nec_float y) {
  return std::sqrt(x*x + y*y);
}

inline nec_float norm2(const nec_float x, const nec_float y, const nec_float z) {
  return (x*x + y*y + z*z);
}

inline nec_float norm(const nec_float x, const nec_float y, const nec_float z) {
  return std::sqrt(norm2(x,y,z));
}

/**!\brief The L1-distance (often called the Manhattan norm) */
inline nec_float normL1(const nec_float x, const nec_float y, const nec_float z) {
  return std::fabs(x) + std::fabs(y) + std::fabs(z);
}




#if USING_EIGEN_3VECT
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wconversion"
    #include <Eigen/Dense>
    #pragma GCC diagnostic pop


    typedef Eigen::Matrix<nec_float,3,1> nec_3vector;
    typedef Eigen::Matrix<nec_complex,3,1> nec_c3vector;
#else

/** \brief A Class for handling 3 dimensional vectors */
class nec_3vector
{
public:
    
  nec_3vector(const nec_float& in_x, const nec_float& in_y, const nec_float& in_z) {
    _v[0] = in_x;
    _v[1] = in_y;
    _v[2] = in_z;
  }
  
  nec_3vector(const nec_3vector&) = default;
  
  /**!\brief The Euclidian norm */
  inline nec_float norm() const {
    return ::norm(_v[0], _v[1], _v[2]);
  }
  
  inline nec_float norm2() const {
    return ::norm2(_v[0], _v[1], _v[2]);
  }
  
  /**!\brief The L1-distance (often called the Manhattan norm) */
  inline nec_float normL1() const {
    return ::normL1(_v[0], _v[1], _v[2]);
  }
  
  inline nec_3vector& operator=(const nec_3vector& copy) {
    _v[0] = copy._v[0];
    _v[1] = copy._v[1];
    _v[2] = copy._v[2];
    return *this;
  }
  
  inline int operator==(const nec_3vector& copy) const {
    if (_v[0] != copy._v[0])
      return 0;
    if (_v[1] != copy._v[1])
      return 0;
    if (_v[2] != copy._v[2])
      return 0;

    return 1;
  }
  
  inline nec_3vector& operator+=(const nec_3vector& a) {
    _v[0] += a._v[0];
    _v[1] += a._v[1];
    _v[2] += a._v[2];
    return *this;
  }
  
  inline nec_3vector operator+(nec_float a) const {
    return nec_3vector(_v[0] + a, _v[1] + a, _v[2] + a);
  }
  
  inline nec_3vector operator+(const nec_3vector& a) const {
    return nec_3vector(_v[0] + a._v[0], _v[1] + a._v[1], _v[2] + a._v[2]);
  }
  
  inline nec_3vector& operator-=(const nec_3vector& a) {
    _v[0] -= a._v[0];
    _v[1] -= a._v[1];
    _v[2] -= a._v[2];
    return *this;
  }
  
  inline nec_3vector operator-(const nec_3vector& a) const {
    return nec_3vector(_v[0] - a._v[0], _v[1] - a._v[1], _v[2] - a._v[2]);
  }
  
  
  inline nec_3vector& operator/=(const nec_float& a) {
    _v[0] /= a;
    _v[1] /= a;
    _v[2] /= a;
    return *this;
  }
  inline nec_3vector operator/(nec_float a) const {
    return nec_3vector(_v[0] / a, _v[1] / a, _v[2] / a);
  }
  
  inline nec_float dot(const nec_3vector& a) const {
    return (_v[0] * a._v[0]) + (_v[1] * a._v[1]) + (_v[2] * a._v[2]);
  }
  
  inline nec_3vector& operator*=(const nec_float& a) {
    _v[0] *= a;
    _v[1] *= a;
    _v[2] *= a;
    return *this;
  }
  inline nec_3vector operator*(nec_float a) const {
    return nec_3vector(_v[0] * a, _v[1] * a, _v[2] * a);
  }
  
  inline nec_float& operator()(int i) {
    return _v[i];
  }
  
  inline const nec_float& operator()(int i) const {
    return _v[i];
  }
  
  /**\brief Cross-product */
  nec_3vector operator*(const nec_3vector& a) const {
    return nec_3vector(
      _v[1]*a._v[2] - _v[2]*a._v[1],
      _v[2]*a._v[0] - _v[0]*a._v[2],
      _v[0]*a._v[1] - _v[1]*a._v[0]);
  }
  
  inline nec_float x() const	{ return _v[0]; }
  inline nec_float y() const	{ return _v[1]; }
  inline nec_float z() const	{ return _v[2]; }
  
  
private:
  nec_float _v[3];
};
#endif

/**!\brief The Euclidian norm */
inline nec_float norm(const nec_3vector& v) {
  return v.norm();
}

/**!\brief The L1-distance (often called the Manhattan norm) */
inline nec_float normL1(const nec_3vector& v) {
  return normL1(v(0), v(1), v(2));
}

#endif /* __math_util__ */
