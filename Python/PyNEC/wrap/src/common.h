#ifndef __common__
#define __common__
/*
  Various Definitions for nec2++
  
  Copyright (C) 2004-2015  Timothy C.A. Molteno
  
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

#include "typesafe_stdint.h"

#include <complex>
/*! \brief Change basic type used by nec2++
  This typedef allows us to use nec2++ with
  a different complex number precision. For example
  float, or long double.
*/
typedef double nec_float;
typedef std::complex<nec_float> nec_complex;

/* Version information */
#ifndef nec_build_date
  #define nec_build_date BUILD_DATE
#endif

#ifdef WIN32
#include "../win32/nec2++/config.h"
#else
#include "config.h"
#endif

#ifndef build_version
  #define nec_version VERSION " [" nec_build_date "]"
#else
  #define nec_version build_version " [" nec_build_date "]"
#endif

#define UNUSED(x) {(void)(x);}
/*
  These are some common constants that should be moved into more appropriate locations
*/

#define ACCS	1.E-12
#define	CONST2	4.771341188

#define	SMIN	1.e-3


/**
  0=E VOLTAGE (A),
  1=LINEAR WAVE (B),
  2= R CIRC WAVE (B)
  3=L CIRC WAVE (B),
  4= CURRENT (C),
  5= VOLTAGE DISC.
*/
enum excitation_type
{
  EXCITATION_VOLTAGE = 0,
  EXCITATION_LINEAR = 1,
  EXCITATION_CIRC_RIGHT = 2,
  EXCITATION_CIRC_LEFT = 3,
  EXCITATION_CURRENT = 4,
  EXCITATION_VOLTAGE_DISC = 5
};

#endif /* __common__ */
