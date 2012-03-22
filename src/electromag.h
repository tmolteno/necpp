#ifndef __electromag__
#define __electromag__

/*
	Various Useful Electromagnetism Utilities for nec2cpp
	
	Copyright (C) 2004-2006  Timothy C.A. Molteno
	
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
 
#include <cmath>
#include <complex>

#include "math_util.h"

namespace em
{

class constants
{
public:
	static nec_float permittivity;
	static nec_float permeability;
};

// Calculate the power from a voltage and a current.
inline nec_float power(nec_complex voltage, nec_complex current)
{
	return 0.5* real(voltage * conj(current));
}

// Electromagnetic Constants

inline nec_float permittivity()
{
	return constants::permittivity;
}

inline nec_float permeability()
{
	return constants::permeability;
}

nec_float speed_of_light();
nec_float impedance();
nec_float inverse_impedance();
nec_float impedance_over_2pi();

}

#endif /* __electromag__ */
