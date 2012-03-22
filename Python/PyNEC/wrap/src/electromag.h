#ifndef __electromag__
#define __electromag__

/*
	Various Useful Electromagnetism Utilities for nec2cpp
	
	Copyright (C) 2004  Timothy C.A. Molteno
	
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

// Calculate the power from a voltage and a current.
inline nec_float power(nec_complex voltage, nec_complex current)
{
	return 0.5* real(voltage * conj(current));
}

}

#endif /* __electromag__ */
