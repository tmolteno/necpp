/*
	Electromagnetic Functions and Definitions for nec2++
	
	Copyright (C) 2006  Timothy C.A. Molteno
	
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
#include "electromag.h"


nec_float em::constants::permittivity = 8.854e-12;
nec_float em::constants::permeability = four_pi() * 1.0e-7;

nec_float em::speed_of_light() // was (CVEL in old nec2) but we have changed it to be in meters per second
{
	static nec_float ret = 1.0 / sqrt(permeability() * permittivity());
	return ret; // 299.8e6;
}

nec_float em::impedance()	// was (ETA in old nec2)
{
	static nec_float ret = sqrt(permeability() / permittivity());
	return ret; // 376.8
}

nec_float em::inverse_impedance()	// was (RETA in old nec2)
{
	static nec_float ret = 1.0 / em::impedance();
	return ret; // 2.654420938E-3;
}

nec_float em::impedance_over_2pi()	// was (59.958 in old nec2)
{
	static nec_float ret = em::impedance() / two_pi();
	return ret; // 59.958
}

