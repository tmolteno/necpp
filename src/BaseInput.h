/*
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
#pragma once

#include "math_util.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

/* Relative difference of two values, normalised by the larger magnitude.
   The scale test must use magnitudes: comparing signed values against the
   1e-8 floor treats any two negative numbers as "effectively zero". */
nec_float diff(nec_float a, nec_float b)
{
	if (a == b)
		return 0;

	nec_float scale = std::max(std::fabs(a), std::fabs(b));
	if (scale < 1e-8)
		return 0;

	nec_float ret = std::fabs(a - b) / scale;
	if (ret > 1e-2)
	{
		cout << "    diff(" << a << "," << b << ") = " << ret << endl;
	}
	return ret;
}

nec_float diff(nec_complex a, nec_complex b)
{
	if (a == b)
		return 0;

	nec_float scale = std::max(abs(a), abs(b));
	if (scale < 1e-8)
		return 0;

	nec_float ret = abs(a - b) / scale;
	if (ret > 1e-2)
	{
		stringstream ss;
		cout << "    diff(" << a << "," << b << ") = " << ret << endl;
	}
	return ret;
}

class BaseInput
{
public:
	BaseInput(const std::string& filename)
		:	m_filename(filename),
			m_stream(filename.c_str())
	{
	}
	
private:
	char linec[512];


protected:

	std::string m_filename;
	std::ifstream m_stream;
	
	string readline()
	{
		m_stream.getline(&linec[0],512);
		return string(&linec[0]);
	}

	double read_sci(std::istream& is)
	{
		double x = 0.0;	/* zero on extraction failure (was returned uninitialized) */
		is.setf(ios_base::skipws);
		is >> x;
		return x;
	}

	double read_fixed(std::istream& is)
	{
		double x = 0.0;	/* zero on extraction failure (was returned uninitialized) */
		is.setf(ios_base::skipws);
		is >> x;
		return x;
	}
};
