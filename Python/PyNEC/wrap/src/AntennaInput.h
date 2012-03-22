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
#ifndef __Antenna_Input__
#define __Antenna_Input__

#include <vector>
#include "BaseInput.h"

/*
                                          - - - ANTENNA INPUT PARAMETERS - - -

   TAG   SEG.    VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)      POWER
   NO.   NO.    REAL        IMAG.       REAL        IMAG.       REAL        IMAG.       REAL        IMAG.     (WATTS)
     0 *   5 1.00000E+00 0.00000E+00 6.64451E-03-3.86651E-03 1.12429E+02 6.54238E+01 6.64451E-03-3.86651E-03 3.32225E-03

                       --------- ANTENNA INPUT PARAMETERS ---------
  TAG   SEG       VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)     POWER
  NO.   NO.     REAL      IMAGINARY     REAL      IMAGINARY     REAL      IMAGINARY    REAL       IMAGINARY   (WATTS)
    0     5  1.0000E+00  0.0000E+00  6.6443E-03 -3.8666E-03  1.1243E+02  6.5428E+01  6.6443E-03 -3.8666E-03  3.3222E-03
*/
class AntennaInput : public BaseInput
{
public:	
	vector<double> tag, seg, vRe, vIm, iRe, iIm, zRe, zIm, power;
	long n_items;

	AntennaInput(std::string& filename)
		: BaseInput(filename)
	{
		n_items = 0;
		string searchString("ANTENNA INPUT PARAMETERS");
		while (m_stream.good())
		{
		string line = readline();

		if (line.find(searchString,0) != string::npos)
		{
			while (line.find("(WATTS)",0) == string::npos)
				line = readline();

			line = readline();

			// get rid of '*' characters
			if (line.find("*",0) != string::npos)
			{
				line.erase(line.find("*",0),1);
			}
			
			stringstream ss(line);
			tag.push_back(read_fixed(ss));
			seg.push_back(read_fixed(ss));
			vRe.push_back(read_fixed(ss));
			vIm.push_back(read_fixed(ss));
			iRe.push_back(read_fixed(ss));
			iIm.push_back(read_fixed(ss));
			zRe.push_back(read_fixed(ss));
			zIm.push_back(read_fixed(ss));
			power.push_back(read_fixed(ss));

			cout << "Impedance : " << zRe[n_items] << " " << zIm[n_items] << endl;
			n_items++;
		}	
		}
	}

	bool equalto(const AntennaInput& ai)
	{
		if (difference(ai) > 1e-4)
			return false;

		return true;
	}

	double difference(const AntennaInput& ai)
	{
		double ret = 0;

		if (n_items != ai.n_items)
			return 1.0;

		try
		{
		for (long i=0; i<n_items; i++)
		{
			ret += diff(ai.vRe[i],vRe[i]);
			ret += diff(ai.vIm[i],vIm[i]);
			ret += diff(ai.iRe[i],iRe[i]);
			ret += diff(ai.iIm[i],iIm[i]);
			ret += diff(ai.zRe[i],zRe[i]);
			ret += diff(ai.zIm[i],zIm[i]);
			ret += diff(ai.power[i],power[i]);
		}
		}
		catch(string message)
		{
			cout << "diff : " << message << endl;
		}
		return ret;
	};
};

#endif /* __Antenna_Input__ */
