/*
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
#ifndef __Radiation_Input__
#define __Radiation_Input__

#include <vector>

#include "BaseInput.h"

class RadiationInput : public BaseInput
{
public:	
	vector<double> theta, phi;
	vector<double>  power_v, power_h, power_t;
	vector<double>  pol_axial_ratio, pol_tilt;
	string pol_sense; // ignore polarization sense
	vector<double> E_theta_mag, E_phi_mag;
	vector<double> E_theta_phase, E_phi_phase;

	long n_items;
	
	RadiationInput(std::string& filename)
		: BaseInput(filename)
	{
		n_items = 0;
		
		string searchString("RADIATION PATTERNS");
		while (m_stream.good())
		{
			string line = readline();
	
			if (line.find(searchString,0) != string::npos)
			{
				
				while (line.find("VOLTS/M",0) == string::npos)
					line = readline();
	
				/*
				- - ANGLES - -           - POWER GAINS -       - - - POLARIZATION - - -    - - - E(THETA) - - -    - - - E(PHI) - - -
				THETA     PHI        VERT.   HOR.    TOTAL      AXIAL     TILT   SENSE     MAGNITUDE    PHASE     MAGNITUDE    PHASE
				DEGREES  DEGREES       DB      DB      DB        RATIO     DEG.              VOLTS/M    DEGREES      VOLTS/M    DEGREES
				   .00      .00    -999.99   -3.23   -3.23     .00000   -90.00  LINEAR    0.00000E+00   -81.39    7.35700E-05  -195.64
				*/
				line = readline();
				
				while (line != "")
				{
					stringstream ss(line);
				
					theta.push_back(read_fixed(ss));
					phi.push_back(read_fixed(ss));
					power_v.push_back(read_fixed(ss)); 
					power_h.push_back(read_fixed(ss));
					power_t.push_back(read_fixed(ss));
					pol_axial_ratio.push_back(read_fixed(ss));
					pol_tilt.push_back(read_fixed(ss));
									
					ss >> pol_sense;
					
					E_theta_mag.push_back(read_sci(ss));
					E_theta_phase.push_back(read_fixed(ss));
	
					E_phi_mag.push_back(read_sci(ss));
					E_phi_phase.push_back(read_fixed(ss));
	
					line = readline();
					n_items++;
				}
				cout << "Radiation pattern: " << n_items << " lines" << endl;
			}	
		}
	}

	bool equalto(const RadiationInput& ai)
	{
		if (difference(ai) > 1e-4)
			return false;

		return true;
	}

	double difference(const RadiationInput& ai)
	{
		double ret = 0.0;
		// compart angles
		if (n_items != ai.n_items)
			return 1;
		
		for (long i=0; i < n_items;i++)
		{
			if (theta[i] != ai.theta[i])
				return 1;
			if (phi[i] != ai.phi[i])
				return 1;
				
			try
			{
				ret += diff(power_v[i], ai.power_v[i]);
				ret += diff(power_h[i], ai.power_h[i]);
				ret += diff(power_t[i], ai.power_t[i]);
				if (power_v[i] > -999.0)
				{
					ret += diff(pol_axial_ratio[i], ai.pol_axial_ratio[i]);
					ret += diff(pol_tilt[i], ai.pol_tilt[i]);
				
					ret += diff(deg_polar(E_theta_mag[i],E_theta_phase[i]), deg_polar(ai.E_theta_mag[i],ai.E_theta_phase[i]));
					ret += diff(deg_polar(E_phi_mag[i],E_phi_phase[i]), deg_polar(ai.E_phi_mag[i],ai.E_phi_phase[i]));
				}
			}
			catch (string message)
			{
				cout << "Diff at [" << theta[i] << "," << phi[i] << "] : " << message << endl;
			}
		}
		return  ret;
	};
};

#endif /* __Radiation_Input__ */
