/*
	Copyright (C) 2004-2008  Timothy C.A. Molteno
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
#include "nec_radiation_pattern.h"
#include "nec_context.h"
#include "c_geometry.h"
#include "nec_exception.h"

int nec_radiation_pattern::get_index(int theta_index, int phi_index) const
{
	if (theta_index >= n_theta) 
		throw new nec_exception("nec_radiation_pattern: Theta index too large");
	if (phi_index >= n_phi) 
		throw new nec_exception("nec_radiation_pattern: Phi index too large");
		
	return phi_index*n_theta + theta_index;
}

nec_radiation_pattern::nec_radiation_pattern(int in_n_theta, int in_n_phi,
	nec_float in_theta_start, nec_float in_phi_start,
	nec_float in_delta_theta, nec_float in_delta_phi,
	nec_float in_range,
	nec_ground& in_ground,
	int in_ifar, nec_float in_wavelength,
	nec_float pinr, nec_float pnlr,
	int in_rp_output_format, int in_rp_normalization, int in_rp_ipd, int in_rp_power_average,
	nec_float in_gnor,
	c_plot_card& in_plot_card)
		: m_ground(in_ground), m_plot_card(in_plot_card)
{
	n_theta = in_n_theta;
	n_phi = in_n_phi;
	
	m_theta_start = in_theta_start;
	m_phi_start = in_phi_start;
	
	delta_theta = in_delta_theta;
	delta_phi = in_delta_phi;
	
	m_range = in_range; // was rfld
	
	m_rp_output_format = in_rp_output_format;
	m_rp_normalization = in_rp_normalization;
	m_rp_ipd = in_rp_ipd;
	m_rp_power_average = in_rp_power_average; // was iavp
	
	m_rp_gnor = in_gnor;
	
	int n_angles = n_theta * n_phi;
	
	_gain.resize(n_angles);
	_power_gain_vert.resize(n_angles);
	_power_gain_horiz.resize(n_angles);
	_power_gain_tot.resize(n_angles);
	_power_gain_rhcp.resize(n_angles);
	_power_gain_lhcp.resize(n_angles);
	_polarization_axial_ratio.resize(n_angles);
	_polarization_tilt.resize(n_angles);
	_polarization_sense_index.resize(n_angles);
	_averaging_scales.resize(n_angles);

	_e_theta.resize(n_angles);
	_e_phi.resize(n_angles);
	_e_r.resize(n_angles);
	
	_ifar = in_ifar;
	_wavelength = in_wavelength;
	_pinr = pinr;
	_pnlr = pnlr;
	
	m_analysis_done = false;
	_maximum_gain = -999.0;
}


/*! \brief Write the analyzed data to a file
*/
void nec_radiation_pattern::write_to_file_aux(ostream& os)
{
	if (false == m_analysis_done)
		throw new nec_exception("Internal Error: Radiation Pattern Analysis not done");
	
	static const char  *hpol[4] = { "LINEAR", "RIGHT ", "LEFT  ", " " };
	static const char  *gain_type[2] = { "----- POWER GAINS ----- ", "--- DIRECTIVE GAINS ---" };
	static const char  *igax[4] = { " MAJOR", " MINOR", " VERTC", " HORIZ" };
	
	int i;
	
	output_helper oh(os,_result_format);
		
	if ( _ifar >= 2) 
	{
		oh.section_start("FAR FIELD GROUND PARAMETERS");
		
		if ( _ifar > 3)
		{
			os << endl;
			os << "                                        RADIAL WIRE GROUND SCREEN" << endl;
			os << "                                        "; oh.int_out(5, m_ground.radial_wire_count); os << " WIRES" << endl;
			os << "                                        WIRE LENGTH= "; oh.real_out(8,2, m_ground.radial_wire_length,false); os << " METERS" << endl;
			os << "                                        WIRE RADIUS= "; oh.real_out(10,3, m_ground.radial_wire_radius); os << " METERS" << endl;
		} /* if ( _ifar > 3) */
		
		if ( _ifar != 4 ) 
		{
			std::string hclif;
			if ( (_ifar == 2) || (_ifar == 5) )
				hclif = "LINEAR";
			if ( (_ifar == 3) || (_ifar == 6) )
				hclif= "CIRCLE";
		
		
			os << endl;
			os << "                                        " << hclif << " CLIFF" << endl;
			os << "                                        EDGE DISTANCE= "; oh.real_out(9,2,m_ground.cliff_edge_distance,false); os << " METERS" << endl;
			os << "                                        HEIGHT= "; oh.real_out(8,2,m_ground.cliff_height,false); os << " METERS" << endl;
			os << "                                        SECOND MEDIUM -" << endl;
			os << "                                        RELATIVE DIELECTRIC CONST.= "; oh.real_out(7,3,m_ground.epsr2, false); os << endl;
			os << "                                        CONDUCTIVITY= "; oh.real_out(10,3,m_ground.sig2,false); os << " MHOS" << endl;
		} /* if ( _ifar != 4 ) */
		
	} /* if ( _ifar >= 2) */

	if ( _ifar == 1)
	{
		oh.section_start("RADIATED FIELDS NEAR GROUND");
		os << "    ------- LOCATION -------     --- E(THETA) ---     ---- E(PHI) ----    --- E(RADIAL) ---" << endl;
		os << "      RHO    PHI        Z           MAG    PHASE         MAG    PHASE        MAG     PHASE" << endl;
		os << "    METERS DEGREES    METERS      VOLTS/M DEGREES      VOLTS/M DEGREES     VOLTS/M  DEGREES" << endl;
	}
	else // _ifar != 1
	{
		oh.section_start("RADIATION PATTERNS");
		
		if ( m_range >= 1.0e-20)
		{
			nec_float exrm = 1.0 / m_range;
			nec_float exra = m_range/ _wavelength;
			exra = -360.0*(exra - floor(exra));
		
			os << "                             RANGE: "; oh.real_out(13,6,m_range); os << " METERS" << endl;
			os << "                             EXP(-JKR)/R: "; oh.real_out(12,5,exrm); os << " AT PHASE: "; oh.real_out(7,2,exra,false); os << " DEGREES" << endl;
		}
		
		int itmp1 = 2 * m_rp_output_format;
		int itmp2 = itmp1+1;
		
		os << " ---- ANGLES -----     "; oh.string_out(23,gain_type[m_rp_ipd]); os << "      ---- POLARIZATION ----   ---- E(THETA) ----    ----- E(PHI) ------" << endl;
		os << "  THETA      PHI      "; oh.string_out(6,igax[itmp1]); os << "  "; oh.string_out(6,igax[itmp2]); os << "   TOTAL       AXIAL      TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE" << endl;
		os << " DEGREES   DEGREES        DB       DB       DB       RATIO   DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES" << endl;
		
	} /* if ( _ifar == 1) */


	i=0;
	nec_float phi = m_phi_start- delta_phi;

	for (int kph = 1; kph <= n_phi; kph++ )
	{
		phi += delta_phi;
		nec_float thet= m_theta_start- delta_theta;
		
		for(int kth = 1; kth <= n_theta; kth++ )
		{
			thet += delta_theta;
			if ( m_ground.present() && (thet > 90.01) && (_ifar != 1) )
				continue;
		
			/* elliptical polarization */
			if ( _ifar == 1)
			{
				nec_complex e_theta = _e_theta[i];
				nec_complex e_phi = _e_phi[i];
				nec_complex e_r = _e_r[i];
				
				oh.start_record();
				oh.padding(" ");
				oh.real_out(9,2,m_range,false);
				oh.separator(); oh.real_out(7,2,phi,false);
				oh.separator(); oh.real_out(9,2,thet,false);
				oh.separator(); oh.real_out(11,4,abs(e_theta));
				oh.separator(); oh.real_out(7,2,arg_degrees(e_theta),false); 
				oh.separator(); oh.real_out(11,4,abs(e_phi));
				oh.separator(); oh.real_out(7,2,arg_degrees(e_phi),false);
				oh.separator(); oh.real_out(11,4,abs(e_r));
				oh.separator(); oh.real_out(7,2,arg_degrees(e_r),false);
				oh.end_record();
			}
			else
			{
				nec_complex e_theta = _e_theta[i];
				nec_complex e_phi = _e_phi[i];
				
				const char* pol_sense = hpol[_polarization_sense_index[i]];
				
				oh.start_record();
				oh.padding(" ");
				oh.real_out(7,2,thet,false); oh.separator(); oh.real_out(9,2,phi,false); oh.separator();
				oh.padding(" ");
				oh.real_out(8,2,_power_gain_vert[i],false);
				oh.separator(); oh.real_out(8,2,_power_gain_horiz[i],false);
				oh.separator(); oh.real_out(8,2,_power_gain_tot[i],false);
				oh.separator(); oh.real_out(11,4,_polarization_axial_ratio[i],false);
				oh.separator(); oh.real_out(9,2,_polarization_tilt[i],false);
				oh.separator(); oh.string_out(6,pol_sense);
				oh.separator(); oh.real_out(11,4,abs(e_theta));
				oh.separator(); oh.real_out(9,2,arg_degrees(e_theta),false);
				oh.separator(); oh.real_out(11,4,abs(e_phi));
				oh.separator(); oh.real_out(9,2,arg_degrees(e_phi),false);
				oh.end_record();
				
				m_plot_card.plot_patterns(thet, phi,
					e_theta, e_phi,
					_power_gain_vert[i], _power_gain_horiz[i], _power_gain_tot[i]);
				
			} /* if ( _ifar != 1) */
		
			i++;
		} /* for( kth = 1; kth <= n_theta; kth++ ) */
	
	} /* for( kph = 1; kph <= n_phi; kph++ ) */

	if ( m_rp_power_average != 0)
	{		
		oh.section_end();
		os << "  AVERAGE POWER GAIN: "; oh.real_out(11,4,_average_power_gain); 
		os << " - SOLID ANGLE USED IN AVERAGING: ("; oh.real_out(7,4,_average_power_solid_angle,false);
		os << ")*PI STERADIANS" << endl;
	}

	if ( m_rp_normalization != 0)
		write_normalized_gain(os);
}


/*! \brief Generate the data for the radiation pattern
*/
void nec_radiation_pattern::analyze(nec_context* m_context)
{
	if (m_analysis_done)
		return;
		
	int pol_sense_index;
	nec_float exrm=0., exra=0., prad, gcon, gcop;
	nec_float phi, thet;
	nec_float tilta, emajr2, eminr2, pol_axial_ratio;
	nec_float dfaz, dfaz2, cdfaz, tstor1=0., tstor2, stilta;
		

	if ((m_context->m_excitation_type == EXCITATION_VOLTAGE) ||
		(m_context->m_excitation_type == EXCITATION_VOLTAGE_DISC) )
	{
		gcop= _wavelength * _wavelength * 2.0 * pi()/(em::impedance() * _pinr);
		prad= _pinr- m_context->structure_power_loss- _pnlr;
		gcon= gcop;
		if ( m_context->ipd != 0)
			gcon= gcon* _pinr/ prad;
	}
	else 
	if ( m_context->m_excitation_type == EXCITATION_CURRENT)
	{
		_pinr=394.51* m_context->xpr6* m_context->xpr6* _wavelength* _wavelength;
		gcop= _wavelength* _wavelength*2.* pi()/(em::impedance() * _pinr);
		prad= _pinr- m_context->structure_power_loss- _pnlr;
		gcon= gcop;
		if ( m_context->ipd != 0)
			gcon= gcon* _pinr/ prad;
	}
	else
	{
		prad=0.;
		gcon=4.* pi()/(1.0+ m_context->xpr6* m_context->xpr6);
		gcop= gcon;
	}

	if ( m_range >= 1.0e-20)
	{
		exrm=1./ m_range;
		exra= m_range/ _wavelength;
		exra=-360.*( exra- floor( exra));
	}
	
	int result_counter=0;
	nec_float pint = 0.0;
	nec_float delta_phi_rad = degrees_to_rad(delta_phi);
	nec_float tmp2 = 0.5 * degrees_to_rad(delta_theta);
	phi= m_phi_start- delta_phi;

	for (int kph = 1; kph <= n_phi; kph++ )
	{
		phi += delta_phi;
		nec_float pha = degrees_to_rad(phi);
		thet= m_theta_start- delta_theta;
		
		for (int kth = 1; kth <= n_theta; kth++ )
		{
			thet += delta_theta;
			if ( m_ground.present() && (thet > 90.01) && (_ifar != 1) )
				continue;
		
			nec_float tha = degrees_to_rad(thet);
			
			nec_complex  eth, eph;
			if ( 1 == _ifar)
			{
				bool space_wave_only = (false == m_ground.present());
				nec_complex erd;
				
				m_context->gfld(m_range/_wavelength, pha, thet/_wavelength,
					&eth, &eph, &erd, space_wave_only, _wavelength );
			
				_e_theta[result_counter] = eth;
				_e_phi[result_counter] = eph;
				_e_r[result_counter] = erd;
			}
			else
			{
				m_context->ffld(tha, pha, &eth, &eph, _wavelength);
		
				nec_float ethm2= norm(eth);
				nec_float ethm= sqrt(ethm2);
				nec_float etha= arg_degrees(eth);
				nec_float ephm2= norm(eph);
				nec_float ephm= sqrt( ephm2);
				nec_float epha= arg_degrees( eph);
		
				/* elliptical polarization calc. */
				if ( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) )
				{
					tilta=0.;
					emajr2=0.;
					eminr2=0.;
					pol_axial_ratio=0.;
					pol_sense_index = 3;
				}
				else
				{
					dfaz= epha- etha;
					if ( epha >= 0.)
						dfaz2= dfaz-360.;
					else
						dfaz2= dfaz+360.;
				
					if ( fabs(dfaz) > fabs(dfaz2) )
						dfaz= dfaz2;
				
					cdfaz= cos(degrees_to_rad(dfaz));
					tstor1= ethm2- ephm2;
					tstor2=2.* ephm* ethm* cdfaz;
					tilta=.5* atan2( tstor2, tstor1);
					stilta= sin( tilta);
					tstor1= tstor1* stilta* stilta;
					tstor2= tstor2* stilta* cos( tilta);
					emajr2=- tstor1+ tstor2+ ethm2;
					eminr2= tstor1- tstor2+ ephm2;
					if ( eminr2 < 0.)
						eminr2=0.;
				
					pol_axial_ratio= sqrt( eminr2/ emajr2);
					tilta= rad_to_degrees(tilta);
					
					if ( pol_axial_ratio <= 1.0e-5)
						pol_sense_index = POL_LINEAR;
					else
					if ( dfaz <= 0.)
						pol_sense_index = POL_RIGHT;
					else
						pol_sense_index = POL_LEFT;
				
				} /* if ( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) ) */
			
				/* Gain Normalization Variables */
				nec_float gnmj= db10( gcon* emajr2);
				nec_float gnmn= db10( gcon* eminr2);
				nec_float gnv = db10( gcon* ethm2);
				nec_float gnh = db10( gcon* ephm2);
				nec_float gtot= db10( gcon*(ethm2+ ephm2) );
			
				if (m_rp_normalization > 0)
				{
					nec_float temp_gain;
					
					switch(m_rp_normalization )
					{
					case 1:
						temp_gain = gnmj;
						break;
				
					case 2:
						temp_gain = gnmn;
						break;
				
					case 3:
						temp_gain = gnv;
						break;
				
					case 4:
						temp_gain = gnh;
						break;
				
					case 5:
						temp_gain = gtot;
						break;
						
					default:
						throw new nec_exception("Unknown Gain Normalization Encountered.");
					}
				
					_gain[result_counter] = temp_gain;
				
				} /* if ( m_rp_normalization > 0) */
			
				if ( m_rp_power_average != 0)
				{
					// compute the numerical integral of the  power gain in angular co-ordinates
					tstor1= gcop*( ethm2+ ephm2);
					nec_float tmp3 = tha - tmp2;
					nec_float tmp4 = tha + tmp2;
				
					if ( kth == 1)
						tmp3= tha;
					else if ( kth == n_theta)
						tmp4= tha;
				
					nec_float da = fabs( delta_phi_rad*( cos(tmp3)- cos(tmp4)));
					if ( (kph == 1) || (kph == n_phi) )
						da *=.5;
					pint += tstor1 * da;
				
					if ( m_rp_power_average == 2) // do not print the power gain values (just compute the average)
						continue;
				}
			
				if ( m_rp_output_format != 1)
				{
					_power_gain_vert[result_counter] = gnmj;
					_power_gain_horiz[result_counter] = gnmn;
				}
				else
				{
					_power_gain_vert[result_counter] = gnv;
					_power_gain_horiz[result_counter] = gnh;
				}
			
				ethm *= _wavelength;
				ephm *= _wavelength;
			
				if ( m_range >= 1.0e-20 )
				{
					ethm= ethm* exrm;
					etha= etha+ exra;
					ephm= ephm* exrm;
					epha= epha+ exra;
				}
			
				_e_theta[result_counter] = deg_polar(ethm, etha);
				_e_phi[result_counter] = deg_polar(ephm, epha);
				
				_power_gain_tot[result_counter] = gtot;
				
				{
					nec_float a = pol_axial_ratio;
					if (pol_sense_index == POL_RIGHT)
						a = -a;
					nec_float gain = from_db10(gtot);
					nec_float lhcp_f =(1+2*a+a*a)/(2*(1+a*a));
					nec_float rhcp_f =(1-2*a+a*a)/(2*(1+a*a));
					_power_gain_rhcp[result_counter] = db10(gain * rhcp_f);
					_power_gain_lhcp[result_counter] = db10(gain * lhcp_f);
				}
				
				_polarization_axial_ratio[result_counter] = pol_axial_ratio;
				_polarization_tilt[result_counter] = tilta;
				_polarization_sense_index[result_counter] = pol_sense_index;
				_averaging_scales[result_counter] = (sin(tha) * (n_theta - 1) / n_theta) + 1.0/n_theta;
				
			} /* if ( _ifar == 1) */
			result_counter++;
		} /* for( kth = 1; kth <= n_theta; kth++ ) */
	
	} /* for( kph = 1; kph <= n_phi; kph++ ) */

	if ( m_rp_power_average != 0)
	{
/*		nec_float tmp3 = degrees_to_rad(m_theta_start);
		tmp4 = tmp3 + degrees_to_rad(delta_theta) * (nec_float)(n_theta-1);
		tmp3 = fabs( degrees_to_rad(delta_phi) * (nec_float)( n_phi-1)*( cos( tmp3)- cos( tmp4)));
		pint /= tmp3;
		tmp3 /= pi();
		_average_power_gain = pint;
		_average_power_solid_angle = tmp3; */

		// We now compute the solid angle over which the power is averaged
		nec_float theta_start_rad = degrees_to_rad(m_theta_start);
		nec_float theta_range = degrees_to_rad(delta_theta) * (nec_float)(n_theta-1);
		nec_float phi_range = degrees_to_rad(delta_phi) * (nec_float)(n_phi-1);
		nec_float total_theta = theta_start_rad + theta_range;
		nec_float solid_angle = fabs(phi_range * (cos(theta_start_rad) - cos(total_theta)) );
						
		_average_power_gain = pint / solid_angle;
		_average_power_solid_angle = solid_angle / pi(); // We display it as a multiple of pi()
	}

	_maximum_gain = _gain.maxCoeff();
	m_analysis_done = true;
}


nec_float nec_radiation_pattern::get_gain_normalization_factor(nec_float gnor)
{
	if ( fabs(gnor) > 1.0e-20)
		return gnor;
		
	if (false == m_analysis_done)
		throw new nec_exception("Internal Error: Radiation Pattern Analysis not done");
		
	return _maximum_gain;
}

void nec_radiation_pattern::write_normalized_gain(ostream& os)
{
	// if  is non-zero then use it as a normaliztion factor
	
	nec_float normalization_factor = get_gain_normalization_factor(m_rp_gnor);
		
	string norm_type;
	switch (m_rp_normalization)
	{
		case 1:
			norm_type = "MAJOR AXIS"; break;
		case 2:
			norm_type = "MINOR AXIS"; break;
		case 3:
			norm_type = "VERTICAL"; break;
		case 4:
			norm_type = "HORIZONTAL"; break;
		case 5:
			norm_type = "TOTAL "; break;
		
		default: throw new nec_exception("Unknown Gain Normalization Encountered.");
	}
	
	output_helper oh(os,_result_format);
	
	oh.section_start("NORMALIZED GAIN");
	{ stringstream mess;	mess << norm_type << " GAIN";	oh.center_text(mess.str(), ""); }
	{ stringstream mess;	mess << "NORMALIZATION FACTOR: " << normalization_factor << " dB"; oh.center_text(mess.str(), ""); }
/*	os << "                                      " << norm_type << " GAIN" << endl;
	os << "                                   NORMALIZATION FACTOR: "; oh.real_out(7,2,normalization_factor,false); os << " db" << endl << endl; */
	os << "    ---- ANGLES ----                ---- ANGLES ----                ---- ANGLES ----" << endl;
	os << "    THETA      PHI        GAIN      THETA      PHI        GAIN      THETA      PHI       GAIN" << endl;
	os << "   DEGREES   DEGREES        DB     DEGREES   DEGREES        DB     DEGREES   DEGREES       DB" << endl;

	
	int row_count = 0;
	int n_cols = 3;
	
	int item_count = 0;
	for (int p=0;p<n_phi;p++)
	{
		nec_float phi = m_phi_start + p*delta_phi ;
		for (int t=0;t<n_theta;t++)
		{

			nec_float theta = m_theta_start + t * delta_theta;
			nec_float norm_gain = _gain[item_count++] - normalization_factor;
			
			oh.start_record();
			oh.padding(" ");
			
			oh.real_out(9,2,theta,false);
			oh.separator(); oh.real_out(9,2,phi,false);
			oh.separator(); oh.real_out(9,2,norm_gain,false);
			oh.padding(" ");
			
			if (_result_format == RESULT_FORMAT_NEC)
			{
				if (item_count % n_cols == 0)
				{
					row_count++;
					oh.end_record();
				}
			}
			else
			{
				oh.end_record();
			}
		}
	}
	os << endl;
}

