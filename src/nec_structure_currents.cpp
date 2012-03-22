/*
	Copyright (C) 2004-2008  Timothy C.A. Molteno
	
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
#include "nec_structure_currents.h"
#include "nec_context.h"
#include "nec_exception.h"
#include "c_geometry.h"

int nec_structure_currents::get_n()
{
	return m_geometry->n;
}

int nec_structure_currents::get_m()
{
	return m_geometry->m;
}
			

nec_structure_currents::nec_structure_currents(nec_context * in_context, enum excitation_type in_pattype,
	int in_nload,
	nec_float in_xpr3, nec_float in_xpr6)
{
	m_context = in_context;
	m_geometry = m_context->m_geometry;
	pattype = in_pattype;
	
	iptflg = m_context->iptflg;
	iptflq = m_context->iptflq;
	iptag = m_context->iptag;
	iptagf = m_context->iptagf;
	iptagt = m_context->iptagt;
	iptaq = m_context->iptaq;
	iptaqf = m_context->iptaqf;
	iptaqt = m_context->iptaqt;
	
	nload = in_nload;
		
	xpr3 = in_xpr3;
	xpr6 = in_xpr6;
	
	wavelength = m_context->wavelength;
	freq_mhz = m_context->freq_mhz;
	
	structure_power_loss=0;
		
	current_nb_elements = 0;
	q_density_nb_elements = 0;
	q_density_last_printed = 0;
	patch_nb_elements = 0;
}

/*
	EXCITATION_VOLTAGE = 0,
	EXCITATION_LINEAR = 1,
	EXCITATION_CIRC_RIGHT = 2,
	EXCITATION_CIRC_LEFT = 3,
	EXCITATION_CURRENT = 4,
	EXCITATION_VOLTAGE_DISC = 5
*/
std::string nec_structure_currents::hpol(enum excitation_type e)
{
	switch(e)
	{
		case EXCITATION_LINEAR:		return std::string("LINEAR");
		case EXCITATION_CIRC_RIGHT:	return std::string("RIGHT");
		case EXCITATION_CIRC_LEFT:	return std::string("LEFT");
		default:
		{	nec_exception* nex = new nec_exception("Unknown Excitation type");
			throw nex;
		}
	}
}

void nec_structure_currents::analyze()
{
	int jump;
	nec_float cmag;
//	nec_complex curi;
	nec_float fr;
	nec_complex eth, eph, ex, ey, ez;
	
	if (m_geometry->n != 0)
	{
		if (iptflg!= -1)
		{
			int itmp1=0;
			jump= iptflg+1;
	
			for (int i = 0; i < m_geometry->n; i++ )
			{
				nec_complex curi= m_context->current_vector[i]* wavelength;
				cmag= abs(curi);
		
				if ( (nload != 0) && (fabs(real(m_context->zarray[i])) >= 1.e-20) )
					structure_power_loss += 0.5*cmag*cmag*real( m_context->zarray[i]) * m_geometry->segment_length[i];
	
				if ( jump == 0)
				continue;

				if ( jump > 0 )
				{
					if ( (iptag != 0) && (m_geometry->segment_tags[i] != iptag) )
						continue;

					itmp1++;
					if ( (itmp1 < iptagf) || (itmp1 > iptagt) )
						continue;

					if ( iptflg != 0)
					{
						if ( iptflg >= 2 )
						{
							m_context->fnorm[m_context->get_inc()-1]= cmag;
							m_context->set_isave(i+1);
						}

						if ( iptflg != 3)
						{
							/*m_output.nec_printf("          %7.2f  %7.2f   %11.4E  %7.2f  %5d",
								xpr1, xpr2, cmag, ph, i+1 );*/
							
							current_nb_elements++;						
							_current_theta.push_back(m_context->get_xpr1());
							_current_phi.push_back(m_context->get_xpr2());
							_current.push_back(curi);
							_current_segment_number.push_back(i+1);
						
							continue;										
						}
					} /* if ( iptflg != 0) */
			
					else /* iptflg == 0, only the currents specified will be printed, using the standard format*/
					{
						/*m_output.nec_printf(
							" %5d %4d %9.4f %9.4f %9.4f %9.5f"
							" %11.4E %11.4E %11.4E %8.3f",
							i+1, m_geometry->segment_tags[i],
							m_geometry->x[i], m_geometry->y[i], m_geometry->z[i], m_geometry->segment_length[i],
							real(curi), imag(curi), cmag, ph );*/
			
						current_nb_elements++;
						_current_segment_number.push_back(i+1);
						_current_segment_tag.push_back(m_geometry->segment_tags[i]);
						_current_segment_center_x.push_back(m_geometry->x[i]);
						_current_segment_center_y.push_back(m_geometry->y[i]);
						_current_segment_center_z.push_back(m_geometry->z[i]);
						_current_segment_length.push_back(m_geometry->segment_length[i]);
						_current.push_back(curi);
										
						// added test for plot_card.is_valid()
						if (m_context->plot_card.is_valid() && m_context->plot_card.currents())
						{
							m_context->plot_card.plot_complex(curi);
							m_context->plot_card.plot_endl();
						}
					} /* iptflg == 0, only the currents specified will be printed, using the standard format*/
				}
			
				else /* iptflg == -2, all currents will be printed, using the standard format*/
				{
					/*m_output.nec_printf(
						" %5d %4d %9.4f %9.4f %9.4f %9.5f"
						" %11.4E %11.4E %11.4E %8.3f",
						i+1, m_geometry->segment_tags[i],
						m_geometry->x[i], m_geometry->y[i], m_geometry->z[i], m_geometry->segment_length[i],
						real(curi), imag(curi), cmag, ph );*/
			
					current_nb_elements++;
					_current_segment_number.push_back(i+1);
					_current_segment_tag.push_back(m_geometry->segment_tags[i]);
					_current_segment_center_x.push_back(m_geometry->x[i]);	
					_current_segment_center_y.push_back(m_geometry->y[i]);
					_current_segment_center_z.push_back(m_geometry->z[i]);
					_current_segment_length.push_back(m_geometry->segment_length[i]);
					_current.push_back(curi);
					
					// added test for plot_card.is_valid()
					if (m_context->plot_card.is_valid() && m_context->plot_card.currents())
					{
						m_context->plot_card.plot_complex(curi);
						m_context->plot_card.plot_endl();
					}
				} /* iptflg == -2, all currents will be printed, using the standard format*/			

			} /* for( i = 0; i < n; i++ ) */
		
		m_context->structure_power_loss = structure_power_loss;		
		
		}/* if (iptflg != -1) */
	
		if (iptflq != -1)
		{
			int itmp1 = 0;
			fr = 1.e-6/(freq_mhz);

			for(int i = 0; i < m_geometry->n; i++ )
			{
				if ( iptflq != -2 )
				{
					if ( (iptaq != 0) && (m_geometry->segment_tags[i] != iptaq) )
						continue;
	
					itmp1++;
					if ( (itmp1 < iptaqf) || (itmp1 > iptaqt) )
						continue;

				} /* if ( iptflq == -2) */
	
				nec_complex curi = fr * nec_complex(- m_context->bii[i], m_context->bir[i]);
			
				/*m_output.nec_printf(
					" %5d %4d %9.4f %9.4f %9.4f %9.5f"
					" %11.4E %11.4E %11.4E %9.3f",
					i+1, m_geometry->segment_tags[i], m_geometry->x[i], m_geometry->y[i], m_geometry->z[i], m_geometry->segment_length[i],
					real(curi), imag(curi), cmag, ph );*/
			
				q_density_nb_elements++;
				_q_density_segment_number.push_back(i+1);
				_q_density_segment_tag.push_back(m_geometry->segment_tags[i]);
				_q_density_segment_center_x.push_back(m_geometry->x[i]);
				_q_density_segment_center_y.push_back(m_geometry->y[i]);
				_q_density_segment_center_z.push_back(m_geometry->z[i]);
				_q_density_segment_length.push_back(m_geometry->segment_length[i]);
				_q_density.push_back(curi);
						
			} /* for(int i = 0; i < m_geometry->n; i++ ) */
	
		} /* if (iptflq != -1) */	
	} /* if (m_geometry->n != 0) */
	
	if (m_geometry->m != 0)
	{
		int j = m_geometry->n-3;
		int itmp1 = -1;

		for(int i = 0; i < m_geometry->m; i++ )
		{
			j += 3;
			itmp1++;
			ASSERT(itmp1 == i);
			
			ex= m_context->current_vector[j];
			ey= m_context->current_vector[j+1];
			ez= m_context->current_vector[j+2];
			eth= ex* m_geometry->t1x[itmp1]+ ey* m_geometry->t1y[itmp1]+ ez* m_geometry->t1z[itmp1];
			eph= ex* m_geometry->t2x[itmp1]+ ey* m_geometry->t2y[itmp1]+ ez* m_geometry->t2z[itmp1];			

			/*m_output.nec_printf(
			      " %4d %7.3f %7.3f %7.3f %11.4E "
			      "%8.2f %11.4E %8.2f"
			      " %9.2E %9.2E %9.2E %9.2E %9.2E %9.2E",
			      i+1, m_geometry->px[itmp1], m_geometry->py[itmp1], m_geometry->pz[itmp1],
			      ethm, etha, ephm, epha, real(ex), imag(ex),
			      real(ey), imag(ey), real(ez), imag(ez));*/
			      
			patch_nb_elements++;
			_patch_number.push_back(i+1);
			_patch_center_x.push_back(m_geometry->px[itmp1]);
			_patch_center_y.push_back(m_geometry->py[itmp1]);
			_patch_center_z.push_back(m_geometry->pz[itmp1]);
			_patch_tangent_vector1.push_back(eth);      
			_patch_tangent_vector2.push_back(eph);
			_patch_e_x.push_back(ex);
			_patch_e_y.push_back(ey);
			_patch_e_z.push_back(ez);
			
		  	m_context->plot_card.plot_currents(ex,ey,ez);
		} /* for( i=0; i<m; i++ ) */
			
	}/* if (m_geometry->m != 0) */ 
}

void nec_structure_currents::write_to_file_aux(ostream& os)
{
	
	output_helper oh(os,_result_format);

	if ( m_geometry->n != 0)
	{
		if ( iptflg != -1)
		{
			if ( iptflg <= 0)
			{
				oh.section_start("CURRENTS AND LOCATION");
				oh.center_text("DISTANCES IN WAVELENGTHS","");
				os << endl;
				os << "   SEG  TAG    COORDINATES OF SEGM CENTER     SEGM    ------------- CURRENT (AMPS) -------------" << endl;
				os << "   No:  No:       X         Y         Z      LENGTH     REAL      IMAGINARY    MAGN        PHASE" << endl;
				
				for(int i=0; i<current_nb_elements; i++)
				{
					oh.start_record();
					oh.padding(" ");
					oh.int_out(5, _current_segment_number[i]); oh.separator();
					oh.int_out(4, _current_segment_tag[i]);	oh.separator();
					oh.real_out(9, 4, _current_segment_center_x[i], false); oh.separator(); 
					oh.real_out(9, 4, _current_segment_center_y[i], false); oh.separator(); 
					oh.real_out(9, 4, _current_segment_center_z[i], false); oh.separator(); 
					oh.real_out(9, 5, _current_segment_length[i], false); oh.separator(); 
					oh.real_out(11,4, real(_current[i]), true); oh.separator();
					oh.real_out(11,4, imag(_current[i]), true); oh.separator();
					oh.real_out(11,4, abs(_current[i]), true); oh.separator();
					oh.real_out(8,3, arg_degrees(_current[i]), false); oh.separator();				
					oh.end_record(); 
				}
			}
			else if (iptflg != 3)
			{
				if (m_context->get_inc() <= 1)
				{
					oh.section_start("RECEIVING PATTERN PARAMETERS");
					os << "                      ETA: "; oh.real_out(7,2,xpr3,false); os << " DEGREES" << endl;
					os << "                      TYPE: "; oh.string_out(6, nec_structure_currents::hpol(pattype)); os << endl;
					os << "                      AXIAL RATIO: "; oh.real_out(6,3,xpr6,false); os << endl << endl;
					os << "            THETA     PHI      ----- CURRENT ----    SEG" << endl;
					os << "            (DEG)    (DEG)     MAGNITUDE    PHASE    No:" << endl;
				
				}
				
				int i = current_nb_elements-1;
				oh.start_record();
				oh.padding("          ");
				oh.real_out(7, 2, _current_theta[i], false); oh.separator(); 
				oh.real_out(7, 2, _current_phi[i], false); oh.separator();
				oh.padding("   "); 
				oh.real_out(11,4, abs(_current[i]), true); oh.separator();
				oh.padding(" ");
				oh.real_out(7,2, arg_degrees(_current[i]), false); oh.separator();				
				oh.padding(" ");
				oh.int_out(5, _current_segment_number[i]); oh.separator();
				oh.end_record();
				
				
			} /* if ( iptflg <= 0) */
		} /* if ( iptflg != -1) */
		
		if ( iptflq != -1)
		{
			oh.section_start("CHARGE DENSITIES");
			oh.center_text("DISTANCES IN WAVELENGTHS","");
			os << endl;
			os << "   SEG   TAG    COORDINATES OF SEG CENTER     SEG          CHARGE DENSITY (COULOMBS/METER)" << endl;
			os << "   NO:   NO:     X         Y         Z       LENGTH     REAL      IMAGINARY     MAGN        PHASE" << endl;
			
			for(int i=q_density_last_printed; i< q_density_nb_elements; i++)
			{
				oh.start_record();
				oh.padding(" ");
				oh.int_out(5, _q_density_segment_number[i]); oh.separator();
				oh.int_out(4, _q_density_segment_tag[i]);	oh.separator();
				oh.real_out(9, 4,_q_density_segment_center_x[i], false); oh.separator(); 
				oh.real_out(9, 4, _q_density_segment_center_y[i], false); oh.separator(); 
				oh.real_out(9, 4, _q_density_segment_center_z[i], false); oh.separator(); 
				oh.real_out(9, 5, _q_density_segment_length[i], false); oh.separator(); 
				oh.real_out(11,4, real(_q_density[i]), true); oh.separator();
				oh.real_out(11,4, imag(_q_density[i]), true); oh.separator();
				oh.real_out(11,4, abs(_q_density[i]), true); oh.separator();
				oh.real_out(9,3, arg_degrees(_q_density[i]), false); oh.separator();				
				oh.end_record(); 
			}
			q_density_last_printed = q_density_nb_elements;
		}
	} /*if ( m_geometry->n != 0) */
	
	if ( m_geometry->m != 0)
	{
		oh.section_start("SURFACE PATCH CURRENTS");
		oh.center_text("DISTANCES IN WAVELENGTHS","");
		oh.center_text("CURRENT IN AMPS/METER","");
		os << endl << endl;
		oh.center_text("SURFACE COMPONENTS","---------");
		oh.center_text("RECTANGULAR COMPONENTS","---------");
		os << "  PCH   --- PATCH CENTER ---     TANGENT VECTOR 1     TANGENT VECTOR 2";
		os << "    ------- X ------    ------- Y ------	  ------- Z ------" << endl;
		os << "  No:    X       Y       Z       MAG.       PHASE     MAG.       PHASE";
		os << "    REAL   IMAGINARY    REAL   IMAGINARY    REAL   IMAGINARY" << endl;
		
		for(int i=0; i<patch_nb_elements; i++)
		{
			oh.start_record();
			oh.padding(" ");
			oh.int_out(4, _patch_number[i]); oh.separator();
			oh.real_out(7, 3,_patch_center_x[i], false); oh.separator(); 
			oh.real_out(7, 3, _patch_center_y[i], false); oh.separator(); 
			oh.real_out(7, 3, _patch_center_z[i], false); oh.separator(); 
			oh.real_out(11,4, abs(_patch_tangent_vector1[i]), true); oh.separator(); 
			oh.real_out(8, 2, arg_degrees(_patch_tangent_vector1[i]), false); oh.separator();
			oh.real_out(11,4, abs(_patch_tangent_vector2[i]), true); oh.separator(); 
			oh.real_out(8, 2, arg_degrees(_patch_tangent_vector2[i]), false); oh.separator();
			oh.real_out(9,2, real(_patch_e_x[i]), true); oh.separator();
			oh.real_out(9,2, imag(_patch_e_x[i]), true); oh.separator();
			oh.real_out(9,2, real(_patch_e_y[i]), true); oh.separator();
			oh.real_out(9,2, imag(_patch_e_y[i]), true); oh.separator();
			oh.real_out(9,2, real(_patch_e_z[i]), true); oh.separator();
			oh.real_out(9,2, imag(_patch_e_z[i]), true); oh.separator();
			oh.end_record(); 
		} /* for(int i=0; i<patch_nb_elements; i++) */
		
	} /* if ( m_geometry->m != 0) */
			    
} /* write_to_file_aux */

