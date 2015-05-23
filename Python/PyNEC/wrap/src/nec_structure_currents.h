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
#ifndef __nec_structure_currents__
#define __nec_structure_currents__

#include "nec_results.h"
#include "math_util.h"
#include "nec_context.h"


class nec_context;
class c_geometry;

class nec_structure_currents : public nec_base_result
{
public:
	/*Structure currents*/
	nec_structure_currents(nec_context * in_context, enum excitation_type in_pattype,
			int in_nload,
			nec_float in_xpr3, nec_float in_xpr6);
	
	static std::string hpol(enum excitation_type e);

	void analyze();
	
	virtual ~nec_structure_currents()
	{
	}
		 
	virtual void write_to_file(ostream& os)
	{
		write_to_file_aux(os);
	}
	
	virtual enum nec_result_type get_result_type()
	{
		return RESULT_STRUCTURE_CURRENTS;
	} 
	
	int get_iptflg()
	{
		return iptflg;
	}
	
	int get_iptflq()
	{
		return iptflq;
	}
	
	int get_n();
	
	int get_m();
			
	vector<int> get_current_segment_number()
	{
		return _current_segment_number;
	}

	vector<int> get_current_segment_tag()
	{
		return _current_segment_tag;
	}
	
	vector<nec_float> get_current_segment_center_x()
	{
		return _current_segment_center_x;
	}
	
	vector<nec_float> get_current_segment_center_y()
	{
		return _current_segment_center_y;
	}
	
	vector<nec_float> get_current_segment_center_z()
	{
		return _current_segment_center_z;
	}
	
	vector<nec_float> get_current_segment_length()
	{
		return _current_segment_length;
	}
	
	vector<nec_float> get_current_theta()
	{
		return _current_theta;
	}
		
	vector<nec_float> get_current_phi()
	{
		return _current_phi;
	}
	
	vector<nec_complex> get_current()
	{
		return _current;
	}
	
	vector<int> get_q_density_segment_number()
	{
		return _q_density_segment_number;
	}
	
	vector<int> get_q_density_segment_tag()
	{
		return _q_density_segment_tag;
	}
	
	vector<nec_float> get_q_density_segment_center_x()
	{
		return _q_density_segment_center_x;
	}
	
	vector<nec_float> get_q_density_segment_center_y()
	{
		return _q_density_segment_center_y;
	}
	
	vector<nec_float> get_q_density_segment_center_z()
	{
		return _q_density_segment_center_z;
	}
	
	vector<nec_float> get_q_density_segment_length()
	{
		return _q_density_segment_length;
	}
	
	vector<nec_complex> get_q_density()
	{
		return _q_density;
	}
	
	vector<int> get_patch_number()
	{
		return _patch_number;
	}
	
	vector<nec_float> get_patch_center_x()
	{
		return _patch_center_x;
	}
	
	vector<nec_float> get_patch_center_y()
	{
		return _patch_center_y;
	}
	
	vector<nec_float> get_patch_center_z()
	{
		return _patch_center_z;
	}
	
	vector<nec_complex> get_patch_tangent_vector1()
	{
		return _patch_tangent_vector1;
	}
	
	vector<nec_complex> get_patch_tangent_vector2()
	{
		return _patch_tangent_vector2;
	}
	
	vector<nec_complex> get_patch_e_x()
	{
		return _patch_e_x;
	}
	
	vector<nec_complex> get_patch_e_y()
	{
		return _patch_e_y;
	}
	
	vector<nec_complex> get_patch_e_z()
	{
		return _patch_e_z;
	}

private:
	
	nec_context *m_context;
	c_geometry * m_geometry;
	enum excitation_type pattype;
	
	int iptflg;
	int iptag, iptagf, iptagt;
	int iptflq;
	int iptaq, iptaqf, iptaqt;
	
	int nload;
	
	nec_float xpr3, xpr6;
	nec_float wavelength;
	nec_float freq_mhz;
	
	nec_float structure_power_loss;
		
	int current_nb_elements;
	int q_density_nb_elements;
	int q_density_last_printed;
	int patch_nb_elements;
	
	complex_array zarray;
	real_array fnorm;	
		
	vector<int> _current_segment_number;
	vector<int> _current_segment_tag;
	
	vector<nec_float> _current_segment_center_x, _current_segment_center_y, _current_segment_center_z;
	vector<nec_float> _current_segment_length;
	vector<nec_float> _current_theta, _current_phi;
	
	vector<nec_complex> _current;
	
	vector<int> _q_density_segment_number;
	vector<int> _q_density_segment_tag;
	
	vector<nec_float> _q_density_segment_center_x, _q_density_segment_center_y, _q_density_segment_center_z;
	vector<nec_float> _q_density_segment_length;
	
	vector<nec_complex> _q_density;
	
	vector<int> _patch_number;
	
	vector<nec_float> _patch_center_x, _patch_center_y, _patch_center_z;
	
	vector<nec_complex> _patch_tangent_vector1, _patch_tangent_vector2;
	vector<nec_complex> _patch_e_x, _patch_e_y, _patch_e_z;
	
	void write_to_file_aux(ostream& os);
};

#endif
