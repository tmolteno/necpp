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
#ifndef __c_plot_card__
#define __c_plot_card__


#include "math_util.h"
#include <string>
#include <stdio.h>

/*!
	All the logic for handling the plot card is contained in this class.
	There appear to be many variants of the PL card. I am using the documentation
	below as a guide.
	
\verbatim
PL		PLOT DATA STORAGE
	I1- 0=NO STORE, 1=CURRENTS, 2=NEAR FIELD, 3=PATTERNS
	I2(1,2)- 0=NO, 1=REAL,IMAG, 3=MAG, PHASE
	I2(3)- 1=TH, 2=PHI, 3=RHO ANGLES
	I3(1)- 0=N0, 1=IX, 2=IY, 3=IZ, 4=IX IY IZ
	I3(2)- 0-N0, 1=X, 2=Y, 3=Z, 4=X Y Z, 5=TOTAL  COMPONENT
	I3(3)- 0=NO, 1=TH, 2=PHI, 3=RHO E-FIELD COMPOMENT
 	I4(1)- BLANK
	I4(2)- 1=X, 2=Y, 3=Z CORDINATE VALUES
	I4(3)- 1=V, 2=H, 3=TOTAL, 4=V H T     GAINS DB
\endverbatim
*/
class c_plot_card
{
public:
	c_plot_card();
	c_plot_card(const c_plot_card& p);
	c_plot_card(int itmp1, int itmp2, int itmp3, int itmp4, std::string& filename);
	
	virtual ~c_plot_card();
	
	c_plot_card& operator=(const c_plot_card&) = default;
	
	bool is_valid()	const;
	
	bool storing()	const;
	bool currents()	const;
	bool near_field() const	;
	bool patterns()	const;
	
	bool realimag()	const;
	bool magphase()	const;

	void set_plot_real_imag_currents();
	
	void plot_endl() const;
	
	void plot_double(nec_float x) const;

	void plot_complex(nec_complex x) const;
	
	void plot_complex_2d(nec_complex x, nec_complex y, nec_complex z) const;

	void plot_currents(nec_complex ex, nec_complex ey, nec_complex ez) const;
	
	void plot_segments(int i,
		real_array& x, real_array& y, real_array& z, real_array& si,
		nec_float xw2, nec_float yw2,
		real_array& bi, int_array& icon1, int_array& icon2) const;
		
	void plot_fields(
		nec_complex ex, nec_complex ey, nec_complex ez, 
		nec_float xob, nec_float yob, nec_float zob);

	void plot_patterns(nec_float theta, nec_float phi,
		nec_complex e_theta, nec_complex e_phi,
		nec_float g_vert, nec_float g_horiz, nec_float g_tot);
private:
	int p1, p2, p3, p4;
	FILE* plot_fp;
};

#endif /* __c_plot_card__ */
