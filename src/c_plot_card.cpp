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
#include "c_plot_card.h"
#include "math_util.h"


#include "math_util.h"
#include <string>
#include <stdio.h>

using namespace std;

c_plot_card::c_plot_card()
{
	p1 = 0;
	p2 = 0;
	p3 = 0;
	p4 = 0;
	
	plot_fp = NULL;
}
	
c_plot_card::c_plot_card(const c_plot_card& p)
{
	p1 = p.p1;
	p2 = p.p2;
	p3 = p.p3;
	p4 = p.p4;
	
	plot_fp = p.plot_fp;
}
	
c_plot_card::c_plot_card(int itmp1, int itmp2, int itmp3, int itmp4, string& filename)
{
	p1 = itmp1;
	p2 = itmp2;
	p3 = itmp3;
	p4 = itmp4;
	
	plot_fp = NULL;
	
	/* Open plot file */
	if ( (plot_fp = fopen(filename.c_str(), "w")) == NULL )
	{
		throw 100;
	}
}

c_plot_card::~c_plot_card()
{
	if ( plot_fp != NULL )
		fclose( plot_fp );
}

bool c_plot_card::is_valid()	const	{	return (NULL != plot_fp); 	}

bool c_plot_card::storing()	const	{	return (p1 != 0); 	}
bool c_plot_card::currents()	const	{	return (p1 == 1); 	}
bool c_plot_card::near_field() const	{	return (p1 == 2); 	}
bool c_plot_card::patterns()	const	{	return (p1 == 3); 	}

bool c_plot_card::realimag()	const	{ 	return ( ((p1 == 1) || (p1 == 2)) && (p2 == 1) );	}
bool c_plot_card::magphase()	const	{ 	return ( ((p1 == 1) || (p1 == 2)) && (p2 == 3) );	}

void c_plot_card::set_plot_real_imag_currents()
{
	p1 = 1;
	p2 = 1;
}


void c_plot_card::plot_endl() const
{
	if (NULL == plot_fp)
		throw 100;
	
	fprintf( plot_fp, "\n");
}

void c_plot_card::plot_double(nec_float x) const
{
	if (NULL == plot_fp)
		throw 100;
	
	fprintf( plot_fp, "%12.4E ", x );
}

void c_plot_card::plot_complex(nec_complex x) const
{
	if (NULL == plot_fp)
		throw 100;
	
	switch( p2 )
	{
		case 2:
			plot_double(real(x));
			plot_double(imag(x));
		case 3:
			plot_double(abs(x));
			plot_double(arg_degrees(x));
	}
}

void c_plot_card::plot_complex_2d(nec_complex x, nec_complex y, nec_complex z) const
{
	switch( p3 )
	{
		case 1:
			plot_complex(x);
			break;
		case 2:
			plot_complex(y);
			break;
		case 3:
			plot_complex(z);
			break;
		case 4:
			plot_complex(x);
			plot_complex(y);
			plot_complex(z);
			break;
	}
}

void c_plot_card::plot_currents(nec_complex ex, nec_complex ey, nec_complex ez) const
{
	if (false == currents())
		return;
		
	plot_complex_2d(ex, ey, ez);
	plot_endl();
	
}

void c_plot_card::plot_segments(int i,
	real_array& x, real_array& y, real_array& z, real_array& si,
	nec_float xw2, nec_float yw2,
	real_array& bi, int_array& icon1, int_array& icon2) const
{
	if (false == near_field())
		return;
		
	fprintf( plot_fp, "%12.4E %12.4E %12.4E "
		"%12.4E %12.4E %12.4E %12.4E %5d %5d %5d\n",
		x[i],y[i],z[i],si[i],xw2,yw2,bi[i],icon1[i],i+1,icon2[i] );
}

	
void c_plot_card::plot_fields(
	nec_complex ex, nec_complex ey, nec_complex ez, 
	nec_float xob, nec_float yob, nec_float zob)
{
	if ( p1 != 2)
		return;

	nec_float xxx;
	
	if ( p4 < 0 )
		xxx = xob;
	else if ( p4 == 0 )
		xxx= yob;
	else
		xxx= zob;

	plot_double(xxx);
	plot_complex_2d(ex,ey,ez);
	plot_endl();		
}

void c_plot_card::plot_patterns(nec_float theta, nec_float phi,
	nec_complex e_theta, nec_complex e_phi,
	nec_float g_vert, nec_float g_horiz, nec_float g_tot)
{
	if ( false == patterns())
		return;
	
	switch (p2)
	{
		case 1:
			plot_double(theta);
			plot_complex(e_theta);
			plot_endl();
			break;
		case 2:
			plot_double(phi);
			plot_complex(e_phi);
			plot_endl();
			break;
	}
	
	
	if ( p4 == 0 )
		return;

	// plot gains I4(3)- 1=V, 2=H, 3=TOTAL, 4=V H T     GAINS DB
	
	switch (p2)
	{
		case 1:
			plot_double(theta);
			break;
		case 2:
			plot_double(phi);
			break;
	}
	switch( p4 )
	{
		case 1:
			plot_double(g_vert); // vertical
			break;
		case 2:
			plot_double(g_horiz); // horizontal
			break;
		case 3:
			plot_double(g_tot); // total
			break;
		case 4:
			plot_double(g_vert);
			plot_double(g_horiz);
			plot_double(g_tot);
			break;
	}
	plot_endl();
}
