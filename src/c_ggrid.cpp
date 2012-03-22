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
#include "c_ggrid.h"
#include "common.h"
#include "nec_output.h" // for DEBUG_TRACE()

#include <cstdlib>

//#define	CONST4	nec_complex(0.0,em::impedance() / 2.0)

int	c_ggrid::m_nxa[3] = {11,17,9};
int	c_ggrid::m_nya[3] = {10,5,8};

nec_float c_ggrid::m_dxa[3] = {.02,.05,.1}; 
nec_float c_ggrid::m_dya[3] = {.1745329252,.0872664626,.1745329252};

nec_float c_ggrid::m_xsa[3] = {0.,.2,.2};
nec_float c_ggrid::m_ysa[3] = {0.,0.,.3490658504};

/*! \brief interpolate (was intrp) uses bivariate cubic interpolation to obtain the values of 4 functions at the point (x,y).
*/
void c_ggrid::interpolate( nec_float x, nec_float y, nec_complex *f1,
    nec_complex *f2, nec_complex *f3, nec_complex *f4 )
{
	static int ix, iy, ixs=-10, iys=-10, igrs=-10, ixeg=0, iyeg=0;
	static int nxm2, nym2, nxms, nyms, nd, ndp;
	static nec_float dx = 1., dy = 1., xs = 0., ys = 0., xz, yz;
	static nec_complex a[4][4], b[4][4], c[4][4], d[4][4];
	static int nda[3] = {11,17,9}, ndpa[3] = {110, 85, 72};
	
	nec_complex p1, p2, p3, p4, fx1, fx2, fx3, fx4;
	
	bool skip_recalculation = false;
	if( (x < xs) || (y < ys) )
		skip_recalculation = true;
	else
	{
		ix = (int)((x-xs) / dx)+1;
		iy = (int)((y-ys) / dy)+1;
	}
	
	/* if point lies in same 4 by 4 point region */
	/* as previous point, old values are reused. */
	if( (ix < ixeg) ||
		(iy < iyeg) ||
		(std::abs(ix - ixs) >= 2) ||
		(std::abs(iy - iys) >= 2) ||
		(skip_recalculation == false) )
	{
		/* determine correct grid and grid region */
		int igr;
		
		if( x <= m_xsa[1])
			igr=0;
		else
		{
			if( y > m_ysa[2])
				igr=2;
			else
				igr=1;
		}
	
		if( igr != igrs)
		{
			igrs= igr;
			dx= m_dxa[igrs];
			dy= m_dya[igrs];
			xs= m_xsa[igrs];
			ys= m_ysa[igrs];
			nxm2= m_nxa[igrs]-2;
			nym2= m_nya[igrs]-2;
			nxms=(( nxm2+1)/3)*3+1;
			nyms=(( nym2+1)/3)*3+1;
			nd= nda[igrs];
			ndp= ndpa[igrs];
			ix= (int)(( x- xs)/ dx)+1;
			iy= (int)(( y- ys)/ dy)+1;
		} /* if( igr != igrs) */
	
		ixs=(( ix-1)/3)*3+2;
		if( ixs < 2)
			ixs=2;
		ixeg=-10000;
	
		if( ixs > nxm2)
		{
			ixs= nxm2;
			ixeg= nxms;
		}
	
		iys=(( iy-1)/3)*3+2;
		if( iys < 2)
			iys=2;
		iyeg=-10000;
	
		if( iys > nym2)
		{
			iys= nym2;
			iyeg= nyms;
		}
	
		/* compute coefficients of 4 cubic polynomials in x for */
		/* the 4 grid values of y for each of the 4 functions */
		int iadz= ixs+( iys-3)* nd- ndp;
		for (int k = 0; k < 4; k++ )
		{
			iadz += ndp;
			int iadd = iadz;
		
			for (int i = 0; i < 4; i++ )
			{
				iadd += nd;
			
				switch( igrs )
				{
				case 0:
					p1= m_ar1[iadd-2];
					p2= m_ar1[iadd-1];
					p3= m_ar1[iadd];
					p4= m_ar1[iadd+1];
					break;
			
				case 1:
					p1= m_ar2[iadd-2];
					p2= m_ar2[iadd-1];
					p3= m_ar2[iadd];
					p4= m_ar2[iadd+1];
					break;
			
				case 2:
					p1= m_ar3[iadd-2];
					p2= m_ar3[iadd-1];
					p3= m_ar3[iadd];
					p4= m_ar3[iadd+1];
				} /* switch( igrs ) */
			
				a[i][k]=( p4- p1+3.*( p2- p3))*.1666666667;
				b[i][k]=( p1-2.* p2+ p3)*.5;
				c[i][k]= p3-(2.* p1+3.* p2+ p4)*.1666666667;
				d[i][k]= p2;
			
			} /* for ( i = 0; i < 4; i++ ) */
		
		} /* for ( k = 0; k < 4; k++ ) */
	
		xz=( ixs-1)* dx+ xs;
		yz=( iys-1)* dy+ ys;
	
	} /* if( (abs(ix- ixs) >= 2) || */
	
	/* evaluate polymomials in x and use cubic */
	/* interpolation in y for each of the 4 functions. */
	nec_float xx=( x- xz)/ dx;
	nec_float yy=( y- yz)/ dy;
	fx1=(( a[0][0]* xx+ b[0][0])* xx+ c[0][0])* xx+ d[0][0];
	fx2=(( a[1][0]* xx+ b[1][0])* xx+ c[1][0])* xx+ d[1][0];
	fx3=(( a[2][0]* xx+ b[2][0])* xx+ c[2][0])* xx+ d[2][0];
	fx4=(( a[3][0]* xx+ b[3][0])* xx+ c[3][0])* xx+ d[3][0];
	p1= fx4- fx1+3.*( fx2- fx3);
	p2=3.*( fx1-2.* fx2+ fx3);
	p3=6.* fx3-2.* fx1-3.* fx2- fx4;
	*f1=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
	fx1=(( a[0][1]* xx+ b[0][1])* xx+ c[0][1])* xx+ d[0][1];
	fx2=(( a[1][1]* xx+ b[1][1])* xx+ c[1][1])* xx+ d[1][1];
	fx3=(( a[2][1]* xx+ b[2][1])* xx+ c[2][1])* xx+ d[2][1];
	fx4=(( a[3][1]* xx+ b[3][1])* xx+ c[3][1])* xx+ d[3][1];
	p1= fx4- fx1+3.*( fx2- fx3);
	p2=3.*( fx1-2.* fx2+ fx3);
	p3=6.* fx3-2.* fx1-3.* fx2- fx4;
	*f2=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
	fx1=(( a[0][2]* xx+ b[0][2])* xx+ c[0][2])* xx+ d[0][2];
	fx2=(( a[1][2]* xx+ b[1][2])* xx+ c[1][2])* xx+ d[1][2];
	fx3=(( a[2][2]* xx+ b[2][2])* xx+ c[2][2])* xx+ d[2][2];
	fx4=(( a[3][2]* xx+ b[3][2])* xx+ c[3][2])* xx+ d[3][2];
	p1= fx4- fx1+3.*( fx2- fx3);
	p2=3.*( fx1-2.* fx2+ fx3);
	p3=6.* fx3-2.* fx1-3.* fx2- fx4;
	*f3=(( p1* yy+ p2)* yy+ p3)* yy*.1666666667+ fx2;
	fx1=(( a[0][3]* xx+ b[0][3])* xx+ c[0][3])* xx+ d[0][3];
	fx2=(( a[1][3]* xx+ b[1][3])* xx+ c[1][3])* xx+ d[1][3];
	fx3=(( a[2][3]* xx+ b[2][3])* xx+ c[2][3])* xx+ d[2][3];
	fx4=(( a[3][3]* xx+ b[3][3])* xx+ c[3][3])* xx+ d[3][3];
	p1= fx4- fx1+3.*( fx2- fx3);
	p2=3.*( fx1-2.* fx2+ fx3);
	p3=6.* fx3-2.* fx1-3.* fx2- fx4;
	*f4=(( p1* yy+ p2)* yy+ p3)* yy*.16666666670+ fx2;
}

#include "electromag.h"

/*! was SOMNEC in the original NEC-2 source code
*/
void c_ggrid::sommerfeld( nec_float epr, nec_float sig, nec_float freq_mhz )
{
	static nec_complex const1_neg = - nec_complex(0.0,4.771341189);
	static nec_complex CONST4(0.0, em::impedance() / 2.0);
	
	nec_float wavelength, dr, dth, r, rk, thet, tfac1, tfac2;
	nec_complex erv, ezv, erh, eph, cl1, cl2, con;
	
	if(sig >= 0.0)
	{
		wavelength = em::speed_of_light() / (1.0e6 * freq_mhz);
		m_epscf = nec_complex(epr,-sig*wavelength*em::impedance_over_2pi());
	}
	else
		m_epscf=nec_complex(epr,sig);
	
	m_evlcom.m_ck2 = two_pi();
	m_evlcom.m_ck2sq = m_evlcom.m_ck2*m_evlcom.m_ck2;
	
	/*
	Sommerfeld integral evaluation uses exp(-jwt), NEC uses exp(+jwt),
	hence need conjg(epscf).  Conjugate of fields occurs in subroutine 
	evlua. */
	
	m_evlcom.m_ck1sq=m_evlcom.m_ck2sq*conj(m_epscf);
	m_evlcom.m_ck1=sqrt(m_evlcom.m_ck1sq);
	m_evlcom.m_ck1r=real(m_evlcom.m_ck1);
	m_evlcom.m_tkmag=100.0*abs(m_evlcom.m_ck1);
	m_evlcom.m_tsmag=100.0*norm(m_evlcom.m_ck1); // TCAM changed from previous line
	m_evlcom.m_cksm=m_evlcom.m_ck2sq/(m_evlcom.m_ck1sq+m_evlcom.m_ck2sq);
	m_evlcom.m_ct1=.5*(m_evlcom.m_ck1sq-m_evlcom.m_ck2sq);
	erv=m_evlcom.m_ck1sq*m_evlcom.m_ck1sq;
	ezv=m_evlcom.m_ck2sq*m_evlcom.m_ck2sq;
	m_evlcom.m_ct2=.125*(erv-ezv);
	erv *= m_evlcom.m_ck1sq;
	ezv *= m_evlcom.m_ck2sq;
	m_evlcom.m_ct3=.0625*(erv-ezv);
	
	/* loop over 3 grid regions */
	for (int k = 0; k < 3; k++ )
	{
		int nr = m_nxa[k];
		int nth = m_nya[k];
		dr = m_dxa[k];
		dth = m_dya[k];
		r = m_xsa[k]-dr;
		int irs=1;
		
		if(k == 0)
		{
			r=m_xsa[k];
			irs=2;
		}
	
		/*  loop over r.  (r=sqrt(m_evlcom.m_rho**2 + (z+h)**2)) */
		for (int ir = irs-1; ir < nr; ir++ )
		{
			r += dr;
			thet = m_ysa[k]-dth;
		
			/* loop over theta.  (theta=atan((z+h)/m_evlcom.m_rho)) */
			for (int ith = 0; ith < nth; ith++ )
			{
				thet += dth;
				m_evlcom.m_rho=r*cos(thet);
				m_evlcom.m_zph=r*sin(thet);
				if(m_evlcom.m_rho < 1.e-7)
					m_evlcom.m_rho=1.e-8;
				if(m_evlcom.m_zph < 1.e-7)
					m_evlcom.m_zph=0.;
			
				m_evlcom.evlua( &erv, &ezv, &erh, &eph );
			
				rk=m_evlcom.m_ck2*r;
				con=const1_neg*r/nec_complex(cos(rk),-sin(rk));
			
				switch( k )
				{
				case 0:
					m_ar1[ir+ith*11+  0]=erv*con;
					m_ar1[ir+ith*11+110]=ezv*con;
					m_ar1[ir+ith*11+220]=erh*con;
					m_ar1[ir+ith*11+330]=eph*con;
					break;
			
				case 1:
					m_ar2[ir+ith*17+  0]=erv*con;
					m_ar2[ir+ith*17+ 85]=ezv*con;
					m_ar2[ir+ith*17+170]=erh*con;
					m_ar2[ir+ith*17+255]=eph*con;
					break;
			
				case 2:
					m_ar3[ir+ith*9+  0]=erv*con;
					m_ar3[ir+ith*9+ 72]=ezv*con;
					m_ar3[ir+ith*9+144]=erh*con;
					m_ar3[ir+ith*9+216]=eph*con;
			
				} /* switch( k ) */
			
			} /* for ( ith = 0; ith < nth; ith++ ) */
		} /* for ( ir = irs-1; ir < nr; ir++; ) */
	} /* for ( k = 0; k < 3; k++; ) */
	
	/* fill grid 1 for r equal to zero. */
	cl2 = -CONST4*(m_epscf-1.)/(m_epscf+1.);
	cl1 = cl2/(m_epscf+1.);
	ezv = m_epscf*cl1;
	thet=-dth;
	int nth = m_nya[0];
	
	for (int ith = 0; ith < nth; ith++ )
	{
		thet += dth;
		if( (ith+1) != nth )
		{
			tfac2=cos(thet);
			tfac1=(1.-sin(thet))/tfac2;
			tfac2=tfac1/tfac2;
			erv=m_epscf*cl1*tfac1;
			erh=cl1*(tfac2-1.)+cl2;
			eph=cl1*tfac2-cl2;
		}
		else
		{
			erv=0.;
			erh=cl2-.5*cl1;
			eph=-erh;
		}
	
		m_ar1[0+ith*11+  0]=erv;
		m_ar1[0+ith*11+110]=ezv;
		m_ar1[0+ith*11+220]=erh;
		m_ar1[0+ith*11+330]=eph;
	}
}



/* fbar is the Sommerfeld attenuation function for numerical distance . */
nec_complex  fbar(const nec_complex& p );
nec_complex  fbar(const nec_complex& p )
{
	static nec_float TOSP = 2.0 / sqrt_pi();

	int minus;
	nec_float tms, sms;
	nec_complex z, zs, sum, pow, term, fbar;
	
	z= cplx_01()* sqrt( p);
	if ( abs( z) <= 3.)
	{
		/* series expansion */
		zs= z* z;
		sum= z;
		pow= z;
	
		for (int i = 1; i <= 100; i++ )
		{
			pow=- pow* zs/ (nec_float)i;
			term= pow/(2.* i+1.);
			sum = sum + term;
			tms = norm(term);
			sms = norm(sum);
			
			if ( tms/sms < ACCS)
				break;
		}
	
		fbar=1.-(1.- sum* TOSP)* z* exp( zs)* sqrt_pi();
		return( fbar );
	
	} /* if ( abs( z) <= 3.) */
	
	/* asymptotic expansion */
	if ( real( z) < 0.)
	{
		minus=1;
		z=- z;
	}
	else
		minus=0;
	
	zs=.5/( z* z);
	sum=cplx_00();
	term=cplx_10();
	
	for (int i = 1; i <= 6; i++ )
	{
		term =- term*(2.*i -1.)* zs;
		sum += term;
	}
	
	if ( minus == 1)
		sum -= 2.0 * sqrt_pi() * z* exp( z* z);
	fbar=- sum;
	
	return( fbar );
}



/* gwave computes the electric field, including ground wave, of a */
/* current element over a ground plane using formulas of k.a. norton */
/* (proc. ire, sept., 1937, pp.1203,1236.) */

void gwave(nec_complex& erv, nec_complex& ezv,
	nec_complex& erh, nec_complex& ezh, nec_complex& eph,
	c_ground_wave& ground_wave)
{
	static nec_complex CONST4(0.0,em::impedance() / 2.0);

	nec_float sppp, sppp2, cppp2, cppp, spp, spp2, cpp2, cpp;
	nec_complex rk1, rk2, t1, t2, t3, t4, p1, rv;
	nec_complex omr, w, f, q1, rh, v, g, xr1, xr2;
	nec_complex x1, x2, x3, x4, x5, x6, x7;
	
	sppp= ground_wave.zmh/ ground_wave.r1;
	sppp2= sppp* sppp;
	cppp2=1.- sppp2;
	
	if ( cppp2 < 1.0e-20)
		cppp2=1.0e-20;
	
	cppp= sqrt( cppp2);
	spp= ground_wave.zph/ ground_wave.r2;
	spp2= spp* spp;
	cpp2=1.- spp2;
	
	if ( cpp2 < 1.0e-20)
		cpp2=1.0e-20;
	
	cpp= sqrt( cpp2);
	rk1=- two_pi_j()* ground_wave.r1;
	rk2=- two_pi_j()* ground_wave.r2;
	t1=1. -ground_wave.u2* cpp2;
	t2= sqrt( t1);
	t3=(1. -1./ rk1)/ rk1;
	t4=(1. -1./ rk2)/ rk2;
	p1= rk2* ground_wave.u2* t1/(2.* cpp2);
	rv=( spp- ground_wave.u* t2)/( spp+ ground_wave.u* t2);
	omr=1.- rv;
	w=1./ omr;
	w= nec_complex(4.0,0.0)* p1* w* w;
	f= fbar( w);
	q1= rk2* t1/(2.* ground_wave.u2* cpp2);
	rh=( t2- ground_wave.u* spp)/( t2+ ground_wave.u* spp);
	v=1./(1.+ rh);
	v=nec_complex(4.0,0.0)* q1* v* v;
	g= fbar( v);
	xr1= ground_wave.xx1/ ground_wave.r1;
	xr2= ground_wave.xx2/ ground_wave.r2;
	x1= cppp2* xr1;
	x2= rv* cpp2* xr2;
	x3= omr* cpp2* f* xr2;
	x4= ground_wave.u* t2* spp*2.* xr2/ rk2;
	x5= xr1* t3*(1.-3.* sppp2);
	x6= xr2* t4*(1.-3.* spp2);
	ezv=( x1+ x2+ x3- x4- x5- x6)* (-CONST4);

	x1= sppp* cppp* xr1;
	x2= rv* spp* cpp* xr2;
	x3= cpp* omr* ground_wave.u* t2* f* xr2;
	x4= spp* cpp* omr* xr2/ rk2;
	x5=3.* sppp* cppp* t3* xr1;
	x6= cpp* ground_wave.u* t2* omr* xr2/ rk2*.5;
	x7=3.* spp* cpp* t4* xr2;
	erv=-( x1+ x2- x3+ x4- x5+ x6- x7)* (-CONST4);

	ezh=-( x1- x2+ x3- x4- x5- x6+ x7)* (-CONST4);

	x1= sppp2* xr1;
	x2= rv* spp2* xr2;
	x4= ground_wave.u2* t1* omr* f* xr2;
	x5= t3*(1.-3.* cppp2)* xr1;
	x6= t4*(1.-3.* cpp2)*(1.- ground_wave.u2*(1.+ rv)- ground_wave.u2* omr* f)* xr2;
	x7= ground_wave.u2* cpp2* omr*(1.-1./ rk2)*( f*( ground_wave.u2* t1- spp2-1./ rk2)+1./rk2)* xr2;
	erh=( x1- x2- x4- x5+ x6+ x7)* (-CONST4);

	x1= xr1;
	x2= rh* xr2;
	x3=( rh+1.)* g* xr2;
	x4= t3* xr1;
	x5= t4*(1.- ground_wave.u2*(1.+ rv)- ground_wave.u2* omr* f)* xr2;
	x6=.5* ground_wave.u2* omr*( f*( ground_wave.u2* t1- spp2-1./ rk2)+1./ rk2)* xr2/ rk2;
	eph=-( x1- x2+ x3- x4+ x5+ x6)* (-CONST4);
}



