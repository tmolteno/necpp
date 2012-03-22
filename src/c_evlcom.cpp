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
#include "c_evlcom.h"

#include "matrix_algebra.h" // for test()
#include "nec_exception.h"

using namespace std;

/* 
	What on earth is this NM thing?
*/
#define NM	131072
#define MAXH	20
#define CRIT	1.0E-4
#define	PTP	.6283185308
#define NTS	4

/* compute integration parameter xlam=lambda from parameter t. */
void c_evlcom::lambda( nec_float t, nec_complex *xlam, nec_complex *dxlam ) const
{
	*dxlam = m_contour_b - m_contour_a;
	*xlam = m_contour_a + *dxlam*t;
}


/*! \brief gshank integrates the 6 Sommerfeld integrals from start to 
	infinity (until convergence) in lambda.  At the break point, bk, 
	the step increment may be changed from dela to delb.  Shank's 
	algorithm to accelerate convergence of a slowly converging series
	is used. */
void c_evlcom::gshank( nec_complex start, nec_complex dela, complex_array& sum,
	int nans, complex_array& seed, int ibk, nec_complex bk, nec_complex delb )
{
	bool brk = false;
	int ibx, jm;
	static nec_float rbk, amg, den, denm;
	nec_complex a1, a2, as1, as2, del, aa;
	nec_complex q1[6][20], q2[6][20];
	 
	complex_array ans1(6), ans2(6);
	
	rbk=real(bk);
	del=dela;
	if (ibk == 0)
		ibx=1;
	else
		ibx=0;
	
	// I believe that this is a spurious generalization for this routine. Hence the
	// following assert.
	ASSERT(nans == 6);
	
	for (int i = 0; i < nans; i++ )
		ans2[i]=seed[i];
	
	m_contour_b=start;
	for (int intx = 1; intx <= MAXH; intx++ )
	{
		int inx=intx-1;
		m_contour_a = m_contour_b;
		m_contour_b += del;
		
		if ( (ibx == 0) && (real(m_contour_b) >= rbk) )
		{
			/* hit break point.  reset seed and start over. */
			ibx=1;
			m_contour_b=bk;
			del=delb;
			rom1(nans,sum,2);
			if ( ibx != 2 )
			{
				for (int i = 0; i < nans; i++ )
					ans2[i] += sum[i];
				intx = 0;
				continue;
			}
			
			for (int i = 0; i < nans; i++ )
				ans2[i]=ans1[i]+sum[i];
			intx = 0;
			continue;
		} /* if ( (ibx == 0) && (real(m_contour_b) >= rbk) ) */
		
		rom1(nans,sum,2);
		for (int i = 0; i < nans; i++ )
			ans1[i] = ans2[i]+sum[i];
			
		m_contour_a = m_contour_b;
		m_contour_b += del;
		
		if ( (ibx == 0) && (real(m_contour_b) >= rbk) )
		{
			/* hit break point.  reset seed and start over. */
			ibx=2;
			m_contour_b=bk;
			del=delb;
			rom1(nans,sum,2);
			if ( ibx != 2 )
			{
				for (int i = 0; i < nans; i++ )
					ans2[i] += sum[i];
				intx = 0;
				continue;
			}
			
			for (int i = 0; i < nans; i++ )
				ans2[i] = ans1[i]+sum[i];
			intx = 0;
			continue;
		} /* if ( (ibx == 0) && (real(m_contour_b) >= rbk) ) */
		
		rom1(nans,sum,2);
		for (int i = 0; i < nans; i++ )
			ans2[i]=ans1[i]+sum[i];
		
		den=0.;
		for (int i = 0; i < nans; i++ )
		{
			as1=ans1[i];
			as2=ans2[i];
			
			if (intx >= 2)
			{
				for (int j = 1; j < intx; j++ )
				{
					jm=j-1;
					aa=q2[i][jm];
					a1=q1[i][jm]+as1-2.*aa;
				
					if ( (real(a1) != 0.) || (imag(a1) != 0.) )
					{
						a2=aa-q1[i][jm];
						a1=q1[i][jm]-a2*a2/a1;
					}
					else
						a1=q1[i][jm];
				
					a2=aa+as2-2.*as1;
					if ( (real(a2) != 0.) || (imag(a2) != 0.) )
						a2=aa-(as1-aa)*(as1-aa)/a2;
					else
						a2=aa;
				
					q1[i][jm]=as1;
					q2[i][jm]=as2;
					as1=a1;
					as2=a2;	
				}
			}
			
			q1[i][intx-1]=as1;
			q2[i][intx-1]=as2;
			amg=fabs(real(as2))+fabs(imag(as2));
			if (amg > den)
				den=amg;		
		} /* for ( i = 0; i < nans; i++ ) */
		
		denm=1.e-3*den*CRIT;
		jm=intx-3;
		if (jm < 1)
			jm = 1;
		
		for (int j = jm-1; j < intx; j++ )
		{
			brk = false;
			for (int i = 0; i < nans; i++ )
			{
				a1=q2[i][j];
				den=(fabs(real(a1))+fabs(imag(a1)))*CRIT;
				if (den < denm)
					den=denm;
				a1=q1[i][j]-a1;
				amg=fabs(real(a1)+fabs(imag(a1)));
				if (amg > den)
				{
					brk = true;
					break;
				}
			
			} /* for ( i = 0; i < nans; i++ ) */
			
			if ( brk )
				break;	
		} /* for ( j = jm-1; j < intx; j++ ) */
		
		if ( false == brk )
		{
			for (int i = 0; i < nans; i++ )
				sum[i]=.5*(q1[i][inx]+q2[i][inx]);
			return;
		}
		
	} /* for ( intx = 1; intx <= maxh; intx++ ) */
	
	/* No convergence */
	throw new nec_exception("No convergence in gshank() - aborting");
}


/*! \brief rom1 integrates the 6 Sommerfeld integrals from m_contour_a to m_contour_b in lambda.
	The method of variable interval width Romberg integration is used. */
void c_evlcom::rom1( int n, complex_array& sum, int nx )
{
	int ns, nt;
	static nec_float z, ze, s, ep, zend, dz=0.0, dzot=0.0, tr, ti;
	static nec_complex t00, t11, t02;
	
	static complex_array g1(6), g2(6), g3(6), g4(6), g5(6), t01(6), t10(6), t20(6);
	
	ASSERT(n == 6);
	
	z = 0.0;
	ze = 1.0;
	s = 1.0;
	ep = s / (1.0e4 * NM);
	zend=ze-ep;
	
	nec_complex cmplx_zero(0.0,0.0);
	for (int i = 0; i < n; i++ )
		sum[i]=cmplx_zero;
		
	ns=nx;
	nt=0;
	saoa(z,g1);
	
	bool jump = false;
	bool lstep = false;
	
	while( true )
	{
		if ( false == jump )
		{
			dz = s/ns;
			if ( (z+dz) > ze )
			{
				dz=ze-z;
				if ( dz <= ep )
					return;
			}
			
			dzot=dz*.5;
			saoa(z+dzot,g3);
			saoa(z+dz,g5);	
		} /* if ( false == jump ) */
		
		bool nogo = false;
		for (int i = 0; i < n; i++ )
		{
			t00=(g1[i]+g5[i])*dzot;
			t01[i]=(t00+dz*g3[i])*.5;
			t10[i]=(4.*t01[i]-t00)/3.;
			
			/* test convergence of 3 point romberg result */
			test( real(t01[i]), real(t10[i]), &tr, imag(t01[i]), imag(t10[i]), &ti, 0. );
			if ( (tr > CRIT) || (ti > CRIT) )
				nogo = true;
		}
		
		if ( false == nogo )
		{
			for (int i = 0; i < n; i++ )
				sum[i] += t10[i];
			
			nt += 2;
			z += dz;
			if (z > zend)
				return;
			
			for (int i = 0; i < n; i++ )
				g1[i]=g5[i];
			
			if ( (nt >= NTS) && (ns > nx) )
			{
				ns=ns/2;
				nt=1;
			}
			
			jump = false;
			continue;
		} /* if ( false == nogo ) */
		
		saoa(z+dz*.25,g2);
		saoa(z+dz*.75,g4);
		nogo = false;
		for (int i = 0; i < n; i++ )
		{
			t02=(t01[i]+dzot*(g2[i]+g4[i]))*.5;
			t11=(4.*t02-t01[i])/3.;
			t20[i]=(16.*t11-t10[i])/15.;
			
			/* test convergence of 5 point Romberg result */
			test( real(t11), real(t20[i]), &tr, imag(t11), imag(t20[i]), &ti, 0.0 );
			if ( (tr > CRIT) || (ti > CRIT) )
				nogo = true;
		}
		
		if ( false == nogo )
		{
			for (int i = 0; i < n; i++ )
				sum[i] += t20[i];
			
			nt++;
			z += dz;
			if (z > zend)
				return;
			
			for (int i = 0; i < n; i++ )
				g1[i]=g5[i];
			
			if ( (nt >= NTS) && (ns > nx) )
			{
				ns=ns/2;
				nt=1;
			}
			
			jump = false;
			continue;	
		} /* if ( false == nogo ) */
		
		nt=0;
		if (ns < NM)
		{
			ns *= 2;
			dz=s/ns;
			dzot=dz*.5;
			
			for (int i = 0; i < n; i++ )
			{
				g5[i]=g3[i];
				g3[i]=g2[i];
			}
			
			jump = true;
			continue;	
		} /* if (ns < NM) */
		
		if ( false == lstep )
		{
			lstep = true;
			lambda( z, &t00, &t11 );
		}
		
		for (int i = 0; i < n; i++ )
			sum[i] += t20[i];
		
		nt++;
		z += dz;
		if (z > zend)
			return;
		
		for (int i = 0; i < n; i++ )
			g1[i]=g5[i];
		
		if ( (nt >= NTS) && (ns > nx) )
		{
			ns /= 2;
			nt=1;
		}
		
		jump = false;
	} /* while( TRUE ) */
}

/*! \brief saoa computes the integrand for each of the 6 Sommerfeld
	integrals for source and observer above ground. */
void c_evlcom::saoa( nec_float t, complex_array& ans)
{
	static nec_complex xl, dxl, cgam1, cgam2, b0, b0p, com, dgam, den1, den2;
	
	lambda(t, &xl, &dxl);
	if ( m_bessel_flag == true )
	{
		/* Bessel function form */
		bessel(xl*m_rho, &b0, &b0p);
		b0  *=2.;
		b0p *=2.;
		cgam1=sqrt(xl*xl-m_ck1sq);
		cgam2=sqrt(xl*xl-m_ck2sq);
		if (real(cgam1) == 0.0)
			cgam1=nec_complex(0.0,-fabs(imag(cgam1)));
		if (real(cgam2) == 0.)
			cgam2=nec_complex(0.0,-fabs(imag(cgam2)));
	}
	else
	{
		/* Hankel function form */
		hankel(xl*m_rho, &b0, &b0p);
		com=xl-m_ck1;
		cgam1=sqrt(xl+m_ck1)*sqrt(com);
		if (real(com) < 0. && imag(com) >= 0.)
			cgam1=-cgam1;
		com=xl-m_ck2;
		cgam2=sqrt(xl+m_ck2)*sqrt(com);
		if (real(com) < 0. && imag(com) >= 0.)
			cgam2=-cgam2;
	}
	
	if (norm(xl) >= m_tsmag)
	{
		if (imag(xl) >= 0.0)
		{
			nec_float xlr = real(xl);
			if (xlr >= m_ck2)
			{
				if (xlr <= m_ck1r)
					dgam=cgam2-cgam1;
				else
				{
					nec_float sign =1.0;
					dgam=1.0/(xl*xl);
					dgam=sign*((m_ct3*dgam+m_ct2)*dgam+m_ct1)/xl;
				}
			}
			else
			{
				nec_float sign=-1.0;
				dgam=1.0/(xl*xl);
				dgam=sign*((m_ct3*dgam+m_ct2)*dgam+m_ct1)/xl;
			} /* if (xlr >= m_ck2) */
		} /* if (imag(xl) >= 0.) */
		else
		{
			nec_float sign=1.0;
			dgam=1.0/(xl*xl);
			dgam=sign*((m_ct3*dgam+m_ct2)*dgam+m_ct1)/xl;
		}
	} /* if (norm(xl) < m_tsmag) */
	else
	{
		dgam=cgam2-cgam1;
	}
	
#if 0
	nec_float xlr = real(xl);
	if ( (xlr >= m_ck2) && (xlr <= m_ck1r))
	{
		dgam=cgam2-cgam1;
	}
	else
	{
		sign = 1.0;
		if ((imag(xl) >= 0.0) && (xlr < m_ck2))
			sign = -1.0;
		
		nec_floaf temp = 1.0/(xl*xl);
		dgam=sign*((m_ct3*temp+m_ct2)*temp+m_ct1)/xl;
	}
#endif	
	
	den2=m_cksm*dgam/(cgam2*(m_ck1sq*cgam2+m_ck2sq*cgam1));
	den1=1./(cgam1+cgam2)-m_cksm/cgam2;
	com=dxl*xl*exp(-cgam2*m_zph);
	ans[5]=com*b0*den1/m_ck1;
	com *= den2;
	
	if (m_rho != 0.)
	{
		b0p=b0p/m_rho;
		ans[0]=-com*xl*(b0p+b0*xl);
		ans[3]=com*xl*b0p;
	}
	else
	{
		ans[0]=-com*xl*xl*.5;
		ans[3]=ans[0];
	}
	
	ans[1]=com*cgam2*cgam2*b0;
	ans[2]=-ans[3]*cgam2*m_rho;
	ans[4]=com*b0;
}

/* evlua controls the integration contour in the complex */
/* lambda plane for evaluation of the sommerfeld integrals */
void c_evlcom::evlua( nec_complex *erv, nec_complex *ezv,
	nec_complex *erh, nec_complex *eph )
{
	static nec_float del, slope, rmis;
	static nec_complex cp1, cp2, cp3, bk, delta, delta2;
	
	complex_array sum(6), ans(6);
	
	del=m_zph;
	if ( m_rho > del )
		del=m_rho;
	
	if (m_zph >= 2.*m_rho)
	{
		/* Bessel function form of Sommerfeld integrals */
		m_bessel_flag=true;
		m_contour_a=nec_complex(0.0,0.0);
		del=1.0/del;
		
		if ( del > m_tkmag)
		{
			m_contour_b=nec_complex(0.1*m_tkmag,-0.1*m_tkmag);
			rom1(6,sum,2);
			m_contour_a=m_contour_b;
			m_contour_b=nec_complex(del,-del);
			rom1 (6,ans,2);
			for (int i = 0; i < 6; i++ )
				sum[i] += ans[i];
		}
		else
		{
			m_contour_b=nec_complex(del,-del);
			rom1(6,sum,2);
		}
		
		delta=PTP*del;
		gshank(m_contour_b,delta,ans,6,sum,0,m_contour_b,m_contour_b);
		ans[5] *= m_ck1;
		
		/* conjugate since nec uses exp(+jwt) */
		*erv=conj(m_ck1sq*ans[2]);
		*ezv=conj(m_ck1sq*(ans[1]+m_ck2sq*ans[4]));
		*erh=conj(m_ck2sq*(ans[0]+ans[5]));
		*eph=-conj(m_ck2sq*(ans[3]+ans[5]));
		
		return;	
	} /* if (m_zph >= 2.*m_rho) */
	
	/* Hankel function form of Sommerfeld integrals */
	m_bessel_flag=false;
	cp1=nec_complex(0.0,.4*m_ck2);
	cp2=nec_complex(.6*m_ck2,-.2*m_ck2);
	cp3=nec_complex(1.02*m_ck2,-.2*m_ck2);
	m_contour_a=cp1;
	m_contour_b=cp2;
	rom1(6,sum,2);
	m_contour_a=cp2;
	m_contour_b=cp3;
	rom1(6,ans,2);
	
	for (int i = 0; i < 6; i++ )
		sum[i]=-(sum[i]+ans[i]);
	
	/* path from imaginary axis to -infinity */
	if (m_zph > .001*m_rho)
		slope=m_rho/m_zph;
	else
		slope=1000.;
	
	del=PTP/del;
	delta=nec_complex(-1.0,slope)*del/sqrt(1.+slope*slope);
	delta2=-conj(delta);
	gshank(cp1,delta,ans,6,sum,0,bk,bk);
	rmis=m_rho*(real(m_ck1)-m_ck2);
	
	bool jump = false;
	if ( (rmis >= 2.*m_ck2) && (m_rho >= 1.e-10) )
	{
		if (m_zph >= 1.e-10)
		{
			bk=nec_complex(-m_zph,m_rho)*(m_ck1-cp3);
			rmis=-real(bk)/fabs(imag(bk));
			if (rmis > 4.*m_rho/m_zph)
				jump = true;
		}
		
		if ( false == jump )
		{
			/* integrate up between branch cuts, then to + infinity */
			cp1=m_ck1- nec_complex(0.1,+0.2);
			cp2=cp1+.2;
			bk=nec_complex(0.,del);
			gshank(cp1,bk,sum,6,ans,0,bk,bk);
			m_contour_a=cp1;
			m_contour_b=cp2;
			rom1(6,ans,1);
			for (int i = 0; i < 6; i++ )
				ans[i] -= sum[i];
			
			gshank(cp3,bk,sum,6,ans,0,bk,bk);
			gshank(cp2,delta2,ans,6,sum,0,bk,bk);
		}
		
		jump = true;	
	} /* if ( (rmis >= 2.*m_ck2) || (m_rho >= 1.e-10) ) */
	else
		jump = false;
	
	if ( false == jump )
	{
		/* integrate below branch points, then to + infinity */
		for (int i = 0; i < 6; i++ )
			sum[i]=-ans[i];
		
		rmis=real(m_ck1)*1.01;
		if ( (m_ck2+1.) > rmis )
			rmis=m_ck2+1.;
		
		bk=nec_complex(rmis,.99*imag(m_ck1));
		delta=bk-cp3;
		delta *= del/abs(delta);
		gshank(cp3,delta,ans,6,sum,1,bk,delta2);
	} /* if ( false == jump ) */
	
	ans[5] *= m_ck1;
	
	/* conjugate since nec uses exp(+jwt) */
	*erv=conj(m_ck1sq*ans[2]);
	*ezv=conj(m_ck1sq*(ans[1]+m_ck2sq*ans[4]));
	*erh=conj(m_ck2sq*(ans[0]+ans[5]));
	*eph=-conj(m_ck2sq*(ans[3]+ans[5]));
}

/*-----------------------------------------------------------------------*/

#define	GAMMA	.5772156649
#define C3	.7978845608
#define P10	.0703125
#define P20	.1121520996
#define Q10	.125
#define Q20	.0732421875
#define P11	.1171875
#define P21	.1441955566
#define Q11	.375
#define Q21	.1025390625
#define POF	.7853981635


/* bessel evaluates the zero-order bessel function */
/* and its derivative for complex argument z. */
void bessel( nec_complex z, nec_complex *j0, nec_complex *j0p )
{
	static int m[101];
	static nec_float a1[25], a2[25];
	static nec_complex cplx_01(0.0,1.0);
	static nec_complex cplx_10(1.0,0.0);
	
	/* initialization of constants */
	static bool bessel_init = false;
	
	if ( false == bessel_init )
	{
		for (int k = 1; k <= 25; k++ )
		{
			int index = k-1;
			a1[index] = -0.25/(k*k);
			a2[index] = 1.0/(k+1.0);
		}
		
		for (int i = 1; i <= 101; i++ )
		{
			nec_float tst=1.0;
			int init;
			for (int k = 0; k < 24; k++ )
			{
				init = k;
				tst *= -i*a1[k];
				if ( tst < 1.0e-6 )
					break;
			}
			
			m[i-1] = init+1;
		} /* for (int i = 1; i<= 101; i++ ) */
		
		bessel_init = true;
	} /* if (false == bessel_init) */
	
	nec_float zms = norm(z);
	
	if (zms <= 1.e-12)
	{
		*j0=cplx_10;
		*j0p=-0.5*z;
		return;
	}
	
	nec_complex j0x, j0px;
	int ib=0;
	if (zms <= 37.21)
	{
		if (zms > 36.0)
			ib=1;
		
		/* series expansion */
		#pragma message("Some strange code below. Why use the norm of a vector as an index?")
		int iz = int(zms); 	// TCAM : conversion of nec_float to int here! 
				// Using int() I think that this is the same as the fortran implicit coercion
				// but perhaps we should be doing an explicit rounding operation?
	
		int miz=m[iz];
		*j0 = cplx_10;
		*j0p = cplx_10;
		nec_complex zk = cplx_10;
		nec_complex zi = z*z;
		
		for (int k = 0; k < miz; k++ )
		{
			zk *= a1[k]*zi;
			*j0 += zk;
			*j0p += a2[k]*zk;
		}
		*j0p *= -0.5*z;
		
		if (ib == 0)
			return;
		
		j0x=*j0;
		j0px=*j0p;
	}
	
	/* asymptotic expansion */
	nec_complex zi = 1.0/z;
	nec_complex zi2 = zi*zi;
	nec_complex p0z = 1.0 + (P20*zi2-P10)*zi2;
	nec_complex p1z = 1.0 +(P11-P21*zi2)*zi2;
	nec_complex q0z = (Q20*zi2-Q10)*zi;
	nec_complex q1z = (Q11-Q21*zi2)*zi;
	nec_complex zk = exp(cplx_01 * (z-POF));
	
	zi2 = 1.0/zk;
	nec_complex cz = 0.5*(zk+zi2);
	nec_complex sz = cplx_01 * 0.5 * (zi2-zk);
	zk = C3*sqrt(zi);
	*j0 = zk*(p0z*cz-q0z*sz);
	*j0p = -zk*(p1z*sz+q1z*cz);
	
	if (ib == 0)
		return;
	
	nec_float pi_10 = pi() * 10.0;
	zms = cos((sqrt(zms)-6.0)*pi_10);
	*j0 = 0.5*(j0x*(1.0+zms)+ *j0*(1.0-zms));
	*j0p = 0.5*(j0px*(1.0+zms)+ *j0p*(1.0-zms));
}


#define C1	-.02457850915
#define C2	.3674669052

/* hankel evaluates hankel function of the first kind,   */
/* order zero, and its derivative for complex argument z */
void hankel( nec_complex z, nec_complex *h0, nec_complex *h0p )
{
	static int m[101];
	static nec_float a1[25], a2[25], a3[25], a4[25];
	nec_complex clogz, p0z, p1z, q0z, q1z, zi, zi2, zk;
	static nec_complex cplx_01(0.0,1.0);

	static bool hankel_init = false;
	/* initialization of constants */
	if ( ! hankel_init )
	{
		nec_float psi=-GAMMA;
		for (int k = 1; k <= 25; k++ )
		{
			int i = k-1;
			a1[i]=-.25/(k*k);
			a2[i]=1.0/(k+1.0);
			psi += 1.0/k;
			a3[i]=psi+psi;
			a4[i]=(psi+psi+1.0/(k+1.0))/(k+1.0);
		}

		for (int i = 1; i <= 101; i++ )
		{
			int init;
			nec_float test=1.0;
			for (int k = 0; k < 24; k++ )
			{
				init = k;
				test *= -i*a1[k];
				if ((test*a3[k]) < 1.e-6)
					break;
			}
			m[i-1]=init+1;
		}
		hankel_init = true;
	} /* if ( ! hankel_init ) */

	nec_float zms = norm(z);

	if (zms == 0.0)
		throw new nec_exception("hankel not valid for z=0.");

	nec_complex y0(0,0);
	nec_complex y0p(0,0);
	int ib=0;
	if (zms <= 16.81)
	{
		if (zms > 16.)
			ib=1;

		/* series expansion */
		int iz = int(zms); // TCAM using explicit int() coercion
		int miz = m[iz];
		nec_complex j0(1.0,0.0);
		nec_complex j0p(1.0,0.0);
		zk = j0;
		zi = z*z;

		for (int k = 0; k < miz; k++ )
		{
			zk *= a1[k]*zi;
			j0 += zk;
			j0p += a2[k]*zk;
			y0 += a3[k]*zk;
			y0p += a4[k]*zk;
		}

		j0p *= -.5*z;
		clogz= log(.5*z);
		y0 = (2.0 * j0 * clogz-y0)/pi() + C2;
		y0p = (2.0/z +2.0*j0p*clogz + 0.5*y0p*z)/pi() + C1*z;
		*h0=j0+cplx_01*y0;
		*h0p=j0p+cplx_01*y0p;

		if (ib == 0)
			return;

		y0=*h0;
		y0p=*h0p;
	} /* if (zms <= 16.81) */

	/* asymptotic expansion */
	zi=1./z;
	zi2=zi*zi;
	p0z=1.+(P20*zi2-P10)*zi2;
	p1z=1.+(P11-P21*zi2)*zi2;
	q0z=(Q20*zi2-Q10)*zi;
	q1z=(Q11-Q21*zi2)*zi;
	zk=exp(cplx_01*(z-POF))*sqrt(zi)*C3;
	*h0=zk*(p0z+cplx_01*q0z);
	*h0p=cplx_01*zk*(p1z+cplx_01*q1z);

	if (ib == 0)
		return;

	zms=cos((sqrt(zms)-4.)*31.41592654);
	*h0=.5*(y0*(1.+zms)+ *h0*(1.-zms));
	*h0p=.5*(y0p*(1.+zms)+ *h0p*(1.-zms));
}


