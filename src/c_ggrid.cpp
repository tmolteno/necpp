/*
    Copyright (C) 2004-2008    Timothy C.A. Molteno
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA    02111-1307    USA
*/
#include "c_ggrid.h"
#include "common.h"
#include "electromag.h"
#include "nec_output.h" // for DEBUG_TRACE()

#include <cstdlib>

//#define    CONST4    nec_complex(0.0,em::impedance() / 2.0)

int    c_ggrid::m_nxa[3] = {11,17,9};
int    c_ggrid::m_nya[3] = {10,5,8};

nec_float c_ggrid::m_dxa[3] = {.02,.05,.1}; 
nec_float c_ggrid::m_dya[3] = {.1745329252,.0872664626,.1745329252};

nec_float c_ggrid::m_xsa[3] = {0.,.2,.2};
nec_float c_ggrid::m_ysa[3] = {0.,0.,.3490658504};

/*! \brief interpolate (was intrp) uses bivariate cubic interpolation to obtain the values of 4 functions at the point (x,y).
*/
void c_ggrid::interpolate( nec_float x, nec_float y, nec_complex *f1,
        nec_complex *f2, nec_complex *f3, nec_complex *f4 )
{
    static const int nda[3] = {11,17,9}, ndpa[3] = {110, 85, 72};

    bool recalculate = true;

    if( (x >= m_ip_xs) && (y >= m_ip_ys) ) {
        m_ip_ix = (int)((x - m_ip_xs) / m_ip_dx) + 1;
        m_ip_iy = (int)((y - m_ip_ys) / m_ip_dy) + 1;
    } else {
        /* if point lies in same 4 by 4 point region */
        /* as previous point, old values are reused. */
        if ( ((m_ip_ix >= m_ip_ixeg) && (m_ip_iy >= m_ip_iyeg)) &&
             ((std::abs(m_ip_ix - m_ip_ixs) < 2) &&  (std::abs(m_ip_iy - m_ip_iys) < 2)) )
            recalculate = false;
    }

    if (true == recalculate) {
        /* determine correct grid and grid region */
        int igr;

        if( x <= m_xsa[1])
            igr=0;
        else {
            if( y > m_ysa[2])
                igr=2;
            else
                igr=1;
        }

        if( igr != m_ip_igrs) {
            m_ip_igrs = igr;
            m_ip_dx = m_dxa[m_ip_igrs];
            m_ip_dy = m_dya[m_ip_igrs];
            m_ip_xs = m_xsa[m_ip_igrs];
            m_ip_ys = m_ysa[m_ip_igrs];
            m_ip_nxm2 = m_nxa[m_ip_igrs] - 2;
            m_ip_nym2 = m_nya[m_ip_igrs] - 2;
            m_ip_nxms = (( m_ip_nxm2 + 1) / 3) * 3 + 1;
            m_ip_nyms = (( m_ip_nym2 + 1) / 3) * 3 + 1;
            m_ip_nd = nda[m_ip_igrs];
            m_ip_ndp = ndpa[m_ip_igrs];
            m_ip_ix = (int)(( x - m_ip_xs) / m_ip_dx) + 1;
            m_ip_iy = (int)(( y - m_ip_ys) / m_ip_dy) + 1;
        } /* if( igr != m_ip_igrs) */

        m_ip_ixs = (( m_ip_ix - 1) / 3) * 3 + 2;
        if( m_ip_ixs < 2)
            m_ip_ixs = 2;
        m_ip_ixeg = -10000;

        if( m_ip_ixs > m_ip_nxm2) {
            m_ip_ixs = m_ip_nxm2;
            m_ip_ixeg = m_ip_nxms;
        }

        m_ip_iys = (( m_ip_iy - 1) / 3) * 3 + 2;
        if( m_ip_iys < 2)
            m_ip_iys = 2;
        m_ip_iyeg = -10000;

        if( m_ip_iys > m_ip_nym2) {
            m_ip_iys = m_ip_nym2;
            m_ip_iyeg = m_ip_nyms;
        }

        // Select the correct array once, outside the inner loop.
        complex_array& ar = (m_ip_igrs == 0) ? m_ar1 : (m_ip_igrs == 1) ? m_ar2 : m_ar3;

        /* compute coefficients of 4 cubic polynomials in x for */
        /* the 4 grid values of y for each of the 4 functions */
        int iadz = m_ip_ixs + (m_ip_iys - 3) * m_ip_nd - m_ip_ndp;
        for (int k = 0; k < 4; k++ ) {
            iadz += m_ip_ndp;
            int iadd = iadz;

            for (int i = 0; i < 4; i++ ) {
                iadd += m_ip_nd;

                nec_complex p1 = ar[iadd-2];
                nec_complex p2 = ar[iadd-1];
                nec_complex p3 = ar[iadd];
                nec_complex p4 = ar[iadd+1];

                m_ip_a[i][k] = ( p4 - p1 + 3.0 * ( p2 - p3)) * 0.1666666667;
                m_ip_b[i][k] = ( p1 - 2.0 * p2 + p3) * 0.5;
                m_ip_c[i][k] = p3 - (2.0 * p1 + 3.0 * p2 + p4) * 0.1666666667;
                m_ip_d[i][k] = p2;

            } /* for ( i = 0; i < 4; i++ ) */

        } /* for ( k = 0; k < 4; k++ ) */

        m_ip_xz = ( m_ip_ixs - 1) * m_ip_dx + m_ip_xs;
        m_ip_yz = ( m_ip_iys - 1) * m_ip_dy + m_ip_ys;

    } /* if (true == recalculate) */

    /* evaluate polynomials in x and use cubic */
    /* interpolation in y for each of the 4 functions. */
    nec_float xx = ( x - m_ip_xz) / m_ip_dx;
    nec_float yy = ( y - m_ip_yz) / m_ip_dy;

    // Helper lambda to evaluate cubic for one function index k
    auto eval_func = [&](int k, nec_complex* fout) {
        nec_complex l_fx1 = ((m_ip_a[0][k] * xx + m_ip_b[0][k]) * xx + m_ip_c[0][k]) * xx + m_ip_d[0][k];
        nec_complex l_fx2 = ((m_ip_a[1][k] * xx + m_ip_b[1][k]) * xx + m_ip_c[1][k]) * xx + m_ip_d[1][k];
        nec_complex l_fx3 = ((m_ip_a[2][k] * xx + m_ip_b[2][k]) * xx + m_ip_c[2][k]) * xx + m_ip_d[2][k];
        nec_complex l_fx4 = ((m_ip_a[3][k] * xx + m_ip_b[3][k]) * xx + m_ip_c[3][k]) * xx + m_ip_d[3][k];
        nec_complex l_p1 = l_fx4 - l_fx1 + 3.0 * (l_fx2 - l_fx3);
        nec_complex l_p2 = 3.0 * (l_fx1 - 2.0 * l_fx2 + l_fx3);
        nec_complex l_p3 = 6.0 * l_fx3 - 2.0 * l_fx1 - 3.0 * l_fx2 - l_fx4;
        *fout = ((l_p1 * yy + l_p2) * yy + l_p3) * yy * 0.1666666667 + l_fx2;
    };

    eval_func(0, f1);
    eval_func(1, f2);
    eval_func(2, f3);
    eval_func(3, f4);
}

#include "electromag.h"

/*! was SOMNEC in the original NEC-2 source code
*/
void c_ggrid::sommerfeld( nec_float epr, nec_float sig, nec_float wavelength )
{
    static nec_complex const1_neg = - nec_complex(0.0,4.771341189);
    static nec_complex CONST4(0.0, em::impedance() / 2.0);
    
    nec_float dr, dth, r, rk, thet, tfac1, tfac2;
    nec_complex erv, ezv, erh, eph, cl1, cl2, con;
    
    if(sig >= 0.0) {
        m_epscf = nec_complex(epr,-sig*wavelength*em::impedance_over_2pi());
    } else {
        m_epscf=nec_complex(epr,sig);
    }
    m_evlcom.m_ck2 = two_pi();
    m_evlcom.m_ck2sq = m_evlcom.m_ck2*m_evlcom.m_ck2;
    
    /*
    Sommerfeld integral evaluation uses exp(-jwt), NEC uses exp(+jwt),
    hence need conjg(epscf).    Conjugate of fields occurs in subroutine 
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
    for (int k = 0; k < 3; k++ ) {
        int nr = m_nxa[k];
        int nth = m_nya[k];
        dr = m_dxa[k];
        dth = m_dya[k];
        r = m_xsa[k]-dr;
        int irs=1;
        
        if(k == 0) {
            r=m_xsa[k];
            irs=2;
        }
    
        /*    loop over r.    (r=sqrt(m_evlcom.m_rho**2 + (z+h)**2)) */
        for (int ir = irs-1; ir < nr; ir++ ) {
            r += dr;
            thet = m_ysa[k]-dth;
        
            /* loop over theta.    (theta=atan((z+h)/m_evlcom.m_rho)) */
            for (int ith = 0; ith < nth; ith++ ) {
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
            
                switch( k ) {
                case 0:
                    m_ar1[ir+ith*11+    0]=erv*con;
                    m_ar1[ir+ith*11+110]=ezv*con;
                    m_ar1[ir+ith*11+220]=erh*con;
                    m_ar1[ir+ith*11+330]=eph*con;
                    break;
            
                case 1:
                    m_ar2[ir+ith*17+    0]=erv*con;
                    m_ar2[ir+ith*17+ 85]=ezv*con;
                    m_ar2[ir+ith*17+170]=erh*con;
                    m_ar2[ir+ith*17+255]=eph*con;
                    break;
            
                case 2:
                    m_ar3[ir+ith*9+    0]=erv*con;
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
    thet = -dth;
    int nth = m_nya[0];
    
    for (int ith = 0; ith < nth; ith++ ) {
        thet += dth;
        if( (ith+1) != nth ) {
            tfac2=cos(thet);
            tfac1=(1.-sin(thet))/tfac2;
            tfac2=tfac1/tfac2;
            erv=m_epscf*cl1*tfac1;
            erh=cl1*(tfac2-1.)+cl2;
            eph=cl1*tfac2-cl2;
        } else {
            erv=0.;
            erh=cl2-.5*cl1;
            eph = -erh;
        }
    
        m_ar1[0+ith*11+    0]=erv;
        m_ar1[0+ith*11+110]=ezv;
        m_ar1[0+ith*11+220]=erh;
        m_ar1[0+ith*11+330]=eph;
    }
}



/* fbar is the Sommerfeld attenuation function for numerical distance . */
nec_complex    fbar(const nec_complex& p );
nec_complex    fbar(const nec_complex& p )    {
    static nec_float TOSP = 2.0 / sqrt_pi();

    int minus;
    nec_float tms, sms;
    nec_complex z, zs, sum, pow, term, fbar;
    
    z= cplx_01()* sqrt( p);
    if ( abs( z) <= 3.) {
        /* series expansion */
        zs= z* z;
        sum= z;
        pow= z;
    
        for (int i = 1; i <= 100; i++ ) {
            pow = - pow* zs/ (nec_float)i;
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
    if ( real( z) < 0.) {
        minus=1;
        z = - z;
    }
    else
        minus=0;
    
    zs=.5/( z* z);
    sum=cplx_00();
    term=cplx_10();
    
    for (int i = 1; i <= 6; i++ ) {
        term = - term*(2.*i -1.)* zs;
        sum += term;
    }
    
    if ( minus == 1)
        sum -= 2.0 * sqrt_pi() * z* exp( z* z);
    fbar = - sum;
    
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
    rk1 = - two_pi_j()* ground_wave.r1;
    rk2 = - two_pi_j()* ground_wave.r2;
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
    erv = -( x1+ x2- x3+ x4- x5+ x6- x7)* (-CONST4);

    ezh = -( x1- x2+ x3- x4- x5- x6+ x7)* (-CONST4);

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
    eph = -( x1- x2+ x3- x4+ x5+ x6)* (-CONST4);
}

