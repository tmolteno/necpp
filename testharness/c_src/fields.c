/******* Translated to the C language by N. Kyriazis  20 Aug 2003 *******/
/*									*/
/* Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,	*/
/* tape15,tape16,tape20,tape21)						*/
/*									*/
/* Numerical Electromagnetics Code (NEC2)  developed at Lawrence	*/
/* Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414	*/
/* for problems with the NEC code. For problems with the vax implem- 	*/
/* entation, contact J. Breakall at 415-422-8196 or E. Domning at 415 	*/
/* 422-5936) 								*/
/* file created 4/11/80. 						*/
/*									*/
/*                ***********Notice********** 				*/
/* This computer code material was prepared as an account of work 	*/
/* sponsored by the United States government.  Neither the United 	*/
/* States nor the United States Department Of Energy, nor any of 	*/
/* their employees, nor any of their contractors, subcontractors, 	*/
/* or their employees, makes any warranty, express or implied, or	*/
/* assumes any legal liability or responsibility for the accuracy, 	*/
/* completeness or usefulness of any information, apparatus, product 	*/
/* or process disclosed, or represents that its use would not infringe 	*/
/* privately-owned rights. 						*/
/*									*/
/************************************************************************/

#include "nec2c.h"

/* common  /dataj/ */
extern int iexk, ind1, indd1, ind2, indd2, ipgnd;
extern long double s, b, xj, yj, zj, cabj, sabj, salpj, rkh;
extern long double t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
extern complex long double  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;

/* common  /gnd/ */
extern int ksymp, ifar, iperf, nradl;
extern long double t2, cl, ch, scrwl, scrwr;
extern complex long double zrati, zrati2, t1, frati;

/* common  /incom/ */
extern int isnor;
extern long double xo, yo, zo, sn, xsn, ysn;

/* common  /tmi/ */
extern int ija; /* changed to ija to avoid conflict */
extern long double zpk, rkb2;

/*common  /tmh/ */
extern long double zpka, rhks;

/* common  /gwav/ */
extern long double r1, r2, zmh, zph;
extern complex long double u, u2, xx1, xx2;

/* common  /data/ */
extern int n, np, m, mp, ipsym, npm, np2m, np3m; /* n+m,n+2m,n+3m */
extern int *icon1, *icon2, *itag;
extern long double *x, *y, *z, *si, *bi;
extern long double *x2, *y2, *z2, *cab, *sab, *salp;
extern long double *t1x, *t1y, *t1z, *t2x, *t2y, *t2z;
extern long double *px, *py, *pz, *pbi, *psalp;
extern long double wlam;

/* common  /crnt/ */
extern long double *air, *aii, *bir, *bii, *cir, *cii;
extern complex long double *cur;

/* common  /fpat/ */
extern int near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp, inor, iax, ixtyp;
extern long double thets, phis, dth, dph, rfld, gnor, clt, cht, epsr2, sig2;
extern long double xpr6, pinr, pnlr, ploss, xnr, ynr, znr, dxnr, dynr, dznr;

/* common  /plot/ */
extern int iplp1, iplp2, iplp3, iplp4;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/*-------------------------------------------------------------------*/

/* compute near e fields of a segment with sine, cosine, and */
/* constant currents.  ground effect included. */
void efld( long double xi, long double yi, long double zi, long double ai, int ij )
{
#define	txk	egnd[0]
#define	tyk	egnd[1]
#define	tzk	egnd[2]
#define	txs	egnd[3]
#define	tys	egnd[4]
#define	tzs	egnd[5]
#define	txc	egnd[6]
#define	tyc	egnd[7]
#define	tzc	egnd[8]

  int ip;
  long double xij, yij, ijx, rfl, salpr, zij, zp, rhox;
  long double rhoy, rhoz, rh, r, rmag, cth, px, py;
  long double xymag, xspec, yspec, rhospc, dmin, shaf;
  complex long double epx, epy, refs, refps, zrsin, zratx, zscrn;
  complex long double tezs, ters, tezc, terc, tezk, terk, egnd[9];

  xij= xi- xj;
  yij= yi- yj;
  ijx= ij;
  rfl=-1.;

  for( ip = 0; ip < ksymp; ip++ )
  {
    if( ip == 1)
      ijx=1;
    rfl=- rfl;
    salpr= salpj* rfl;
    zij= zi- rfl* zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;

    rh= sqrtl( rhox* rhox+ rhoy* rhoy+ rhoz* rhoz+ ai* ai);
    if( rh <= 1.e-10)
    {
      rhox=0.;
      rhoy=0.;
      rhoz=0.;
    }
    else
    {
      rhox= rhox/ rh;
      rhoy= rhoy/ rh;
      rhoz= rhoz/ rh;
    }

    /* lumped current element approx. for large separations */
    r= sqrtl( zp* zp+ rh* rh);
    if( r >= rkh)
    {
      rmag= TP* r;
      cth= zp/ r;
      px= rh/ r;
      txk= cmplx( cosl( rmag),- sinl( rmag));
      py= TP* r* r;
      tyk= ETA* cth* txk* cmplx(1.0,-1.0/ rmag)/ py;
      tzk= ETA* px* txk* cmplx(1.0, rmag-1.0/ rmag)/(2.* py);
      tezk= tyk* cth- tzk* px;
      terk= tyk* px+ tzk* cth;
      rmag= sinl( PI* s)/ PI;
      tezc= tezk* rmag;
      terc= terk* rmag;
      tezk= tezk* s;
      terk= terk* s;
      txs=CPLX_00;
      tys=CPLX_00;
      tzs=CPLX_00;

    } /* if( r >= rkh) */

    if( r < rkh)
    {
      /* eksc for thin wire approx. or ekscx for extended t.w. approx. */
      if( iexk != 1)
	eksc( s, zp, rh, TP, ijx, &tezs, &ters,
	    &tezc, &terc, &tezk, &terk );
      else
	ekscx( b, s, zp, rh, TP, ijx, ind1, ind2,
	    &tezs, &ters, &tezc, &terc, &tezk, &terk);

      txs= tezs* cabj+ ters* rhox;
      tys= tezs* sabj+ ters* rhoy;
      tzs= tezs* salpr+ ters* rhoz;

    } /* if( r < rkh) */

    txk= tezk* cabj+ terk* rhox;
    tyk= tezk* sabj+ terk* rhoy;
    tzk= tezk* salpr+ terk* rhoz;
    txc= tezc* cabj+ terc* rhox;
    tyc= tezc* sabj+ terc* rhoy;
    tzc= tezc* salpr+ terc* rhoz;

    if( ip == 1)
    {
      if( iperf <= 0)
      {
	zratx= zrati;
	rmag= r;
	xymag= sqrtl( xij* xij+ yij* yij);

	/* set parameters for radial wire ground screen. */
	if( nradl != 0)
	{
	  xspec=( xi* zj+ zi* xj)/( zi+ zj);
	  yspec=( yi* zj+ zi* yj)/( zi+ zj);
	  rhospc= sqrtl( xspec* xspec+ yspec* yspec+ t2* t2);

	  if( rhospc <= scrwl)
	  {
	    zscrn= t1* rhospc* logl( rhospc/ t2);
	    zratx=( zscrn* zrati)/( ETA* zrati+ zscrn);
	  }
	} /* if( nradl != 0) */

	/* calculation of reflection coefficients when ground is specified. */
	if( xymag <= 1.0e-6)
	{
	  px=0.;
	  py=0.;
	  cth=1.;
	  zrsin=CPLX_10;
	}
	else
	{
	  px=- yij/ xymag;
	  py= xij/ xymag;
	  cth= zij/ rmag;
	  zrsin= csqrtl(1.0 - zratx*zratx*(1.0 - cth*cth) );

	} /* if( xymag <= 1.0e-6) */

	refs=( cth- zratx* zrsin)/( cth+ zratx* zrsin);
	refps=-( zratx* cth- zrsin)/( zratx* cth+ zrsin);
	refps= refps- refs;
	epy= px* txk+ py* tyk;
	epx= px* epy;
	epy= py* epy;
	txk= refs* txk+ refps* epx;
	tyk= refs* tyk+ refps* epy;
	tzk= refs* tzk;
	epy= px* txs+ py* tys;
	epx= px* epy;
	epy= py* epy;
	txs= refs* txs+ refps* epx;
	tys= refs* tys+ refps* epy;
	tzs= refs* tzs;
	epy= px* txc+ py* tyc;
	epx= px* epy;
	epy= py* epy;
	txc= refs* txc+ refps* epx;
	tyc= refs* tyc+ refps* epy;
	tzc= refs* tzc;

      } /* if( iperf <= 0) */

      exk= exk- txk* frati;
      eyk= eyk- tyk* frati;
      ezk= ezk- tzk* frati;
      exs= exs- txs* frati;
      eys= eys- tys* frati;
      ezs= ezs- tzs* frati;
      exc= exc- txc* frati;
      eyc= eyc- tyc* frati;
      ezc= ezc- tzc* frati;
      continue;

    } /* if( ip == 1) */

    exk= txk;
    eyk= tyk;
    ezk= tzk;
    exs= txs;
    eys= tys;
    ezs= tzs;
    exc= txc;
    eyc= tyc;
    ezc= tzc;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  if( iperf != 2)
    return;

  /* field due to ground using sommerfeld/norton */
  sn= sqrtl( cabj* cabj+ sabj* sabj);
  if( sn >= 1.0e-5)
  {
    xsn= cabj/ sn;
    ysn= sabj/ sn;
  }
  else
  {
    sn=0.;
    xsn=1.;
    ysn=0.;
  }

  /* displace observation point for thin wire approximation */
  zij= zi+ zj;
  salpr=- salpj;
  rhox= sabj* zij- salpr* yij;
  rhoy= salpr* xij- cabj* zij;
  rhoz= cabj* yij- sabj* xij;
  rh= rhox* rhox+ rhoy* rhoy+ rhoz* rhoz;

  if( rh <= 1.e-10)
  {
    xo= xi- ai* ysn;
    yo= yi+ ai* xsn;
    zo= zi;
  }
  else
  {
    rh= ai/ sqrtl( rh);
    if( rhoz < 0.)
      rh=- rh;
    xo= xi+ rh* rhox;
    yo= yi+ rh* rhoy;
    zo= zi+ rh* rhoz;

  } /* if( rh <= 1.e-10) */

  r= xij* xij+ yij* yij+ zij* zij;
  if( r <= .95)
  {
    /* field from interpolation is integrated over segment */
    isnor=1;
	dmin= exk* conjl( exk)+ eyk* conjl( eyk)+ ezk* conjl( ezk);
    dmin=.01* sqrtl( dmin);
    shaf=.5* s;
    rom2(- shaf, shaf, egnd, dmin);
  }
  else
  {
    /* norton field equations and lumped current element approximation */
    isnor=2;
    sflds(0., egnd);
  } /* if( r <= .95) */

  if( r > .95)
  {
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rh= r- zp* zp;
    if( rh <= 1.e-10)
      dmin=0.;
    else
      dmin= sqrtl( rh/( rh+ ai* ai));

    if( dmin <= .95)
    {
      px=1.- dmin;
      terk=( txk* cabj+ tyk* sabj+ tzk* salpr)* px;
      txk= dmin* txk+ terk* cabj;
      tyk= dmin* tyk+ terk* sabj;
      tzk= dmin* tzk+ terk* salpr;
      ters=( txs* cabj+ tys* sabj+ tzs* salpr)* px;
      txs= dmin* txs+ ters* cabj;
      tys= dmin* tys+ ters* sabj;
      tzs= dmin* tzs+ ters* salpr;
      terc=( txc* cabj+ tyc* sabj+ tzc* salpr)* px;
      txc= dmin* txc+ terc* cabj;
      tyc= dmin* tyc+ terc* sabj;
      tzc= dmin* tzc+ terc* salpr;

    } /* if( dmin <= .95) */

  } /* if( r > .95) */

  exk= exk+ txk;
  eyk= eyk+ tyk;
  ezk= ezk+ tzk;
  exs= exs+ txs;
  eys= eys+ tys;
  ezs= ezs+ tzs;
  exc= exc+ txc;
  eyc= eyc+ tyc;
  ezc= ezc+ tzc;

  return;
}

/*-----------------------------------------------------------------------*/

/* compute e field of sine, cosine, and constant */
/* current filaments by thin wire approximation. */
void eksc( long double s, long double z, long double rh, long double xk, int ij,
    complex long double *ezs, complex long double *ers, complex long double *ezc,
    complex long double *erc, complex long double *ezk, complex long double *erk )
{
  long double rhk, sh, shk, ss, cs, z1a, z2a, cint, sint;
  complex long double gz1, gz2, gp1, gp2, gzp1, gzp2;

  ija= ij;
  zpk= xk* z;
  rhk= xk* rh;
  rkb2= rhk* rhk;
  sh=.5* s;
  shk= xk* sh;
  ss= sinl( shk);
  cs= cosl( shk);
  z2a= sh- z;
  z1a=-( sh+ z);
  gx( z1a, rh, xk, &gz1, &gp1);
  gx( z2a, rh, xk, &gz2, &gp2);
  gzp1= gp1* z1a;
  gzp2= gp2* z2a;
  *ezs=  CONST1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *ezc=- CONST1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *erk= CONST1*( gp2- gp1)* rh;
  intx(- shk, shk, rhk, ij, &cint, &sint);
  *ezk=- CONST1*( gzp2- gzp1+ xk* xk* cmplx( cint,- sint));
  gzp1= gzp1* z1a;
  gzp2= gzp2* z2a;

  if( rh >= 1.0e-10)
  {
    *ers=- CONST1*(( gzp2+ gzp1+ gz2+ gz1)*
	ss-( z2a* gz2- z1a* gz1)* cs*xk)/ rh;
    *erc=- CONST1*(( gzp2- gzp1+ gz2- gz1)*
	cs+( z2a* gz2+ z1a* gz1)* ss*xk)/ rh;
    return;
  }

  *ers = CPLX_00;
  *erc = CPLX_00;

  return;
}

/*-----------------------------------------------------------------------*/

/* compute e field of sine, cosine, and constant current */
/* filaments by extended thin wire approximation. */
void ekscx( long double bx, long double s, long double z,
    long double rhx, long double xk, int ij, int inx1, int inx2,
    complex long double *ezs, complex long double *ers, complex long double *ezc,
    complex long double *erc, complex long double *ezk, complex long double *erk )
{
  int ira;
  long double b, rh, sh, rhk, shk, ss, cs, z1a;
  long double z2a, a2, bk, bk2, cint, sint;
  complex long double gz1, gz2, gzp1, gzp2, gr1, gr2;
  complex long double grp1, grp2, grk1, grk2, gzz1, gzz2;

  if( rhx >= bx)
  {
    rh= rhx;
    b= bx;
    ira=0;
  }
  else
  {
    rh= bx;
    b= rhx;
    ira=1;
  }

  sh=.5* s;
  ija= ij;
  zpk= xk* z;
  rhk= xk* rh;
  rkb2= rhk* rhk;
  shk= xk* sh;
  ss= sinl( shk);
  cs= cosl( shk);
  z2a= sh- z;
  z1a=-( sh+ z);
  a2= b* b;

  if( inx1 != 2)
    gxx( z1a, rh, b, a2, xk, ira, &gz1,
	&gzp1, &gr1, &grp1, &grk1, &gzz1);
  else
  {
    gx( z1a, rhx, xk, &gz1, &grk1);
    gzp1= grk1* z1a;
    gr1= gz1/ rhx;
    grp1= gzp1/ rhx;
    grk1= grk1* rhx;
    gzz1= CPLX_00;
  }

  if( inx2 != 2)
    gxx( z2a, rh, b, a2, xk, ira, &gz2,
	&gzp2, &gr2, &grp2, &grk2, &gzz2);
  else
  {
    gx( z2a, rhx, xk, &gz2, &grk2);
    gzp2= grk2* z2a;
    gr2= gz2/ rhx;
    grp2= gzp2/ rhx;
    grk2= grk2* rhx;
    gzz2= CPLX_00;
  }

  *ezs= CONST1*(( gz2- gz1)* cs* xk-( gzp2+ gzp1)* ss);
  *ezc=- CONST1*(( gz2+ gz1)* ss* xk+( gzp2- gzp1)* cs);
  *ers=- CONST1*(( z2a* grp2+ z1a* grp1+ gr2+ gr1)*ss
      -( z2a* gr2- z1a* gr1)* cs* xk);
  *erc=- CONST1*(( z2a* grp2- z1a* grp1+ gr2- gr1)*cs
      +( z2a* gr2+ z1a* gr1)* ss* xk);
  *erk= CONST1*( grk2- grk1);
  intx(- shk, shk, rhk, ij, &cint, &sint);
  bk= b* xk;
  bk2= bk* bk*.25;
  *ezk=- CONST1*( gzp2- gzp1+ xk* xk*(1.- bk2)*
      cmplx( cint,- sint)-bk2*( gzz2- gzz1));

  return;
}

/*-----------------------------------------------------------------------*/

/* integrand for h field of a wire */
void gh( long double zk, long double *hr, long double *hi)
{
  long double rs, r, ckr, skr, rr2, rr3;

  rs= zk- zpka;
  rs= rhks+ rs* rs;
  r= sqrtl( rs);
  ckr= cosl( r);
  skr= sinl( r);
  rr2=1./ rs;
  rr3= rr2/ r;
  *hr= skr* rr2+ ckr* rr3;
  *hi= ckr* rr2- skr* rr3;

  return;
}

/*-----------------------------------------------------------------------*/

/* gwave computes the electric field, including ground wave, of a */
/* current element over a ground plane using formulas of k.a. norton */
/* (proc. ire, sept., 1937, pp.1203,1236.) */

void gwave( complex long double *erv, complex long double *ezv,
    complex long double *erh, complex long double *ezh, complex long double *eph )
{
  long double sppp, sppp2, cppp2, cppp, spp, spp2, cpp2, cpp;
  complex long double rk1, rk2, t1, t2, t3, t4, p1, rv;
  complex long double omr, w, f, q1, rh, v, g, xr1, xr2;
  complex long double x1, x2, x3, x4, x5, x6, x7;

  sppp= zmh/ r1;
  sppp2= sppp* sppp;
  cppp2=1.- sppp2;

  if( cppp2 < 1.0e-20)
    cppp2=1.0e-20;

  cppp= sqrtl( cppp2);
  spp= zph/ r2;
  spp2= spp* spp;
  cpp2=1.- spp2;

  if( cpp2 < 1.0e-20)
    cpp2=1.0e-20;

  cpp= sqrtl( cpp2);
  rk1=- TPJ* r1;
  rk2=- TPJ* r2;
  t1=1. -u2* cpp2;
  t2= csqrtl( t1);
  t3=(1. -1./ rk1)/ rk1;
  t4=(1. -1./ rk2)/ rk2;
  p1= rk2* u2* t1/(2.* cpp2);
  rv=( spp- u* t2)/( spp+ u* t2);
  omr=1.- rv;
  w=1./ omr;
  w=(4.0 + 0.0fj)* p1* w* w;
  f= fbar( w);
  q1= rk2* t1/(2.* u2* cpp2);
  rh=( t2- u* spp)/( t2+ u* spp);
  v=1./(1.+ rh);
  v=(4.0 + 0.0fj)* q1* v* v;
  g= fbar( v);
  xr1= xx1/ r1;
  xr2= xx2/ r2;
  x1= cppp2* xr1;
  x2= rv* cpp2* xr2;
  x3= omr* cpp2* f* xr2;
  x4= u* t2* spp*2.* xr2/ rk2;
  x5= xr1* t3*(1.-3.* sppp2);
  x6= xr2* t4*(1.-3.* spp2);
  *ezv=( x1+ x2+ x3- x4- x5- x6)* (-CONST4);
  x1= sppp* cppp* xr1;
  x2= rv* spp* cpp* xr2;
  x3= cpp* omr* u* t2* f* xr2;
  x4= spp* cpp* omr* xr2/ rk2;
  x5=3.* sppp* cppp* t3* xr1;
  x6= cpp* u* t2* omr* xr2/ rk2*.5;
  x7=3.* spp* cpp* t4* xr2;
  *erv=-( x1+ x2- x3+ x4- x5+ x6- x7)* (-CONST4);
  *ezh=-( x1- x2+ x3- x4- x5- x6+ x7)* (-CONST4);
  x1= sppp2* xr1;
  x2= rv* spp2* xr2;
  x4= u2* t1* omr* f* xr2;
  x5= t3*(1.-3.* cppp2)* xr1;
  x6= t4*(1.-3.* cpp2)*(1.- u2*(1.+ rv)- u2* omr* f)* xr2;
  x7= u2* cpp2* omr*(1.-1./ rk2)*( f*( u2* t1- spp2-1./ rk2)+1./rk2)* xr2;
  *erh=( x1- x2- x4- x5+ x6+ x7)* (-CONST4);
  x1= xr1;
  x2= rh* xr2;
  x3=( rh+1.)* g* xr2;
  x4= t3* xr1;
  x5= t4*(1.- u2*(1.+ rv)- u2* omr* f)* xr2;
  x6=.5* u2* omr*( f*( u2* t1- spp2-1./ rk2)+1./ rk2)* xr2/ rk2;
  *eph=-( x1- x2+ x3- x4+ x5+ x6)* (-CONST4);

  return;
}

/*-----------------------------------------------------------------------*/

/* segment end contributions for thin wire approx. */
void gx( long double zz, long double rh, long double xk,
    complex long double *gz, complex long double *gzp)
{
  long double r, r2, rkz;

  r2= zz* zz+ rh* rh;
  r= sqrtl( r2);
  rkz= xk* r;
  *gz= cmplx( cosl( rkz),- sinl( rkz))/ r;
  *gzp=- cmplx(1.0, rkz)* *gz/ r2;

  return;
}

/*-----------------------------------------------------------------------*/

/* segment end contributions for ext. thin wire approx. */
void gxx( long double zz, long double rh, long double a, long double a2, long double xk, int ira,
    complex long double *g1, complex long double *g1p, complex long double *g2,
    complex long double *g2p, complex long double *g3, complex long double *gzp )
{
  long double r, r2, r4, rk, rk2, rh2, t1, t2;
  complex long double  gz, c1, c2, c3;

  r2= zz* zz+ rh* rh;
  r= sqrtl( r2);
  r4= r2* r2;
  rk= xk* r;
  rk2= rk* rk;
  rh2= rh* rh;
  t1=.25* a2* rh2/ r4;
  t2=.5* a2/ r2;
  c1= cmplx(1.0, rk);
  c2=3.* c1- rk2;
  c3= cmplx(6.0, rk)* rk2-15.* c1;
  gz= cmplx( cosl( rk),- sinl( rk))/ r;
  *g2= gz*(1.+ t1* c2);
  *g1= *g2- t2* c1* gz;
  gz= gz/ r2;
  *g2p= gz*( t1* c3- c1);
  *gzp= t2* c2* gz;
  *g3= *g2p+ *gzp;
  *g1p= *g3* zz;

  if( ira != 1)
  {
    *g3=( *g3+ *gzp)* rh;
    *gzp=- zz* c1* gz;

    if( rh <= 1.0e-10)
    {
      *g2=0.;
      *g2p=0.;
      return;
    }

    *g2= *g2/ rh;
    *g2p= *g2p* zz/ rh;
    return;

  } /* if( ira != 1) */

  t2=.5* a;
  *g2=- t2* c1* gz;
  *g2p= t2* gz* c2/ r2;
  *g3= rh2* *g2p- a* gz* c1;
  *g2p= *g2p* zz;
  *gzp=- zz* c1* gz;

  return;
}

/*-----------------------------------------------------------------------*/

/* hfk computes the h field of a uniform current */
/* filament by numerical integration */
void hfk( long double el1, long double el2, long double rhk,
    long double zpkx, long double *sgr, long double *sgi )
{
  int nx = 1, nma = 65536, nts = 4;
  int ns, nt;
  int flag = TRUE;
  long double rx = 1.0e-4;
  long double z, ze, s, ep, zend, dz=0., zp, dzot=0., t00r, g1r, g5r, t00i;
  long double g1i, g5i, t01r, g3r, t01i, g3i, t10r, t10i, te1i, te1r, t02r;
  long double g2r, g4r, t02i, g2i, g4i, t11r, t11i, t20r, t20i, te2i, te2r;

  zpka= zpkx;
  rhks= rhk* rhk;
  z= el1;
  ze= el2;
  s= ze- z;
  ep= s/(10.* nma);
  zend= ze- ep;
  *sgr=0.0;
  *sgi=0.0;
  ns= nx;
  nt=0;
  gh( z, &g1r, &g1i);

  while( TRUE )
  {
    if( flag )
    {
      dz= s/ ns;
      zp= z+ dz;

      if( zp > ze )
      {
	dz= ze- z;
	if( fabsl(dz) <= ep )
	{
	  *sgr= *sgr* rhk*.5;
	  *sgi= *sgi* rhk*.5;
	  return;
	}
      }

      dzot= dz*.5;
      zp= z+ dzot;
      gh( zp, &g3r, &g3i);
      zp= z+ dz;
      gh( zp, &g5r, &g5i);

    } /* if( flag ) */

    t00r=( g1r+ g5r)* dzot;
    t00i=( g1i+ g5i)* dzot;
    t01r=( t00r+ dz* g3r)*0.5;
    t01i=( t00i+ dz* g3i)*0.5;
    t10r=(4.0* t01r- t00r)/3.0;
    t10i=(4.0* t01i- t00i)/3.0;

    test( t01r, t10r, &te1r, t01i, t10i, &te1i, 0.);
    if( (te1i <= rx) && (te1r <= rx) )
    {
      *sgr= *sgr+ t10r;
      *sgi= *sgi+ t10i;
      nt += 2;

      z += dz;
      if( z >= zend)
      {
	*sgr= *sgr* rhk*.5;
	*sgi= *sgi* rhk*.5;
	return;
      }

      g1r= g5r;
      g1i= g5i;
      if( nt >= nts)
	if( ns > nx)
	{
	  ns= ns/2;
	  nt=1;
	}
      flag = TRUE;
      continue;

    } /* if( (te1i <= rx) && (te1r <= rx) ) */

    zp= z+ dz*0.25;
    gh( zp, &g2r, &g2i);
    zp= z+ dz*0.75;
    gh( zp, &g4r, &g4i);
    t02r=( t01r+ dzot*( g2r+ g4r))*0.5;
    t02i=( t01i+ dzot*( g2i+ g4i))*0.5;
    t11r=(4.0* t02r- t01r)/3.0;
    t11i=(4.0* t02i- t01i)/3.0;
    t20r=(16.0* t11r- t10r)/15.0;
    t20i=(16.0* t11i- t10i)/15.0;

    test( t11r, t20r, &te2r, t11i, t20i, &te2i, 0.);
    if( (te2i > rx) || (te2r > rx) )
    {
      nt=0;
      if( ns >= nma)
	fprintf( output_fp, "\n  STEP SIZE LIMITED AT Z= %10.5LF", z );
      else
      {
	ns= ns*2;
	dz= s/ ns;
	dzot= dz*0.5;
	g5r= g3r;
	g5i= g3i;
	g3r= g2r;
	g3i= g2i;

	flag = FALSE;
	continue;
      }

    } /* if( (te2i > rx) || (te2r > rx) ) */

    *sgr= *sgr+ t20r;
    *sgi= *sgi+ t20i;
    nt++;

    z += dz;
    if( z >= zend)
    {
      *sgr= *sgr* rhk*.5;
      *sgi= *sgi* rhk*.5;
      return;
    }

    g1r= g5r;
    g1i= g5i;
    if( nt >= nts)
      if( ns > nx)
      {
	ns= ns/2;
	nt=1;
      }
    flag = TRUE;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* hintg computes the h field of a patch current */
void hintg( long double xi, long double yi, long double zi )
{
  int ip;
  long double rx, ry, rfl, xymag, pxx, pyy, cth;
  long double rz, rsq, r, rk, cr, sr, t1zr, t2zr;
  complex long double  gam, f1x, f1y, f1z, f2x, f2y, f2z, rrv, rrh;

  rx= xi- xj;
  ry= yi- yj;
  rfl=-1.;
  exk=CPLX_00;
  eyk=CPLX_00;
  ezk=CPLX_00;
  exs=CPLX_00;
  eys=CPLX_00;
  ezs=CPLX_00;

  for( ip = 1; ip <= ksymp; ip++ )
  {
    rfl=- rfl;
    rz= zi- zj* rfl;
    rsq= rx* rx+ ry* ry+ rz* rz;

    if( rsq < 1.0e-20)
      continue;

    r = sqrtl( rsq );
    rk= TP* r;
    cr= cosl( rk);
    sr= sinl( rk);
    gam=-( cmplx(cr,-sr)+rk*cmplx(sr,cr) )/( FPI*rsq*r )* s;
    exc= gam* rx;
    eyc= gam* ry;
    ezc= gam* rz;
    t1zr= t1zj* rfl;
    t2zr= t2zj* rfl;
    f1x= eyc* t1zr- ezc* t1yj;
    f1y= ezc* t1xj- exc* t1zr;
    f1z= exc* t1yj- eyc* t1xj;
    f2x= eyc* t2zr- ezc* t2yj;
    f2y= ezc* t2xj- exc* t2zr;
    f2z= exc* t2yj- eyc* t2xj;

    if( ip != 1)
    {
      if( iperf == 1)
      {
	f1x=- f1x;
	f1y=- f1y;
	f1z=- f1z;
	f2x=- f2x;
	f2y=- f2y;
	f2z=- f2z;
      }
      else
      {
	xymag= sqrtl( rx* rx+ ry* ry);
	if( xymag <= 1.0e-6)
	{
	  pxx=0.;
	  pyy=0.;
	  cth=1.;
	  rrv=CPLX_10;
	}
	else
	{
	  pxx=- ry/ xymag;
	  pyy= rx/ xymag;
	  cth= rz/ r;
	  rrv= csqrtl(1.- zrati* zrati*(1.- cth* cth));

	} /* if( xymag <= 1.0e-6) */

	rrh= zrati* cth;
	rrh=( rrh- rrv)/( rrh+ rrv);
	rrv= zrati* rrv;
	rrv=-( cth- rrv)/( cth+ rrv);
	gam=( f1x* pxx+ f1y* pyy)*( rrv- rrh);
	f1x= f1x* rrh+ gam* pxx;
	f1y= f1y* rrh+ gam* pyy;
	f1z= f1z* rrh;
	gam=( f2x* pxx+ f2y* pyy)*( rrv- rrh);
	f2x= f2x* rrh+ gam* pxx;
	f2y= f2y* rrh+ gam* pyy;
	f2z= f2z* rrh;

      } /* if( iperf == 1) */

    } /* if( ip != 1) */

    exk += f1x;
    eyk += f1y;
    ezk += f1z;
    exs += f2x;
    eys += f2y;
    ezs += f2z;

  } /* for( ip = 1; ip <= ksymp; ip++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* hsfld computes the h field for constant, sine, and */
/* cosine current on a segment including ground effects. */
void hsfld( long double xi, long double yi, long double zi, long double ai )
{
  int ip;
  long double xij, yij, rfl, salpr, zij, zp, rhox, rhoy, rhoz, rh, phx;
  long double phy, phz, rmag, xymag, xspec, yspec, rhospc, px, py, cth;
  complex long double hpk, hps, hpc, qx, qy, qz, rrv, rrh, zratx;

  xij= xi- xj;
  yij= yi- yj;
  rfl=-1.;

  for( ip = 0; ip < ksymp; ip++ )
  {
    rfl=- rfl;
    salpr= salpj* rfl;
    zij= zi- rfl* zj;
    zp= xij* cabj+ yij* sabj+ zij* salpr;
    rhox= xij- cabj* zp;
    rhoy= yij- sabj* zp;
    rhoz= zij- salpr* zp;
    rh= sqrtl( rhox* rhox+ rhoy* rhoy+ rhoz* rhoz+ ai* ai);

    if( rh <= 1.0e-10)
    {
      exk=0.;
      eyk=0.;
      ezk=0.;
      exs=0.;
      eys=0.;
      ezs=0.;
      exc=0.;
      eyc=0.;
      ezc=0.;
      continue;
    }

    rhox= rhox/ rh;
    rhoy= rhoy/ rh;
    rhoz= rhoz/ rh;
    phx= sabj* rhoz- salpr* rhoy;
    phy= salpr* rhox- cabj* rhoz;
    phz= cabj* rhoy- sabj* rhox;

    hsflx( s, rh, zp, &hpk, &hps, &hpc);

    if( ip == 1 )
    {
      if( iperf != 1 )
      {
	zratx= zrati;
	rmag= sqrtl( zp* zp+ rh* rh);
	xymag= sqrtl( xij* xij+ yij* yij);

	/* set parameters for radial wire ground screen. */
	if( nradl != 0)
	{
	  xspec=( xi* zj+ zi* xj)/( zi+ zj);
	  yspec=( yi* zj+ zi* yj)/( zi+ zj);
	  rhospc= sqrtl( xspec* xspec+ yspec* yspec+ t2* t2);

	  if( rhospc <= scrwl)
	  {
	    rrv= t1* rhospc* logl( rhospc/ t2);
	    zratx=( rrv* zrati)/( ETA* zrati+ rrv);
	  }
	}

	/* calculation of reflection coefficients when ground is specified. */
	if( xymag <= 1.0e-6)
	{
	  px=0.;
	  py=0.;
	  cth=1.;
	  rrv=CPLX_10;
	}
	else
	{
	  px=- yij/ xymag;
	  py= xij/ xymag;
	  cth= zij/ rmag;
	  rrv= csqrtl(1.- zratx* zratx*(1.- cth* cth));
	}

	rrh= zratx* cth;
	rrh=-( rrh- rrv)/( rrh+ rrv);
	rrv= zratx* rrv;
	rrv=( cth- rrv)/( cth+ rrv);
	qy=( phx* px+ phy* py)*( rrv- rrh);
	qx= qy* px+ phx* rrh;
	qy= qy* py+ phy* rrh;
	qz= phz* rrh;
	exk= exk- hpk* qx;
	eyk= eyk- hpk* qy;
	ezk= ezk- hpk* qz;
	exs= exs- hps* qx;
	eys= eys- hps* qy;
	ezs= ezs- hps* qz;
	exc= exc- hpc* qx;
	eyc= eyc- hpc* qy;
	ezc= ezc- hpc* qz;
	continue;

      } /* if( iperf != 1 ) */

      exk= exk- hpk* phx;
      eyk= eyk- hpk* phy;
      ezk= ezk- hpk* phz;
      exs= exs- hps* phx;
      eys= eys- hps* phy;
      ezs= ezs- hps* phz;
      exc= exc- hpc* phx;
      eyc= eyc- hpc* phy;
      ezc= ezc- hpc* phz;
      continue;

    } /* if( ip == 1 ) */

    exk= hpk* phx;
    eyk= hpk* phy;
    ezk= hpk* phz;
    exs= hps* phx;
    eys= hps* phy;
    ezs= hps* phz;
    exc= hpc* phx;
    eyc= hpc* phy;
    ezc= hpc* phz;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* calculates h field of sine cosine, and constant current of segment */
void hsflx( long double s, long double rh, long double zpx,
    complex long double *hpk, complex long double *hps,
    complex long double *hpc )
{
  long double r1, r2, zp, z2a, hss, dh, z1;
  long double rhz, dk, cdk, sdk, hkr, hki, rh2;
  complex long double fjk, ekr1, ekr2, t1, t2, cons;

  fjk = -TPJ;
  if( rh >= 1.0e-10)
  {
    if( zpx >= 0.)
    {
      zp= zpx;
      hss=1.;
    }
    else
    {
      zp=- zpx;
      hss=-1.;
    }

    dh=.5* s;
    z1= zp+ dh;
    z2a= zp- dh;
    if( z2a >= 1.0e-7)
      rhz= rh/ z2a;
    else
      rhz=1.;

    dk= TP* dh;
    cdk= cosl( dk);
    sdk= sinl( dk);
    hfk(- dk, dk, rh* TP, zp* TP, &hkr, &hki);
    *hpk= cmplx( hkr, hki);

    if( rhz >= 1.0e-3)
    {
      rh2= rh* rh;
      r1= sqrtl( rh2+ z1* z1);
      r2= sqrtl( rh2+ z2a* z2a);
      ekr1= cexp( fjk* r1);
      ekr2= cexp( fjk* r2);
      t1= z1* ekr1/ r1;
      t2= z2a* ekr2/ r2;
      *hps=( cdk*( ekr2- ekr1)- CPLX_01* sdk*( t2+ t1))* hss;
      *hpc=- sdk*( ekr2+ ekr1)- CPLX_01* cdk*( t2- t1);
      cons=- CPLX_01/(2.* TP* rh);
      *hps= cons* *hps;
      *hpc= cons* *hpc;
      return;

    } /* if( rhz >= 1.0e-3) */

    ekr1= cmplx( cdk, sdk)/( z2a* z2a);
    ekr2= cmplx( cdk,- sdk)/( z1* z1);
    t1= TP*(1./ z1-1./ z2a);
    t2= cexp( fjk* zp)* rh/ PI8;
    *hps= t2*( t1+( ekr1+ ekr2)* sdk)* hss;
    *hpc= t2*(- CPLX_01* t1+( ekr1- ekr2)* cdk);
    return;

  } /* if( rh >= 1.0e-10) */

  *hps=CPLX_00;
  *hpc=CPLX_00;
  *hpk=CPLX_00;

  return;
}

/*-----------------------------------------------------------------------*/

/* nefld computes the near field at specified points in space after */
/* the structure currents have been computed. */
void nefld( long double xob, long double yob, long double zob,
    complex long double *ex, complex long double *ey, complex long double *ez )
{
  int i, ix, ipr, iprx, jc, ipa;
  long double zp, xi, ax;
  complex long double acx, bcx, ccx;

  *ex=CPLX_00;
  *ey=CPLX_00;
  *ez=CPLX_00;
  ax=0.;

  if( n != 0)
  {
    for( i = 0; i < n; i++ )
    {
      xj= xob- x[i];
      yj= yob- y[i];
      zj= zob- z[i];
      zp= cab[i]* xj+ sab[i]* yj+ salp[i]* zj;

      if( fabsl( zp) > 0.5001* si[i])
	continue;

      zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
      xj= bi[i];

      if( zp > 0.9* xj* xj)
	continue;

      ax= xj;
      break;

    } /* for( i = 0; i < n; i++ ) */

    for( i = 0; i < n; i++ )
    {
      ix = i+1;
      s= si[i];
      b= bi[i];
      xj= x[i];
      yj= y[i];
      zj= z[i];
      cabj= cab[i];
      sabj= sab[i];
      salpj= salp[i];

      if( iexk != 0)
      {
	ipr= icon1[i];

	if( ipr < 0 )
	{
	  ipr = -ipr;
	  iprx = ipr-1;

	  if( -icon1[iprx] != ix )
	    ind1=2;
	  else
	  {
	    xi= fabsl( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabsl(bi[iprx]/b-1.) > 1.0e-6) )
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* if( ipr < 0 ) */
	else
	  if( ipr == 0 )
	    ind1=1;
	  else
	  {
	    iprx = ipr-1;

	    if( ipr != ix )
	    {
	      if( icon2[iprx] != ix )
		ind1=2;
	      else
	      {
		xi= fabsl( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
		if( (xi < 0.999999) || (fabsl(bi[iprx]/b-1.) > 1.0e-6) )
		  ind1=2;
		else
		  ind1=0;
	      }
	    } /* if( ipr != ix ) */
	    else
	    {
	      if( cabj* cabj+ sabj* sabj > 1.0e-8)
		ind1=2;
	      else
		ind1=0;
	    }
	  } /* else */

	ipr= icon2[i];

	if( ipr < 0 )
	{
	  ipr = -ipr;
	  iprx = ipr-1;

	  if( -icon2[iprx] != ix )
	    ind1=2;
	  else
	  {
	    xi= fabsl( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
	    if( (xi < 0.999999) || (fabsl(bi[iprx]/b-1.) > 1.0e-6) )
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* if( ipr < 0 ) */
	else
	  if( ipr == 0 )
	    ind2=1;
	  else
	  {
	    iprx = ipr-1;

	    if( ipr != ix )
	    {
	      if( icon1[iprx] != ix )
		ind2=2;
	      else
	      {
		xi= fabsl( cabj* cab[iprx]+ sabj* sab[iprx]+ salpj* salp[iprx]);
		if( (xi < 0.999999) || (fabsl(bi[iprx]/b-1.) > 1.0e-6) )
		  ind2=2;
		else
		  ind2=0;
	      }
	    } /* if( ipr != (i+1) ) */
	    else
	    {
	      if( cabj* cabj+ sabj* sabj > 1.0e-8)
		ind1=2;
	      else
		ind1=0;
	    }

	  } /* else */

      } /* if( iexk != 0) */

      efld( xob, yob, zob, ax,1);
      acx= cmplx( air[i], aii[i]);
      bcx= cmplx( bir[i], bii[i]);
      ccx= cmplx( cir[i], cii[i]);
      *ex += exk* acx+ exs* bcx+ exc* ccx;
      *ey += eyk* acx+ eys* bcx+ eyc* ccx;
      *ez += ezk* acx+ ezs* bcx+ ezc* ccx;

    } /* for( i = 0; i < n; i++ ) */

    if( m == 0)
      return;

  } /* if( n != 0) */

  jc= n-1;
  for( i = 0; i < m; i++ )
  {
    s= pbi[i];
    xj= px[i];
    yj= py[i];
    zj= pz[i];
    t1xj= t1x[i];
    t1yj= t1y[i];
    t1zj= t1z[i];
    t2xj= t2x[i];
    t2yj= t2y[i];
    t2zj= t2z[i];
    jc += 3;
    acx= t1xj* cur[jc-2]+ t1yj* cur[jc-1]+ t1zj* cur[jc];
    bcx= t2xj* cur[jc-2]+ t2yj* cur[jc-1]+ t2zj* cur[jc];

    for( ipa = 0; ipa < ksymp; ipa++ )
    {
      ipgnd= ipa+1;
      unere( xob, yob, zob);
      *ex= *ex+ acx* exk+ bcx* exs;
      *ey= *ey+ acx* eyk+ bcx* eys;
      *ez= *ez+ acx* ezk+ bcx* ezs;
    }

  } /* for( i = 0; i < m; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* compute near e or h fields over a range of points */
void nfpat( void )
{
  int i, j, kk;
  long double znrt, cth=0., sth=0., ynrt, cph=0., sph=0., xnrt, xob, yob;
  long double zob, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, xxx;
  complex long double ex, ey, ez;

  if( nfeh != 1)
  {
    fprintf( output_fp,	"\n\n\n"
	"                             "
	"-------- NEAR ELECTRIC FIELDS --------\n"
	"     ------- LOCATION -------     ------- EX ------    ------- EY ------    ------- EZ ------\n"
	"      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
	"    METERS    METERS    METERS     VOLTS/M  DEGREES    VOLTS/M   DEGREES     VOLTS/M  DEGREES" );
  }
  else
  {
    fprintf( output_fp,	"\n\n\n"
	"                                   "
	"-------- NEAR MAGNETIC FIELDS ---------\n\n"
	"     ------- LOCATION -------     ------- HX ------    ------- HY ------    ------- HZ ------\n"
	"      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n"
	"    METERS    METERS    METERS      AMPS/M  DEGREES      AMPS/M  DEGREES      AMPS/M  DEGREES" );
  }

  znrt= znr- dznr;
  for( i = 0; i < nrz; i++ )
  {
    znrt += dznr;
    if( near != 0)
    {
      cth= cosl( TA* znrt);
      sth= sinl( TA* znrt);
    }

    ynrt= ynr- dynr;
    for( j = 0; j < nry; j++ )
    {
      ynrt += dynr;
      if( near != 0)
      {
	cph= cosl( TA* ynrt);
	sph= sinl( TA* ynrt);
      }

      xnrt= xnr- dxnr;
      for( kk = 0; kk < nrx; kk++ )
      {
	xnrt += dxnr;
	if( near != 0)
	{
	  xob= xnrt* sth* cph;
	  yob= xnrt* sth* sph;
	  zob= xnrt* cth;
	}
	else
	{
	  xob= xnrt;
	  yob= ynrt;
	  zob= znrt;
	}

	tmp1= xob/ wlam;
	tmp2= yob/ wlam;
	tmp3= zob/ wlam;

	if( nfeh != 1)
	  nefld( tmp1, tmp2, tmp3, &ex, &ey, &ez);
	else
	  nhfld( tmp1, tmp2, tmp3, &ex, &ey, &ez);

	tmp1= cabsl( ex);
	tmp2= cang( ex);
	tmp3= cabsl( ey);
	tmp4= cang( ey);
	tmp5= cabsl( ez);
	tmp6= cang( ez);

	fprintf( output_fp, "\n"
	    " %9.4LF %9.4LF %9.4LF  %11.4LE %7.2LF  %11.4LE %7.2LF  %11.4LE %7.2LF",
	    xob, yob, zob, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );

	if( iplp1 != 2)
	  continue;

	if( iplp4 < 0 )
	  xxx= xob;
	else
	  if( iplp4 == 0 )
	    xxx= yob;
	  else
	    xxx= zob;

	if( iplp2 == 2)
	{
	  switch( iplp3 )
	  {
	    case 1:
	      fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, tmp1, tmp2 );
	      break;
	    case 2:
	      fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, tmp3, tmp4 );
	      break;
	    case 3:
	      fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, tmp5, tmp6 );
	      break;
	    case 4:
	      fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE %12.4LE %12.4LE %12.4LE %12.4LE\n",
		  xxx, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
	  }
	  continue;
	}

	if( iplp2 != 1)
	  continue;

	switch( iplp3 )
	{
	  case 1:
	    fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, creall(ex), cimagl(ex) );
	    break;
	  case 2:
	    fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, creall(ey), cimagl(ey) );
	    break;
	  case 3:
	    fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", xxx, creall(ez), cimagl(ez) );
	    break;
	  case 4:
	    fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE %12.4LE %12.4LE %12.4LE %12.4LE\n",
		xxx,creall(ex),cimagl(ex),creall(ey),cimagl(ey),creall(ez),cimagl(ez) );
	}
      } /* for( kk = 0; kk < nrx; kk++ ) */

    } /* for( j = 0; j < nry; j++ ) */

  } /* for( i = 0; i < nrz; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* nhfld computes the near field at specified points in space after */
/* the structure currents have been computed. */

void nhfld( long double xob, long double yob, long double zob,
    complex long double *hx, complex long double *hy, complex long double *hz )
{
  int i, jc;
  long double ax, zp;
  complex long double acx, bcx, ccx;

  *hx=CPLX_00;
  *hy=CPLX_00;
  *hz=CPLX_00;
  ax=0.;

  if( n != 0)
  {
    for( i = 0; i < n; i++ )
    {
      xj= xob- x[i];
      yj= yob- y[i];
      zj= zob- z[i];
      zp= cab[i]* xj+ sab[i]* yj+ salp[i]* zj;

      if( fabsl( zp) > 0.5001* si[i])
	continue;

      zp= xj* xj+ yj* yj+ zj* zj- zp* zp;
      xj= bi[i];

      if( zp > 0.9* xj* xj)
	continue;

      ax= xj;
      break;
    }

    for( i = 0; i < n; i++ )
    {
      s= si[i];
      b= bi[i];
      xj= x[i];
      yj= y[i];
      zj= z[i];
      cabj= cab[i];
      sabj= sab[i];
      salpj= salp[i];
      hsfld( xob, yob, zob, ax);
      acx= cmplx( air[i], aii[i]);
      bcx= cmplx( bir[i], bii[i]);
      ccx= cmplx( cir[i], cii[i]);
      *hx += exk* acx+ exs* bcx+ exc* ccx;
      *hy += eyk* acx+ eys* bcx+ eyc* ccx;
      *hz += ezk* acx+ ezs* bcx+ ezc* ccx;
    }

    if( m == 0)
      return;

  } /* if( n != 0) */

  jc= n-1;
  for( i = 0; i < m; i++ )
  {
    s= pbi[i];
    xj= px[i];
    yj= py[i];
    zj= pz[i];
    t1xj= t1x[i];
    t1yj= t1y[i];
    t1zj= t1z[i];
    t2xj= t2x[i];
    t2yj= t2y[i];
    t2zj= t2z[i];
    hintg( xob, yob, zob);
    jc += 3;
    acx= t1xj* cur[jc-2]+ t1yj* cur[jc-1]+ t1zj* cur[jc];
    bcx= t2xj* cur[jc-2]+ t2yj* cur[jc-1]+ t2zj* cur[jc];
    *hx= *hx+ acx* exk+ bcx* exs;
    *hy= *hy+ acx* eyk+ bcx* eys;
    *hz= *hz+ acx* ezk+ bcx* ezs;
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* integrate over patches at wire connection point */
void pcint( long double xi, long double yi, long double zi, long double cabi,
    long double sabi, long double salpi, complex long double *e )
{
  int nint, i1, i2;
  long double d, ds, da, gcon, fcon, xxj, xyj, xzj, xs, s1;
  long double xss, yss, zss, s2x, s2, g1, g2, g3, g4, f2, f1;
  complex long double e1, e2, e3, e4, e5, e6, e7, e8, e9;

  nint = 10;
  d= sqrtl( s)*.5;
  ds=4.* d/ (long double) nint;
  da= ds* ds;
  gcon=1./ s;
  fcon=1./(2.* TP* d);
  xxj= xj;
  xyj= yj;
  xzj= zj;
  xs= s;
  s= da;
  s1= d+ ds*.5;
  xss= xj+ s1*( t1xj+ t2xj);
  yss= yj+ s1*( t1yj+ t2yj);
  zss= zj+ s1*( t1zj+ t2zj);
  s1= s1+ d;
  s2x= s1;
  e1=CPLX_00;
  e2=CPLX_00;
  e3=CPLX_00;
  e4=CPLX_00;
  e5=CPLX_00;
  e6=CPLX_00;
  e7=CPLX_00;
  e8=CPLX_00;
  e9=CPLX_00;

  for( i1 = 0; i1 < nint; i1++ )
  {
    s1= s1- ds;
    s2= s2x;
    xss= xss- ds* t1xj;
    yss= yss- ds* t1yj;
    zss= zss- ds* t1zj;
    xj= xss;
    yj= yss;
    zj= zss;

    for( i2 = 0; i2 < nint; i2++ )
    {
      s2= s2- ds;
      xj= xj- ds* t2xj;
      yj= yj- ds* t2yj;
      zj= zj- ds* t2zj;
      unere( xi, yi, zi);
      exk= exk* cabi+ eyk* sabi+ ezk* salpi;
      exs= exs* cabi+ eys* sabi+ ezs* salpi;
      g1=( d+ s1)*( d+ s2)* gcon;
      g2=( d- s1)*( d+ s2)* gcon;
      g3=( d- s1)*( d- s2)* gcon;
      g4=( d+ s1)*( d- s2)* gcon;
      f2=( s1* s1+ s2* s2)* TP;
      f1= s1/ f2-( g1- g2- g3+ g4)* fcon;
      f2= s2/ f2-( g1+ g2- g3- g4)* fcon;
      e1= e1+ exk* g1;
      e2= e2+ exk* g2;
      e3= e3+ exk* g3;
      e4= e4+ exk* g4;
      e5= e5+ exs* g1;
      e6= e6+ exs* g2;
      e7= e7+ exs* g3;
      e8= e8+ exs* g4;
      e9= e9+ exk* f1+ exs* f2;

    } /* for( i2 = 0; i2 < nint; i2++ ) */

  } /* for( i1 = 0; i1 < nint; i1++ ) */

  e[0]= e1;
  e[1]= e2;
  e[2]= e3;
  e[3]= e4;
  e[4]= e5;
  e[5]= e6;
  e[6]= e7;
  e[7]= e8;
  e[8]= e9;
  xj= xxj;
  yj= xyj;
  zj= xzj;
  s= xs;

  return;
}

/*-----------------------------------------------------------------------*/

/* calculates the electric field due to unit current */
/* in the t1 and t2 directions on a patch */
void unere( long double xob, long double yob, long double zob )
{
  long double zr, t1zr, t2zr, rx, ry, rz, r, tt1;
  long double tt2, rt, xymag, px, py, cth, r2;
  complex long double er, q1, q2, rrv, rrh, edp;

  zr= zj;
  t1zr= t1zj;
  t2zr= t2zj;

  if( ipgnd == 2)
  {
    zr=- zr;
    t1zr=- t1zr;
    t2zr=- t2zr;
  }

  rx= xob- xj;
  ry= yob- yj;
  rz= zob- zr;
  r2= rx* rx+ ry* ry+ rz* rz;

  if( r2 <= 1.0e-20)
  {
    exk=CPLX_00;
    eyk=CPLX_00;
    ezk=CPLX_00;
    exs=CPLX_00;
    eys=CPLX_00;
    ezs=CPLX_00;
    return;
  }

  r= sqrtl( r2);
  tt1=- TP* r;
  tt2= tt1* tt1;
  rt= r2* r;
  er= cmplx( sinl( tt1),- cosl( tt1))*( CONST2* s);
  q1= cmplx( tt2-1., tt1)* er/ rt;
  q2= cmplx(3.- tt2,-3.* tt1)* er/( rt* r2);
  er = q2*( t1xj* rx+ t1yj* ry+ t1zr* rz);
  exk= q1* t1xj+ er* rx;
  eyk= q1* t1yj+ er* ry;
  ezk= q1* t1zr+ er* rz;
  er= q2*( t2xj* rx+ t2yj* ry+ t2zr* rz);
  exs= q1* t2xj+ er* rx;
  eys= q1* t2yj+ er* ry;
  ezs= q1* t2zr+ er* rz;

  if( ipgnd == 1)
    return;

  if( iperf == 1)
  {
    exk=- exk;
    eyk=- eyk;
    ezk=- ezk;
    exs=- exs;
    eys=- eys;
    ezs=- ezs;
    return;
  }

  xymag= sqrtl( rx* rx+ ry* ry);
  if( xymag <= 1.0e-6)
  {
    px=0.;
    py=0.;
    cth=1.;
    rrv=CPLX_10;
  }
  else
  {
    px=- ry/ xymag;
    py= rx/ xymag;
    cth= rz/ sqrtl( xymag* xymag+ rz* rz);
    rrv= csqrtl(1.- zrati* zrati*(1.- cth* cth));
  }

  rrh= zrati* cth;
  rrh=( rrh- rrv)/( rrh+ rrv);
  rrv= zrati* rrv;
  rrv=-( cth- rrv)/( cth+ rrv);
  edp=( exk* px+ eyk* py)*( rrh- rrv);
  exk= exk* rrv+ edp* px;
  eyk= eyk* rrv+ edp* py;
  ezk= ezk* rrv;
  edp=( exs* px+ eys* py)*( rrh- rrv);
  exs= exs* rrv+ edp* px;
  eys= eys* rrv+ edp* py;
  ezs= ezs* rrv;

  return;
}

/*-----------------------------------------------------------------------*/


