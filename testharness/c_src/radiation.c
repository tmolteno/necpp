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

/* common  /data/ */
extern int n, np, m, mp, ipsym, npm, np2m, np3m; /* n+m,n+2m,n+3m */
extern int *icon1, *icon2, *itag;
extern long double *x, *y, *z, *si, *bi;
extern long double *x2, *y2, *z2, *cab, *sab, *salp;
extern long double *t1x, *t1y, *t1z, *t2x, *t2y, *t2z;
extern long double *px, *py, *pz, *pbi, *psalp;
extern long double wlam;

/* common  /gnd/ */
extern int ksymp, ifar, iperf, nradl;
extern long double t2, cl, ch, scrwl, scrwr;
extern complex long double zrati, zrati2, t1, frati;

/* common  /crnt/ */
extern long double *air, *aii, *bir, *bii, *cir, *cii;
extern complex long double *cur;

/* common  /gwav/ */
extern long double r1, r2, zmh, zph;
extern complex long double u, u2, xx1, xx2;

/* common  /fpat/ */
extern int near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp, inor, iax, ixtyp;
extern long double thets, phis, dth, dph, rfld, gnor, clt, cht, epsr2, sig2;
extern long double xpr6, pinr, pnlr, ploss, xnr, ynr, znr, dxnr, dynr, dznr;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/* common  /save/ */
extern int *ip;
extern long double epsr, sig, scrwlt, scrwrt, fmhz;

/* common  /plot/ */
extern int iplp1, iplp2, iplp3, iplp4;


/*-----------------------------------------------------------------------*/

/* ffld calculates the far zone radiated electric fields, */
/* the factor exp(j*k*r)/(r/lamda) not included */
void ffld( long double thet, long double phi,
    complex long double *eth, complex long double *eph )
{
  int k, i, ip, jump;
  long double phx, phy, roz, rozs, thx, thy, thz, rox, roy;
  long double tthet=0., darg=0., omega, el, sill, top, bot, a;
  long double too, boo, b, c, d, rr, ri, arg, dr, rfl, rrz;
  complex long double cix, ciy, ciz, exa, ccx, ccy, ccz, cdp;
  complex long double zrsin, rrv, rrh, rrv1, rrh1, rrv2, rrh2;
  complex long double tix, tiy, tiz, zscrn, ex, ey, ez, gx, gy, gz;

  phx=- sinl( phi);
  phy= cosl( phi);
  roz= cosl( thet);
  rozs= roz;
  thx= roz* phy;
  thy=- roz* phx;
  thz=- sinl( thet);
  rox=- thz* phy;
  roy= thz* phx;

  jump = FALSE;
  if( n != 0)
  {
    /* loop for structure image if any */
    /* calculation of reflection coeffecients */
    for( k = 0; k < ksymp; k++ )
    {
      if( k != 0 )
      {
	/* for perfect ground */
	if( iperf == 1)
	{
	  rrv=-CPLX_10;
	  rrh=-CPLX_10;
	}
	else
	{
	  /* for infinite planar ground */
	  zrsin= csqrtl(1.- zrati* zrati* thz* thz);
	  rrv=-( roz- zrati* zrsin)/( roz+ zrati* zrsin);
	  rrh=( zrati* roz- zrsin)/( zrati* roz+ zrsin);

	} /* if( iperf == 1) */

	/* for the cliff problem, two reflction coefficients calculated */
	if( ifar > 1)
	{
	  rrv1= rrv;
	  rrh1= rrh;
	  tthet= tanl( thet);

	  if( ifar != 4)
	  {
	    zrsin= csqrtl(1.- zrati2* zrati2* thz* thz);
	    rrv2=-( roz- zrati2* zrsin)/( roz+ zrati2* zrsin);
	    rrh2=( zrati2* roz- zrsin)/( zrati2* roz+ zrsin);
	    darg=- TP*2.* ch* roz;
	  }
	} /* if( ifar > 1) */

	roz=- roz;
	ccx= cix;
	ccy= ciy;
	ccz= ciz;

      } /* if( k != 0 ) */

      cix=CPLX_00;
      ciy=CPLX_00;
      ciz=CPLX_00;

      /* loop over structure segments */
      for( i = 0; i < n; i++ )
      {
	omega=-( rox* cab[i]+ roy* sab[i]+ roz* salp[i]);
	el= PI* si[i];
	sill= omega* el;
	top= el+ sill;
	bot= el- sill;

	if( fabsl( omega) >= 1.0e-7)
	  a=2.* sinl( sill)/ omega;
	else
	  a=(2.- omega* omega* el* el/3.)* el;

	if( fabsl( top) >= 1.0e-7)
	  too= sinl( top)/ top;
	else
	  too=1.- top* top/6.;

	if( fabsl( bot) >= 1.0e-7)
	  boo= sinl( bot)/ bot;
	else
	  boo=1.- bot* bot/6.;

	b= el*( boo- too);
	c= el*( boo+ too);
	rr= a* air[i]+ b* bii[i]+ c* cir[i];
	ri= a* aii[i]- b* bir[i]+ c* cii[i];
	arg= TP*( x[i]* rox+ y[i]* roy+ z[i]* roz);

	if( (k != 1) || (ifar < 2) )
	{
	  /* summation for far field integral */
	  exa= cmplx( cosl( arg), sinl( arg))* cmplx( rr, ri);
	  cix= cix+ exa* cab[i];
	  ciy= ciy+ exa* sab[i];
	  ciz= ciz+ exa* salp[i];
	  continue;
	}

	/* calculation of image contribution */
	/* in cliff and ground screen problems */

	/* specular point distance */
	dr= z[i]* tthet;

	d= dr* phy+ x[i];
	if( ifar == 2)
	{
	  if(( cl- d) > 0.)
	  {
	    rrv= rrv1;
	    rrh= rrh1;
	  }
	  else
	  {
	    rrv= rrv2;
	    rrh= rrh2;
	    arg= arg+ darg;
	  }
	} /* if( ifar == 2) */
	else
	{
	  d= sqrtl( d*d + (y[i]-dr*phx)*(y[i]-dr*phx) );
	  if( ifar == 3)
	  {
	    if(( cl- d) > 0.)
	    {
	      rrv= rrv1;
	      rrh= rrh1;
	    }
	    else
	    {
	      rrv= rrv2;
	      rrh= rrh2;
	      arg= arg+ darg;
	    }
	  } /* if( ifar == 3) */
	  else
	  {
	    if(( scrwl- d) >= 0.)
	    {
	      /* radial wire ground screen reflection coefficient */
	      d= d+ t2;
	      zscrn= t1* d* logl( d/ t2);
	      zscrn=( zscrn* zrati)/( ETA* zrati+ zscrn);
	      zrsin= csqrtl(1.- zscrn* zscrn* thz* thz);
	      rrv=( roz+ zscrn* zrsin)/(- roz+ zscrn* zrsin);
	      rrh=( zscrn* roz+ zrsin)/( zscrn* roz- zrsin);
	    } /* if(( scrwl- d) < 0.) */
	    else
	    {
	      if( ifar == 4)
	      {
		rrv= rrv1;
		rrh= rrh1;
	      } /* if( ifar == 4) */
	      else
	      {
		if( ifar == 5)
		  d= dr* phy+ x[i];

		if(( cl- d) > 0.)
		{
		  rrv= rrv1;
		  rrh= rrh1;
		}
		else
		{
		  rrv= rrv2;
		  rrh= rrh2;
		  arg= arg+ darg;
		} /* if(( cl- d) > 0.) */

	      } /* if( ifar == 4) */

	    } /* if(( scrwl- d) < 0.) */

	  } /* if( ifar == 3) */

	} /* if( ifar == 2) */

	/* contribution of each image segment modified by */
	/* reflection coef, for cliff and ground screen problems */
	exa= cmplx( cosl( arg), sinl( arg))* cmplx( rr, ri);
	tix= exa* cab[i];
	tiy= exa* sab[i];
	tiz= exa* salp[i];
	cdp=( tix* phx+ tiy* phy)*( rrh- rrv);
	cix= cix+ tix* rrv+ cdp* phx;
	ciy= ciy+ tiy* rrv+ cdp* phy;
	ciz= ciz- tiz* rrv;

      } /* for( i = 0; i < n; i++ ) */

      if( k == 0 )
	continue;

      /* calculation of contribution of structure image for infinite ground */
      if( ifar < 2)
      {
	cdp=( cix* phx+ ciy* phy)*( rrh- rrv);
	cix= ccx+ cix* rrv+ cdp* phx;
	ciy= ccy+ ciy* rrv+ cdp* phy;
	ciz= ccz- ciz* rrv;
      }
      else
      {
	cix= cix+ ccx;
	ciy= ciy+ ccy;
	ciz= ciz+ ccz;
      }

    } /* for( k=0; k < ksymp; k++ ) */

    if( m > 0)
      jump = TRUE;
    else
    {
      *eth=( cix* thx+ ciy* thy+ ciz* thz)* CONST3;
      *eph=( cix* phx+ ciy* phy)* CONST3;
      return;
    }

  } /* if( n != 0) */

  if( ! jump )
  {
    cix=CPLX_00;
    ciy=CPLX_00;
    ciz=CPLX_00;
  }

  /* electric field components */
  roz= rozs;
  rfl=-1.;
  for( ip = 0; ip < ksymp; ip++ )
  {
    rfl=- rfl;
    rrz= roz* rfl;
    fflds( rox, roy, rrz, &cur[n], &gx, &gy, &gz);

    if( ip != 1 )
    {
      ex= gx;
      ey= gy;
      ez= gz;
      continue;
    }

    if( iperf == 1)
    {
      gx=- gx;
      gy=- gy;
      gz=- gz;
    }
    else
    {
      rrv= csqrtl(1.- zrati* zrati* thz* thz);
      rrh= zrati* roz;
      rrh=( rrh- rrv)/( rrh+ rrv);
      rrv= zrati* rrv;
      rrv=-( roz- rrv)/( roz+ rrv);
      *eth=( gx* phx+ gy* phy)*( rrh- rrv);
      gx= gx* rrv+ *eth* phx;
      gy= gy* rrv+ *eth* phy;
      gz= gz* rrv;

    } /* if( iperf == 1) */

    ex= ex+ gx;
    ey= ey+ gy;
    ez= ez- gz;

  } /* for( ip = 0; ip < ksymp; ip++ ) */

  ex= ex+ cix* CONST3;
  ey= ey+ ciy* CONST3;
  ez= ez+ ciz* CONST3;
  *eth= ex* thx+ ey* thy+ ez* thz;
  *eph= ex* phx+ ey* phy;

  return;
}

/*-----------------------------------------------------------------------*/

/* calculates the xyz components of the electric */
/* field due to surface currents */
void fflds( long double rox, long double roy, long double roz,
    complex long double *scur, complex long double *ex,
    complex long double *ey, complex long double *ez )
{
  long double *xs, *ys, *zs, *s;
  int j, i, k;
  long double arg;
  complex long double ct;

  xs = px;
  ys = py;
  zs = pz;
  s = pbi;

  *ex=CPLX_00;
  *ey=CPLX_00;
  *ez=CPLX_00;

  i= -1;
  for( j = 0; j < m; j++ )
  {
    i++;
    arg= TP*( rox* xs[i]+ roy* ys[i]+ roz* zs[i]);
    ct= cmplx( cosl( arg)* s[i], sinl( arg)* s[i]);
    k=3*j;
    *ex += scur[k  ]* ct;
    *ey += scur[k+1]* ct;
    *ez += scur[k+2]* ct;
  }

  ct= rox* *ex+ roy* *ey+ roz* *ez;
  *ex= CONST4*( ct* rox- *ex);
  *ey= CONST4*( ct* roy- *ey);
  *ez= CONST4*( ct* roz- *ez);

  return;
}

/*-----------------------------------------------------------------------*/

/* gfld computes the radiated field including ground wave. */
void gfld( long double rho, long double phi, long double rz,
    complex long double *eth, complex long double *epi,
    complex long double *erd, complex long double ux, int ksymp )
{
  int i, k;
  long double b, r, thet, arg, phx, phy, rx, ry, dx, dy, dz, rix, riy, rhs, rhp;
  long double rhx, rhy, calp, cbet, sbet, cph, sph, el, rfl, riz, thx, thy, thz;
  long double rxyz, rnx, rny, rnz, omega, sill, top, bot, a, too, boo, c, rr, ri;
  complex long double cix, ciy, ciz, exa, erv;
  complex long double ezv, erh, eph, ezh, ex, ey;

  r= sqrtl( rho*rho+ rz*rz );
  if( (ksymp == 1) || (cabs(ux) > .5) || (r > 1.e5) )
  {
    /* computation of space wave only */
    if( rz >= 1.0e-20)
      thet= atanl( rho/ rz);
    else
      thet= PI*.5;

    ffld( thet, phi, eth, epi);
    arg=- TP* r;
    exa= cmplx( cosl( arg), sinl( arg))/ r;
    *eth= *eth* exa;
    *epi= *epi* exa;
    *erd=CPLX_00;
    return;
  } /* if( (ksymp == 1) && (cabs(ux) > .5) && (r > 1.e5) ) */

  /* computation of space and ground waves. */
  u= ux;
  u2= u* u;
  phx=- sinl( phi);
  phy= cosl( phi);
  rx= rho* phy;
  ry=- rho* phx;
  cix=CPLX_00;
  ciy=CPLX_00;
  ciz=CPLX_00;

  /* summation of field from individual segments */
  for( i = 0; i < n; i++ )
  {
    dx= cab[i];
    dy= sab[i];
    dz= salp[i];
    rix= rx- x[i];
    riy= ry- y[i];
    rhs= rix* rix+ riy* riy;
    rhp= sqrtl( rhs);

    if( rhp >= 1.0e-6)
    {
      rhx= rix/ rhp;
      rhy= riy/ rhp;
    }
    else
    {
      rhx=1.;
      rhy=0.;
    }

    calp=1.- dz* dz;
    if( calp >= 1.0e-6)
    {
      calp= sqrtl( calp);
      cbet= dx/ calp;
      sbet= dy/ calp;
      cph= rhx* cbet+ rhy* sbet;
      sph= rhy* cbet- rhx* sbet;
    }
    else
    {
      cph= rhx;
      sph= rhy;
    }

    el= PI* si[i];
    rfl=-1.;

    /* integration of (current)*(phase factor) over segment and image for */
    /* constant, sine, and cosine current distributions */
    for( k = 0; k < 2; k++ )
    {
      rfl=- rfl;
      riz= rz- z[i]* rfl;
      rxyz= sqrtl( rix* rix+ riy* riy+ riz* riz);
      rnx= rix/ rxyz;
      rny= riy/ rxyz;
      rnz= riz/ rxyz;
      omega=-( rnx* dx+ rny* dy+ rnz* dz* rfl);
      sill= omega* el;
      top= el+ sill;
      bot= el- sill;

      if( fabsl( omega) >= 1.0e-7)
	a=2.* sinl( sill)/ omega;
      else
	a=(2.- omega* omega* el* el/3.)* el;

      if( fabsl( top) >= 1.0e-7)
	too= sinl( top)/ top;
      else
	too=1.- top* top/6.;

      if( fabsl( bot) >= 1.0e-7)
	boo= sinl( bot)/ bot;
      else
	boo=1.- bot* bot/6.;

      b= el*( boo- too);
      c= el*( boo+ too);
      rr= a* air[i]+ b* bii[i]+ c* cir[i];
      ri= a* aii[i]- b* bir[i]+ c* cii[i];
      arg= TP*( x[i]* rnx+ y[i]* rny+ z[i]* rnz* rfl);
      exa= cmplx( cosl( arg), sinl( arg))* cmplx( rr, ri)/ TP;

      if( k != 1 )
      {
	xx1= exa;
	r1= rxyz;
	zmh= riz;
	continue;
      }

      xx2= exa;
      r2= rxyz;
      zph= riz;

    } /* for( k = 0; k < 2; k++ ) */

    /* call subroutine to compute the field */
    /* of segment including ground wave. */
    gwave( &erv, &ezv, &erh, &ezh, &eph);
    erh= erh* cph* calp+ erv* dz;
    eph= eph* sph* calp;
    ezh= ezh* cph* calp+ ezv* dz;
    ex= erh* rhx- eph* rhy;
    ey= erh* rhy+ eph* rhx;
    cix= cix+ ex;
    ciy= ciy+ ey;
    ciz= ciz+ ezh;

  } /* for( i = 0; i < n; i++ ) */

  arg=- TP* r;
  exa= cmplx( cosl( arg), sinl( arg));
  cix= cix* exa;
  ciy= ciy* exa;
  ciz= ciz* exa;
  rnx= rx/ r;
  rny= ry/ r;
  rnz= rz/ r;
  thx= rnz* phy;
  thy=- rnz* phx;
  thz=- rho/ r;
  *eth= cix* thx+ ciy* thy+ ciz* thz;
  *epi= cix* phx+ ciy* phy;
  *erd= cix* rnx+ ciy* rny+ ciz* rnz;

  return;
}

/*-----------------------------------------------------------------------*/

/* compute radiation pattern, gain, normalized gain */
void rdpat( void )
{
  char  *hpol[3] = { "LINEAR", "RIGHT ", "LEFT  " };
  char    hcir[] = " CIRCLE";
  char  *igtp[2] = { "----- POWER GAINS ----- ", "--- DIRECTIVE GAINS ---" };
  char  *igax[4] = { " MAJOR", " MINOR", " VERTC", " HORIZ" };
  char *igntp[5] =  { " MAJOR AXIS", "  MINOR AXIS",
    "    VERTICAL", "  HORIZONTAL", "       TOTAL " };

    char *hclif=NULL, *isens;
    int i, j, jump, itmp1, itmp2, kth, kph, itmp3, itmp4;
    long double exrm=0., exra=0., prad, gcon, gcop, gmax, pint, tmp1, tmp2;
    long double phi, pha, thet, tha, erdm=0., erda=0., ethm2, ethm, *gain = NULL;
    long double etha, ephm2, ephm, epha, tilta, emajr2, eminr2, axrat;
    long double dfaz, dfaz2, cdfaz, tstor1=0., tstor2, stilta, gnmj;
    long double gnmn, gnv, gnh, gtot, tmp3, tmp4, da, tmp5, tmp6;
    complex long double  eth, eph, erd;

    /* Allocate memory to gain buffer */
    if( inor > 0 )
      mem_alloc( (void *)&gain, nth*nph * sizeof(long double) );

    if( ifar >= 2)
    {
      fprintf( output_fp, "\n\n\n"
	  "                                 "
	  "------ FAR FIELD GROUND PARAMETERS ------\n\n" );

      jump = FALSE;
      if( ifar > 3)
      {
	fprintf( output_fp, "\n"
	    "                                        "
	    "RADIAL WIRE GROUND SCREEN\n"
	    "                                        "
	    "%5d WIRES\n"
	    "                                        "
	    "WIRE LENGTH= %8.2LF METERS\n"
	    "                                        "
	    "WIRE RADIUS= %10.3LE METERS",
	    nradl, scrwlt, scrwrt );

	if( ifar == 4)
	  jump = TRUE;

      } /* if( ifar > 3) */

      if( ! jump )
      {
	if( (ifar == 2) || (ifar == 5) )
	  hclif= hpol[0];
	if( (ifar == 3) || (ifar == 6) )
	  hclif= hcir;

	cl= clt/ wlam;
	ch= cht/ wlam;
	zrati2= csqrtl(1./ cmplx( epsr2,- sig2* wlam*59.96));

	fprintf( output_fp, "\n"
	    "                                        "
	    "%6s CLIFF\n"
	    "                                        "
	    "EDGE DISTANCE= %9.2LF METERS\n"
	    "                                        "
	    "HEIGHT= %8.2LF METERS\n"
	    "                                        "
	    "SECOND MEDIUM -\n"
	    "                                        "
	    "RELATIVE DIELECTRIC CONST.= %7.3LF\n"
	    "                                        "
	    "CONDUCTIVITY= %10.3LF MHOS",
	    hclif, clt, cht, epsr2, sig2 );

      } /* if( ! jump ) */

    } /* if( ifar >= 2) */

    if( ifar == 1)
    {
      fprintf( output_fp, "\n\n\n"
	  "                             "
	  "------- RADIATED FIELDS NEAR GROUND --------\n\n"
	  "    ------- LOCATION -------     --- E(THETA) ---    "
	  " ---- E(PHI) ----    --- E(RADIAL) ---\n"
	  "      RHO    PHI        Z           MAG    PHASE     "
	  "    MAG    PHASE        MAG     PHASE\n"
	  "    METERS DEGREES    METERS      VOLTS/M DEGREES   "
	  "   VOLTS/M DEGREES     VOLTS/M  DEGREES" );
    }
    else
    {
      itmp1=2* iax;
      itmp2= itmp1+1;

      fprintf( output_fp, "\n\n\n"
	  "                             "
	  "---------- RADIATION PATTERNS -----------\n" );

      if( rfld >= 1.0e-20)
      {
	exrm=1./ rfld;
	exra= rfld/ wlam;
	exra=-360.*( exra- floorl( exra));

	fprintf( output_fp, "\n"
	    "                             "
	    "RANGE: %13.6LE METERS\n"
	    "                             "
	    "EXP(-JKR)/R: %12.5LE AT PHASE: %7.2LF DEGREES\n",
	    rfld, exrm, exra );
      }

      fprintf( output_fp, "\n"
	  " ---- ANGLES -----     %23s      ---- POLARIZATION ----  "
	  " ---- E(THETA) ----    ----- E(PHI) ------\n"
	  "  THETA      PHI      %6s   %6s    TOTAL       AXIAL    "
	  "  TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE\n"
	  " DEGREES   DEGREES        DB       DB       DB       RATIO  "
	  " DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES",
	  igtp[ipd], igax[itmp1], igax[itmp2] );

    } /* if( ifar == 1) */

    if( (ixtyp == 0) || (ixtyp == 5) )
    {
      gcop= wlam* wlam*2.* PI/(376.73* pinr);
      prad= pinr- ploss- pnlr;
      gcon= gcop;
      if( ipd != 0)
	gcon= gcon* pinr/ prad;
    }
    else
      if( ixtyp == 4)
      {
	pinr=394.51* xpr6* xpr6* wlam* wlam;
	gcop= wlam* wlam*2.* PI/(376.73* pinr);
	prad= pinr- ploss- pnlr;
	gcon= gcop;
	if( ipd != 0)
	  gcon= gcon* pinr/ prad;
      }
      else
      {
	prad=0.;
	gcon=4.* PI/(1.+ xpr6* xpr6);
	gcop= gcon;
      }

    i=0;
    gmax=-1.e+10;
    pint=0.;
    tmp1= dph* TA;
    tmp2=.5* dth* TA;
    phi= phis- dph;

    for( kph = 1; kph <= nph; kph++ )
    {
      phi += dph;
      pha= phi* TA;
      thet= thets- dth;

      for( kth = 1; kth <= nth; kth++ )
      {
	thet += dth;
	if( (ksymp == 2) && (thet > 90.01) && (ifar != 1) )
	  continue;

	tha= thet* TA;
	if( ifar != 1)
	  ffld( tha, pha, &eth, &eph);
	else
	{
	  gfld( rfld/wlam, pha, thet/wlam,
	      &eth, &eph, &erd, zrati, ksymp);
	  erdm= cabs( erd);
	  erda= cang( erd);
	}

	ethm2= creal( eth* conjl( eth));
	ethm= sqrtl( ethm2);
	etha= cang( eth);
	ephm2= creal( eph* conjl( eph));
	ephm= sqrtl( ephm2);
	epha= cang( eph);

	/* elliptical polarization calc. */
	if( ifar != 1)
	{
	  if( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) )
	  {
	    tilta=0.;
	    emajr2=0.;
	    eminr2=0.;
	    axrat=0.;
	    isens= " ";
	  }
	  else
	  {
	    dfaz= epha- etha;
	    if( epha >= 0.)
	      dfaz2= dfaz-360.;
	    else
	      dfaz2= dfaz+360.;

	    if( fabsl(dfaz) > fabsl(dfaz2) )
	      dfaz= dfaz2;

	    cdfaz= cosl( dfaz* TA);
	    tstor1= ethm2- ephm2;
	    tstor2=2.* ephm* ethm* cdfaz;
	    tilta=.5* atan2l( tstor2, tstor1);
	    stilta= sinl( tilta);
	    tstor1= tstor1* stilta* stilta;
	    tstor2= tstor2* stilta* cosl( tilta);
	    emajr2=- tstor1+ tstor2+ ethm2;
	    eminr2= tstor1- tstor2+ ephm2;
	    if( eminr2 < 0.)
	      eminr2=0.;

	    axrat= sqrtl( eminr2/ emajr2);
	    tilta= tilta* TD;
	    if( axrat <= 1.0e-5)
	      isens= hpol[0];
	    else
	      if( dfaz <= 0.)
		isens= hpol[1];
	      else
		isens= hpol[2];

	  } /* if( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) ) */

	  gnmj= db10( gcon* emajr2);
	  gnmn= db10( gcon* eminr2);
	  gnv = db10( gcon* ethm2);
	  gnh = db10( gcon* ephm2);
	  gtot= db10( gcon*(ethm2+ ephm2) );

	  if( inor > 0)
	  {
	    i++;
	    switch( inor )
	    {
	      case 1:
		tstor1= gnmj;
		break;

	      case 2:
		tstor1= gnmn;
		break;

	      case 3:
		tstor1= gnv;
		break;

	      case 4:
		tstor1= gnh;
		break;

	      case 5:
		tstor1= gtot;
	    }

	    gain[i-1]= tstor1;
	    if( tstor1 > gmax)
	      gmax= tstor1;

	  } /* if( inor > 0) */

	  if( iavp != 0)
	  {
	    tstor1= gcop*( ethm2+ ephm2);
	    tmp3= tha- tmp2;
	    tmp4= tha+ tmp2;

	    if( kth == 1)
	      tmp3= tha;
	    else
	      if( kth == nth)
		tmp4= tha;

	    da= fabsl( tmp1*( cosl( tmp3)- cosl( tmp4)));
	    if( (kph == 1) || (kph == nph) )
	      da *=.5;
	    pint += tstor1* da;

	    if( iavp == 2)
	      continue;
	  }

	  if( iax != 1)
	  {
	    tmp5= gnmj;
	    tmp6= gnmn;
	  }
	  else
	  {
	    tmp5= gnv;
	    tmp6= gnh;
	  }

	  ethm= ethm* wlam;
	  ephm= ephm* wlam;

	  if( rfld >= 1.0e-20 )
	  {
	    ethm= ethm* exrm;
	    etha= etha+ exra;
	    ephm= ephm* exrm;
	    epha= epha+ exra;
	  }

	  fprintf( output_fp, "\n"
	      " %7.2LF %9.2LF  %8.2LF %8.2LF %8.2LF %11.4LF"
	      " %9.2LF %6s %11.4LE %9.2LF %11.4LE %9.2LF",
	      thet, phi, tmp5, tmp6, gtot, axrat,
	      tilta, isens, ethm, etha, ephm, epha );

	  if( iplp1 != 3)
	    continue;

	  if( iplp3 != 0)
	  {
	    if( iplp2 == 1 )
	    {
	      if( iplp3 == 1 )
		fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", thet, ethm, etha );
	      else
		if( iplp3 == 2 )
		  fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", thet, ephm, epha );
	    }

	    if( iplp2 == 2 )
	    {
	      if( iplp3 == 1 )
		fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", phi, ethm, etha );
	      else
		if( iplp3 == 2 )
		  fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE\n", phi, ephm, epha );
	    }
	  }

	  if( iplp4 == 0 )
	    continue;

	  if( iplp2 == 1 )
	  {
	    switch( iplp4 )
	    {
	      case 1:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", thet, tmp5 );
		break;
	      case 2:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", thet, tmp6 );
		break;
	      case 3:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", thet, gtot );
	    }
	  }

	  if( iplp2 == 2 )
	  {
	    switch( iplp4 )
	    {
	      case 1:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", phi, tmp5 );
		break;
	      case 2:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", phi, tmp6 );
		break;
	      case 3:
		fprintf( plot_fp, "%12.4LE %12.4LE\n", phi, gtot );
	    }
	  }

	  continue;
	} /* if( ifar != 1) */

	fprintf( output_fp, "\n"
	    " %9.2LF %7.2LF %9.2LF  %11.4LE %7.2LF  %11.4LE %7.2LF  %11.4LE %7.2LF",
	    rfld, phi, thet, ethm, etha, ephm, epha, erdm, erda );

      } /* for( kth = 1; kth <= nth; kth++ ) */

    } /* for( kph = 1; kph <= nph; kph++ ) */

    if( iavp != 0)
    {
      tmp3= thets* TA;
      tmp4= tmp3+ dth* TA* (long double)( nth-1);
      tmp3= fabsl( dph* TA* (long double)( nph-1)*( cosl( tmp3)- cosl( tmp4)));
      pint /= tmp3;
      tmp3 /= PI;

      fprintf( output_fp, "\n\n\n"
	  "  AVERAGE POWER GAIN: %11.4LE - SOLID ANGLE"
	  " USED IN AVERAGING: (%+7.4LF)*PI STERADIANS",
	  pint, tmp3 );
    }

    if( inor == 0)
      return;

    if( fabsl( gnor) > 1.0e-20)
      gmax= gnor;
    itmp1=( inor-1);

    fprintf( output_fp,	"\n\n\n"
	"                             "
	" ---------- NORMALIZED GAIN ----------\n"
	"                                      %6s GAIN\n"
	"                                  "
	" NORMALIZATION FACTOR: %.2LF db\n\n"
	"    ---- ANGLES ----                ---- ANGLES ----"
	"                ---- ANGLES ----\n"
	"    THETA      PHI        GAIN      THETA      PHI  "
	"      GAIN      THETA      PHI       GAIN\n"
	"   DEGREES   DEGREES        DB     DEGREES   DEGREES "
	"       DB     DEGREES   DEGREES       DB",
	igntp[itmp1], gmax );

    itmp2= nph* nth;
    itmp1=( itmp2+2)/3;
    itmp2= itmp1*3- itmp2;
    itmp3= itmp1;
    itmp4=2* itmp1;

    if( itmp2 == 2)
      itmp4--;

    for( i = 0; i < itmp1; i++ )
    {
      itmp3++;
      itmp4++;
      j= i/ nth;
      tmp1= thets+ (long double)( i - j*nth )* dth;
      tmp2= phis+ (long double)(j)* dph;
      j=( itmp3-1)/ nth;
      tmp3= thets+ (long double)( itmp3- j* nth-1)* dth;
      tmp4= phis+ (long double)(j)* dph;
      j=( itmp4-1)/ nth;
      tmp5= thets+ (long double)( itmp4- j* nth-1)* dth;
      tmp6= phis+ (long double)(j)* dph;
      tstor1= gain[i]- gmax;

      if( ((i+1) == itmp1) && (itmp2 != 0) )
      {
	if( itmp2 != 2)
	{
	  tstor2= gain[itmp3-1]- gmax;
	  fprintf( output_fp, "\n"
	      " %9.2LF %9.2LF %9.2LF   %9.2LF %9.2LF %9.2LF   ",
	      tmp1, tmp2, tstor1, tmp3, tmp4, tstor2 );
	  return;
	}

	fprintf( output_fp, "\n"
	    " %9.2LF %9.2LF %9.2LF   ",
	    tmp1, tmp2, tstor1 );
	return;

      } /* if( ((i+1) == itmp1) && (itmp2 != 0) ) */

      tstor2= gain[itmp3-1]- gmax;
      pint= gain[itmp4-1]- gmax;

      fprintf( output_fp, "\n"
	  " %9.2LF %9.2LF %9.2LF   %9.2LF %9.2LF %9.2LF   %9.2LF %9.2LF %9.2LF",
	  tmp1, tmp2, tstor1, tmp3, tmp4, tstor2, tmp5, tmp6, pint );

    } /* for( i = 0; i < itmp1; i++ ) */

    free_ptr( (void *)&gain );

    return;
}

/*-----------------------------------------------------------------------*/

