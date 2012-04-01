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

/* common  /segj/ */
extern int *jco, jsno, nscon, maxcon; /* Max. no. connections */
extern long double *ax, *bx, *cx;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/* common  /plot/ */
extern int iplp1, iplp2, iplp3, iplp4;

/*-------------------------------------------------------------------*/

/* arc generates segment geometry data for an arc of ns segments */
void arc( int itg, int ns, long double rada,
    long double ang1, long double ang2, long double rad )
{
  int ist, i, mreq;
  long double ang, dang, xs1, xs2, zs1, zs2;

  ist= n;
  n += ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  if( fabsl( ang2- ang1) < 360.00001)
  {
    /* Reallocate tags buffer */
    mem_realloc( (void *)&itag, (n+m) * sizeof(int) );

    /* Reallocate wire buffers */
    mreq = n * sizeof(long double);
    mem_realloc( (void *)&x, mreq );
    mem_realloc( (void *)&y, mreq );
    mem_realloc( (void *)&z, mreq );
    mem_realloc( (void *)&x2, mreq );
    mem_realloc( (void *)&y2, mreq );
    mem_realloc( (void *)&z2, mreq );
    mem_realloc( (void *)&bi, mreq );

    ang= ang1* TA;
    dang=( ang2- ang1)* TA/ ns;
    xs1= rada* cosl( ang);
    zs1= rada* sinl( ang);

    for( i = ist; i < n; i++ )
    {
      ang += dang;
      xs2= rada* cosl( ang);
      zs2= rada* sinl( ang);
      x[i]= xs1;
      y[i]=0.;
      z[i]= zs1;
      x2[i]= xs2;
      y2[i]=0.;
      z2[i]= zs2;
      xs1= xs2;
      zs1= zs2;
      bi[i]= rad;
      itag[i]= itg;

    } /* for( i = ist; i < n; i++ ) */

  } /* if( fabsl( ang2- ang1) < 360.00001) */
  else
  {
    fprintf( output_fp, "\n  ERROR -- ARC ANGLE EXCEEDS 360 DEGREES");
    stop(-1);
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* connect sets up segment connection data in arrays icon1 and */
/* icon2 by searching for segment ends that are in contact. */
void conect( int ignd )
{
  int i, iz, ic, j, jx, ix, ixx, iseg, iend, jend, nsflg, jump, ipf;
  long double sep=0., xi1, yi1, zi1, xi2, yi2, zi2;
  long double slen, xa, ya, za, xs, ys, zs;

  nscon= -1;
  maxcon = 1;

  if( ignd != 0)
  {
    fprintf( output_fp, "\n\n     GROUND PLANE SPECIFIED." );

    if( ignd > 0)
      fprintf( output_fp,
	  "\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL"
	  " BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n" );

    if( ipsym == 2)
    {
      np=2* np;
      mp=2* mp;
    }

    if( abs( ipsym) > 2 )
    {
      np= n;
      mp= m;
    }

    /*** possibly should be error condition?? **/
    if( np > n)
    {
      fprintf( output_fp,
	  "\n ERROR: NP > N IN CONECT()" );
      stop(-1);
    }

    if( (np == n) && (mp == m) )
      ipsym=0;

  } /* if( ignd != 0) */

  if( n != 0)
  {
    /* Allocate memory to connections */
    mem_alloc( (void *)&icon1, (n+m) * sizeof(int) );
    mem_alloc( (void *)&icon2, (n+m) * sizeof(int) );

    for( i = 0; i < n; i++ )
    {
      iz = i+1;
      xi1= x[i];
      yi1= y[i];
      zi1= z[i];
      xi2= x2[i];
      yi2= y2[i];
      zi2= z2[i];
      slen= sqrtl( (xi2- xi1)*(xi2- xi1) + (yi2- yi1) *
	  (yi2- yi1) + (zi2- zi1)*(zi2- zi1) ) * SMIN;

      /* determine connection data for end 1 of segment. */
      jump = FALSE;
      if( ignd > 0)
      {
	if( zi1 <= -slen)
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR -- SEGMENT"
	      " %d EXTENDS BELOW GROUND", iz );
	  stop(-1);
	}

	if( zi1 <= slen)
	{
	  icon1[i]= iz;
	  z[i]=0.;
	  jump = TRUE;

	} /* if( zi1 <= slen) */

      } /* if( ignd > 0) */

      if( ! jump )
      {
	ic= i;
	for( j = 1; j < n; j++)
	{
	  ic++;
	  if( ic >= n)
	    ic=0;

	  sep= fabsl( xi1- x[ic])+ fabsl(yi1- y[ic])+ fabsl(zi1- z[ic]);
	  if( sep <= slen)
	  {
	    icon1[i]= -(ic+1);
	    break;
	  }

	  sep= fabsl( xi1- x2[ic])+ fabsl(yi1- y2[ic])+ fabsl(zi1- z2[ic]);
	  if( sep <= slen)
	  {
	    icon1[i]= (ic+1);
	    break;
	  }

	} /* for( j = 1; j < n; j++) */

	if( ((iz > 0) || (icon1[i] <= PCHCON)) && (sep > slen) )
	  icon1[i]=0;

      } /* if( ! jump ) */

      /* determine connection data for end 2 of segment. */
      if( (ignd > 0) || jump )
      {
	if( zi2 <= -slen)
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR -- SEGMENT"
	      " %d EXTENDS BELOW GROUND", iz );
	  stop(-1);
	}

	if( zi2 <= slen)
	{
	  if( icon1[i] == iz )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR -- SEGMENT"
		" %d LIES IN GROUND PLANE", iz );
	    stop(-1);
	  }

	  icon2[i]= iz;
	  z2[i]=0.;
	  continue;

	} /* if( zi2 <= slen) */

      } /* if( ignd > 0) */

      ic= i;
      for( j = 1; j < n; j++ )
      {
	ic++;
	if( ic >= n)
	  ic=0;

	sep= fabsl(xi2- x[ic])+ fabsl(yi2- y[ic])+ fabsl(zi2- z[ic]);
	if( sep <= slen)
	{
	  icon2[i]= (ic+1);
	  break;
	}

	sep= fabsl(xi2- x2[ic])+ fabsl(yi2- y2[ic])+ fabsl(zi2- z2[ic]);
	if( sep <= slen)
	{
	  icon2[i]= -(ic+1);
	  break;
	}

      } /* for( j = 1; j < n; j++ ) */

      if( ((iz > 0) || (icon2[i] <= PCHCON)) && (sep > slen) )
	icon2[i]=0;

    } /* for( i = 0; i < n; i++ ) */

    /* find wire-surface connections for new patches */
    if( m != 0)
    {
      ix = -1;
      i = 0;
      while( ++i <= m )
      {
	ix++;
	xs= px[ix];
	ys= py[ix];
	zs= pz[ix];

	for( iseg = 0; iseg < n; iseg++ )
	{
	  xi1= x[iseg];
	  yi1= y[iseg];
	  zi1= z[iseg];
	  xi2= x2[iseg];
	  yi2= y2[iseg];
	  zi2= z2[iseg];

	  /* for first end of segment */
	  slen=( fabsl(xi2- xi1)+ fabsl(yi2- yi1)+ fabsl(zi2- zi1))* SMIN;
	  sep= fabsl(xi1- xs)+ fabsl(yi1- ys)+ fabsl(zi1- zs);

	  /* connection - divide patch into 4 patches at present array loc. */
	  if( sep <= slen)
	  {
	    icon1[iseg]=PCHCON+ i;
	    ic=0;
	    subph( i, ic );
	    break;
	  }

	  sep= fabsl(xi2- xs)+ fabsl(yi2- ys)+ fabsl(zi2- zs);
	  if( sep <= slen)
	  {
	    icon2[iseg]=PCHCON+ i;
	    ic=0;
	    subph( i, ic );
	    break;
	  }

	} /* for( iseg = 0; iseg < n; iseg++ ) */

      } /* while( ++i <= m ) */

    } /* if( m != 0) */

  } /* if( n != 0) */

  fprintf( output_fp, "\n\n"
      "     TOTAL SEGMENTS USED: %d   SEGMENTS IN A"
      " SYMMETRIC CELL: %d   SYMMETRY FLAG: %d",
      n, np, ipsym );

  if( m > 0)
    fprintf( output_fp,	"\n"
	"       TOTAL PATCHES USED: %d   PATCHES"
	" IN A SYMMETRIC CELL: %d",  m, mp );

  iseg=( n+ m)/( np+ mp);
  if( iseg != 1)
  {
    /*** may be error condition?? ***/
    if( ipsym == 0 )
    {
      fprintf( output_fp,
	  "\n  ERROR: IPSYM=0 IN CONECT()" );
      stop(-1);
    }

    if( ipsym < 0 )
      fprintf( output_fp,
	  "\n  STRUCTURE HAS %d FOLD ROTATIONAL SYMMETRY\n", iseg );
    else
    {
      ic= iseg/2;
      if( iseg == 8)
	ic=3;
      fprintf( output_fp,
	  "\n  STRUCTURE HAS %d PLANES OF SYMMETRY\n", ic );
    } /* if( ipsym < 0 ) */

  } /* if( iseg == 1) */

  if( n == 0)
    return;

  /* Allocate to connection buffers */
  mem_alloc( (void *)&jco, maxcon * sizeof(int) );

  /* adjust connected seg. ends to exactly coincide.  print junctions */
  /* of 3 or more seg.  also find old seg. connecting to new seg. */
  iseg = 0;
  ipf = FALSE;
  for( j = 0; j < n; j++ )
  {
    jx = j+1;
    iend=-1;
    jend=-1;
    ix= icon1[j];
    ic=1;
    jco[0]= -jx;
    xa= x[j];
    ya= y[j];
    za= z[j];

    while( TRUE )
    {
      if( (ix != 0) && (ix != (j+1)) && (ix <= PCHCON) )
      {
	nsflg=0;

	do
	{
	  if( ix == 0 )
	  {
	    fprintf( output_fp,
		"\n  CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT: %d", ix );
	    stop(-1);
	  }

	  if( ix < 0 )
	    ix= -ix;
	  else
	    jend= -jend;

	  jump = FALSE;

	  if( ix == jx )
	    break;

	  if( ix < jx )
	  {
	    jump = TRUE;
	    break;
	  }

	  /* Record max. no. of connections */
	  ic++;
	  if( ic >= maxcon )
	  {
	    maxcon = ic+1;
	    mem_realloc( (void *)&jco, maxcon * sizeof(int) );
	  }
	  jco[ic-1]= ix* jend;

	  if( ix > 0)
	    nsflg=1;

	  ixx = ix-1;
	  if( jend != 1)
	  {
	    xa= xa+ x[ixx];
	    ya= ya+ y[ixx];
	    za= za+ z[ixx];
	    ix= icon1[ixx];
	    continue;
	  }

	  xa= xa+ x2[ixx];
	  ya= ya+ y2[ixx];
	  za= za+ z2[ixx];
	  ix= icon2[ixx];

	} /* do */
	while( ix != 0 );

	if( jump && (iend == 1) )
	  break;
	else
	  if( jump )
	  {
	    iend=1;
	    jend=1;
	    ix= icon2[j];
	    ic=1;
	    jco[0]= jx;
	    xa= x2[j];
	    ya= y2[j];
	    za= z2[j];
	    continue;
	  }

	sep= (long double)ic;
	xa= xa/ sep;
	ya= ya/ sep;
	za= za/ sep;

	for( i = 0; i < ic; i++ )
	{
	  ix= jco[i];
	  if( ix <= 0)
	  {
	    ix=- ix;
	    ixx = ix-1;
	    x[ixx]= xa;
	    y[ixx]= ya;
	    z[ixx]= za;
	    continue;
	  }

	  ixx = ix-1;
	  x2[ixx]= xa;
	  y2[ixx]= ya;
	  z2[ixx]= za;

	} /* for( i = 0; i < ic; i++ ) */

	if( ic >= 3)
	{
	  if( ! ipf )
	  {
	    fprintf( output_fp, "\n\n"
		"    ---------- MULTIPLE WIRE JUNCTIONS ----------\n"
		"    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)" );
	    ipf = TRUE;
	  }

	  iseg++;
	  fprintf( output_fp, "\n   %5d      ", iseg );

	  for( i = 1; i <= ic; i++ )
	  {
	    fprintf( output_fp, "%5d", jco[i-1] );
	    if( !(i % 20) )
	      fprintf( output_fp, "\n              " );
	  }

	} /* if( ic >= 3) */

      } /*if( (ix != 0) && (ix != j) && (ix <= PCHCON) ) */

      if( iend == 1)
	break;

      iend=1;
      jend=1;
      ix= icon2[j];
      ic=1;
      jco[0]= jx;
      xa= x2[j];
      ya= y2[j];
      za= z2[j];

    } /* while( TRUE ) */

  } /* for( j = 0; j < n; j++ ) */

  mem_alloc( (void *)&ax, maxcon * sizeof(long double) );
  mem_alloc( (void *)&bx, maxcon * sizeof(long double) );
  mem_alloc( (void *)&cx, maxcon * sizeof(long double) );

  return;
}

/*-----------------------------------------------------------------------*/

/* datagn is the main routine for input of geometry data. */
void datagn( void )
{
  char gm[3];
  char ifx[2] = {'*', 'X'}, ify[2]={'*','Y'}, ifz[2]={'*','Z'};
  char ipt[4] = { 'P', 'R', 'T', 'Q' };

  /* input card mnemonic list */
  /* "XT" stands for "exit", added for testing */
#define GM_NUM  12
  char *atst[GM_NUM] =
  {
    "GW", "GX", "GR", "GS", "GE", "GM", \
    "SP", "SM", "GA", "SC", "GH", "GF"
  };

  int nwire, isct, iphd, i1, i2, itg, iy, iz, mreq;
  int ix, i, ns, gm_num; /* geometry card id as a number */
  long double rad, xs1, xs2, ys1, ys2, zs1, zs2, x4=0, y4=0, z4=0;
  long double x3=0, y3=0, z3=0, xw1, xw2, yw1, yw2, zw1, zw2;
  long double dummy;

  ipsym=0;
  nwire=0;
  n=0;
  np=0;
  m=0;
  mp=0;
  isct=0;
  iphd = FALSE;

  /* read geometry data card and branch to */
  /* section for operation requested */
  do
  {
    readgm( gm, &itg, &ns, &xw1, &yw1, &zw1, &xw2, &yw2, &zw2, &rad);

    /* identify card id mnemonic */
    for( gm_num = 0; gm_num < GM_NUM; gm_num++ )
      if( strncmp( gm, atst[gm_num], 2) == 0 )
	break;

    if( iphd == FALSE )
    {
      fprintf( output_fp, "\n\n\n"
	  "                               "
	  "-------- STRUCTURE SPECIFICATION --------\n"
	  "                                     "
	  "COORDINATES MUST BE INPUT IN\n"
	  "                                     "
	  "METERS OR BE SCALED TO METERS\n"
	  "                                     "
	  "BEFORE STRUCTURE INPUT IS ENDED\n" );

      fprintf( output_fp, "\n"
	  "  WIRE                                           "
	  "                                      SEG FIRST  LAST  TAG\n"
	  "   No:        X1         Y1         Z1         X2      "
	  "   Y2         Z2       RADIUS   No:   SEG   SEG  No:" );

      iphd=1;
    }

    if( gm_num != 10 )
      isct=0;

    switch( gm_num )
    {
      case 0: /* "gw" card, generate segment data for straight wire. */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d  %10.4LF %10.4LF %10.4LF %10.4LF"
	    " %10.4LF %10.4LF %10.4LF %5d %5d %5d %4d",
	    nwire, xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, i1, i2, itg );

	if( rad != 0)
	{
	  xs1=1.;
	  ys1=1.;
	}
	else
	{
	  readgm( gm, &ix, &iy, &xs1, &ys1, &zs1,
	      &dummy, &dummy, &dummy, &dummy);

	  if( strcmp(gm, "GC" ) != 0 )
	  {
	    fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	    stop(-1);
	  }

	  fprintf( output_fp,
	      "\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: %9.5LF\n"
	      "                                 "
	      "RADIUS FROM: %9.5LF TO: %9.5LF", xs1, ys1, zs1 );

	  if( (ys1 == 0) || (zs1 == 0) )
	  {
	    fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	    stop(-1);
	  }

	  rad= ys1;
	  ys1= powl( (zs1/ys1), (1./(ns-1.)) );
	}

	wire( xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1, ns, itg);

	continue;

	/* reflect structure along x,y, or z */
	/* axes or rotate to form cylinder.  */
      case 1: /* "gx" card */

	iy= ns/10;
	iz= ns- iy*10;
	ix= iy/10;
	iy= iy- ix*10;

	if( ix != 0)
	  ix=1;
	if( iy != 0)
	  iy=1;
	if( iz != 0)
	  iz=1;

	fprintf( output_fp,
	    "\n  STRUCTURE REFLECTED ALONG THE AXES %c %c %c"
	    " - TAGS INCREMENTED BY %d\n",
	    ifx[ix], ify[iy], ifz[iz], itg );

	reflc( ix, iy, iz, itg, ns);

	continue;

      case 2: /* "gr" card */

	fprintf( output_fp,
	    "\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
	    " - LABELS INCREMENTED BY %d\n", ns, itg );

	ix=-1;
	iz = 0;
	reflc( ix, iy, iz, itg, ns);

	continue;

      case 3: /* "gs" card, scale structure dimensions by factor xw1. */

	if( n > 0)
	{
	  for( i = 0; i < n; i++ )
	  {
	    x[i]= x[i]* xw1;
	    y[i]= y[i]* xw1;
	    z[i]= z[i]* xw1;
	    x2[i]= x2[i]* xw1;
	    y2[i]= y2[i]* xw1;
	    z2[i]= z2[i]* xw1;
	    bi[i]= bi[i]* xw1;
	  }
	} /* if( n >= n2) */

	if( m > 0)
	{
	  yw1= xw1* xw1;
	  for( i = 0; i < m; i++ )
	  {
	    px[i]= px[i]* xw1;
	    py[i]= py[i]* xw1;
	    pz[i]= pz[i]* xw1;
	    pbi[i]= pbi[i]* yw1;
	  }
	} /* if( m >= m2) */

	fprintf( output_fp,
	    "\n     STRUCTURE SCALED BY FACTOR: %10.5LF", xw1 );

	continue;

      case 4: /* "ge" card, terminate structure geometry input. */

	if( ns != 0)
	{
	  iplp1=1;
	  iplp2=1;
	}

	conect( itg);

	if( n != 0)
	{
	  /* Allocate wire buffers */
	  mreq = n * sizeof(long double);
	  mem_alloc( (void *)&si, mreq );
	  mem_alloc( (void *)&sab, mreq );
	  mem_alloc( (void *)&cab, mreq );
	  mem_alloc( (void *)&salp, mreq );

	  fprintf( output_fp, "\n\n\n"
	      "                              "
	      " ---------- SEGMENTATION DATA ----------\n"
	      "                                       "
	      " COORDINATES IN METERS\n"
	      "                           "
	      " I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n" );

	  fprintf( output_fp, "\n"
	      "   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION"
	      " ANGLES    WIRE    CONNECTION DATA   TAG\n"
	      "   No:       X         Y         Z      LENGTH     ALPHA     "
	      " BETA    RADIUS    I-     I    I+   No:" );

	  for( i = 0; i < n; i++ )
	  {
	    xw1= x2[i]- x[i];
	    yw1= y2[i]- y[i];
	    zw1= z2[i]- z[i];
	    x[i]=( x[i]+ x2[i])*.5;
	    y[i]=( y[i]+ y2[i])*.5;
	    z[i]=( z[i]+ z2[i])*.5;
	    xw2= xw1* xw1+ yw1* yw1+ zw1* zw1;
	    yw2= sqrtl( xw2);
	    yw2=( xw2/ yw2+ yw2)*.5;
	    si[i]= yw2;
	    cab[i]= xw1/ yw2;
	    sab[i]= yw1/ yw2;
	    xw2= zw1/ yw2;

	    if( xw2 > 1.)
	      xw2=1.;
	    if( xw2 < -1.)
	      xw2=-1.;

	    salp[i]= xw2;
	    xw2= asinl( xw2)* TD;
	    yw2= atan2l( yw1, xw1)* TD;

	    fprintf( output_fp, "\n"
		" %5d %9.4LF %9.4LF %9.4LF %9.4LF"
		" %9.4LF %9.4LF %9.4LF %5d %5d %5d %5d",
		i+1, x[i], y[i], z[i], si[i], xw2, yw2,
		bi[i], icon1[i], i+1, icon2[i], itag[i] );

	    if( iplp1 == 1)
	      fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE "
		  "%12.4LE %12.4LE %12.4LE %12.4LE %5d %5d %5d\n",
		  x[i],y[i],z[i],si[i],xw2,yw2,bi[i],icon1[i],i+1,icon2[i] );

	    if( (si[i] <= 1.e-20) || (bi[i] <= 0.) )
	    {
	      fprintf( output_fp, "\n SEGMENT DATA ERROR" );
	      stop(-1);
	    }

	  } /* for( i = 0; i < n; i++ ) */

	} /* if( n != 0) */

	if( m != 0)
	{
	  fprintf( output_fp, "\n\n\n"
	      "                                   "
	      " --------- SURFACE PATCH DATA ---------\n"
	      "                                            "
	      " COORDINATES IN METERS\n\n"
	      " PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      "
	      " PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n"
	      "  No:       X          Y          Z          X        Y        Z      "
	      " AREA         X1       Y1       Z1        X2       Y2      Z2" );

	  for( i = 0; i < m; i++ )
	  {
	    xw1=( t1y[i]* t2z[i]- t1z[i]* t2y[i])* psalp[i];
	    yw1=( t1z[i]* t2x[i]- t1x[i]* t2z[i])* psalp[i];
	    zw1=( t1x[i]* t2y[i]- t1y[i]* t2x[i])* psalp[i];

	    fprintf( output_fp, "\n"
		" %4d %10.5LF %10.5LF %10.5LF  %8.4LF %8.4LF %8.4LF"
		" %10.5LF  %8.4LF %8.4LF %8.4LF  %8.4LF %8.4LF %8.4LF",
		i+1, px[i], py[i], pz[i], xw1, yw1, zw1, pbi[i],
		t1x[i], t1y[i], t1z[i], t2x[i], t2y[i], t2z[i] );

	  } /* for( i = 0; i < m; i++ ) */

	} /* if( m == 0) */

	npm  = n+m;
	np2m = n+2*m;
	np3m = n+3*m;

	return;

	/* "gm" card, move structure or reproduce */
	/* original structure in new positions.   */
      case 5:

	fprintf( output_fp,
	    "\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n"
	    "   %3d %5d %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
	    itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

	xw1= xw1* TA;
	yw1= yw1* TA;
	zw1= zw1* TA;

	move( xw1, yw1, zw1, xw2, yw2, zw2, (int)( rad+.5), ns, itg);
	continue;

      case 6: /* "sp" card, generate single new patch */

	i1= m+1;
	ns++;

	if( itg != 0)
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	fprintf( output_fp, "\n"
	    " %5d%c %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
	    i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

	if( (ns == 2) || (ns == 4) )
	  isct=1;

	if( ns > 1)
	{
	  readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

	  if( (ns == 2) || (itg > 0) )
	  {
	    x4= xw1+ x3- xw2;
	    y4= yw1+ y3- yw2;
	    z4= zw1+ z3- zw2;
	  }

	  fprintf( output_fp, "\n"
	      "      %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
	      x3, y3, z3, x4, y4, z4 );

	  if( strcmp(gm, "SC") != 0 )
	  {
	    fprintf( output_fp, "\n  PATCH DATA ERROR" );
	    stop(-1);
	  }

	} /* if( ns > 1) */
	else
	{
	  xw2= xw2* TA;
	  yw2= yw2* TA;
	}

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 7: /* "sm" card, generate multiple-patch surface */

	i1= m+1;
	fprintf( output_fp, "\n"
	    " %5d%c %10.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF"
	    "     SURFACE - %d BY %d PATCHES",
	    i1, ipt[1], xw1, yw1, zw1, xw2, yw2, zw2, itg, ns );

	if( (itg < 1) || (ns < 1) )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

	if( (ns == 2) || (itg > 0) )
	{
	  x4= xw1+ x3- xw2;
	  y4= yw1+ y3- yw2;
	  z4= zw1+ z3- zw2;
	}

	fprintf( output_fp, "\n"
	    "      %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
	    x3, y3, z3, x4, y4, z4 );

	if( strcmp(gm, "SC" ) != 0 )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 8: /* "ga" card, generate segment data for wire arc */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d  ARC RADIUS: %9.5LF  FROM: %8.3LF TO: %8.3LF DEGREES"
	    "       %11.5LF %5d %5d %5d %4d",
	    nwire, xw1, yw1, zw1, xw2, ns, i1, i2, itg );

	arc( itg, ns, xw1, yw1, zw1, xw2);

	continue;

      case 9: /* "sc" card */

	if( isct == 0)
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	i1= m+1;
	ns++;

	if( (itg != 0) || ((ns != 2) && (ns != 4)) )
	{
	  fprintf( output_fp, "\n  PATCH DATA ERROR" );
	  stop(-1);
	}

	xs1= x4;
	ys1= y4;
	zs1= z4;
	xs2= x3;
	ys2= y3;
	zs2= z3;
	x3= xw1;
	y3= yw1;
	z3= zw1;

	if( ns == 4)
	{
	  x4= xw2;
	  y4= yw2;
	  z4= zw2;
	}

	xw1= xs1;
	yw1= ys1;
	zw1= zs1;
	xw2= xs2;
	yw2= ys2;
	zw2= zs2;

	if( ns != 4)
	{
	  x4= xw1+ x3- xw2;
	  y4= yw1+ y3- yw2;
	  z4= zw1+ z3- zw2;
	}

	fprintf( output_fp, "\n"
	    " %5d%c %10.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
	    i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

	fprintf( output_fp, "\n"
	    "      %11.5LF %11.5LF %11.5LF  %11.5LF %11.5LF %11.5LF",
	    x3, y3, z3, x4, y4, z4 );

	patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

	continue;

      case 10: /* "gh" card, generate helix */

	nwire++;
	i1= n+1;
	i2= n+ ns;

	fprintf( output_fp, "\n"
	    " %5d HELIX STRUCTURE - SPACING OF TURNS: %8.3LF AXIAL"
	    " LENGTH: %8.3LF  %8.3LF %5d %5d %5d %4d\n      "
	    " RADIUS X1:%8.3LF Y1:%8.3LF X2:%8.3LF Y2:%8.3LF ",
	    nwire, xw1, yw1, rad, ns, i1, i2, itg, zw1, xw2, yw2, zw2 );

	helix( xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, itg);

	continue;

      case 11: /* "gf" card, not supported */
	abort_on_error(-5);

      default: /* error message */

	fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
	fprintf( output_fp, "\n"
	    " %2s %3d %5d %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
	    gm, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

	stop(-1);

    } /* switch( gm_num ) */

  } /* do */
  while( TRUE );

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine helix generates segment geometry */
/* data for a helix of ns segments */
void helix( long double s, long double hl, long double a1, long double b1,
    long double a2, long double b2, long double rad, int ns, int itg )
{
  int ist, i, mreq;
  long double turns, zinc, copy, sangle, hdia, turn, pitch, hmaj, hmin;

  ist= n;
  n += ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  turns= fabsl( hl/ s);
  zinc= fabsl( hl/ ns);

  /* Reallocate tags buffer */
  mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = n * sizeof(long double);
  mem_realloc( (void *)&x, mreq );
  mem_realloc( (void *)&y, mreq );
  mem_realloc( (void *)&z, mreq );
  mem_realloc( (void *)&x2, mreq );
  mem_realloc( (void *)&y2, mreq );
  mem_realloc( (void *)&z2, mreq );
  mem_realloc( (void *)&bi, mreq );

  z[ist]=0.;
  for( i = ist; i < n; i++ )
  {
    bi[i]= rad;
    itag[i]= itg;

    if( i != ist )
      z[i]= z[i-1]+ zinc;

    z2[i]= z[i]+ zinc;

    if( a2 == a1)
    {
      if( b1 == 0.)
	b1= a1;

      x[i]= a1* cosl(2.* PI* z[i]/ s);
      y[i]= b1* sinl(2.* PI* z[i]/ s);
      x2[i]= a1* cosl(2.* PI* z2[i]/ s);
      y2[i]= b1* sinl(2.* PI* z2[i]/ s);
    }
    else
    {
      if( b2 == 0.)
	b2= a2;

      x[i]=( a1+( a2- a1)* z[i]/ fabsl( hl))* cosl(2.* PI* z[i]/ s);
      y[i]=( b1+( b2- b1)* z[i]/ fabsl( hl))* sinl(2.* PI* z[i]/ s);
      x2[i]=( a1+( a2- a1)* z2[i]/ fabsl( hl))* cosl(2.* PI* z2[i]/ s);
      y2[i]=( b1+( b2- b1)* z2[i]/ fabsl( hl))* sinl(2.* PI* z2[i]/ s);

    } /* if( a2 == a1) */

    if( hl > 0.)
      continue;

    copy= x[i];
    x[i]= y[i];
    y[i]= copy;
    copy= x2[i];
    x2[i]= y2[i];
    y2[i]= copy;

  } /* for( i = ist; i < n; i++ ) */

  if( a2 != a1)
  {
    sangle= atanl( a2/( fabsl( hl)+( fabsl( hl)* a1)/( a2- a1)));
    fprintf( output_fp,
	"\n       THE CONE ANGLE OF THE SPIRAL IS %10.4LF", sangle );
    return;
  }

  if( a1 == b1)
  {
    hdia=2.* a1;
    turn= hdia* PI;
    pitch= atanl( s/( PI* hdia));
    turn= turn/ cosl( pitch);
    pitch=180.* pitch/ PI;
  }
  else
  {
    if( a1 >= b1)
    {
      hmaj=2.* a1;
      hmin=2.* b1;
    }
    else
    {
      hmaj=2.* b1;
      hmin=2.* a1;
    }

    hdia= sqrtl(( hmaj*hmaj+ hmin*hmin)/2* hmaj);
    turn=2.* PI* hdia;
    pitch=(180./ PI)* atanl( s/( PI* hdia));

  } /* if( a1 == b1) */

  fprintf( output_fp, "\n"
      "       THE PITCH ANGLE IS: %.4LF    THE LENGTH OF WIRE/TURN IS: %.4LF",
      pitch, turn );

  return;
}

/*-----------------------------------------------------------------------*/

/* isegno returns the segment number of the mth segment having the */
/* tag number itagi.  if itagi=0 segment number m is returned. */
int isegno( int itagi, int mx)
{
  int icnt, i, iseg;

  if( mx <= 0)
  {
    fprintf( output_fp,
	"\n  CHECK DATA, PARAMETER SPECIFYING SEGMENT"
	" POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO" );
    stop(-1);
  }

  icnt=0;
  if( itagi == 0)
  {
    iseg = mx;
    return( iseg );
  }

  if( n > 0)
  {
    for( i = 0; i < n; i++ )
    {
      if( itag[i] != itagi )
	continue;

      icnt++;
      if( icnt == mx)
      {
	iseg= i+1;
	return( iseg );
      }

    } /* for( i = 0; i < n; i++ ) */

  } /* if( n > 0) */

  fprintf( output_fp, "\n\n"
      "  NO SEGMENT HAS AN ITAG OF %d",  itagi );
  stop(-1);

  return(0);
}

/*-----------------------------------------------------------------------*/

/* subroutine move moves the structure with respect to its */
/* coordinate system or reproduces structure in new positions. */
/* structure is rotated about x,y,z axes by rox,roy,roz */
/* respectively, then shifted by xs,ys,zs */
void move( long double rox, long double roy, long double roz, long double xs,
    long double ys, long double zs, int its, int nrpt, int itgi )
{
  int nrp, ix, i1, k, ir, i, ii, mreq;
  long double sps, cps, sth, cth, sph, cph, xx, xy;
  long double xz, yx, yy, yz, zx, zy, zz, xi, yi, zi;

  if( fabsl( rox)+ fabsl( roy) > 1.0e-10)
    ipsym= ipsym*3;

  sps= sinl( rox);
  cps= cosl( rox);
  sth= sinl( roy);
  cth= cosl( roy);
  sph= sinl( roz);
  cph= cosl( roz);
  xx= cph* cth;
  xy= cph* sth* sps- sph* cps;
  xz= cph* sth* cps+ sph* sps;
  yx= sph* cth;
  yy= sph* sth* sps+ cph* cps;
  yz= sph* sth* cps- cph* sps;
  zx=- sth;
  zy= cth* sps;
  zz= cth* cps;

  if( nrpt == 0)
    nrp=1;
  else
    nrp= nrpt;

  ix=1;
  if( n > 0)
  {
    i1= isegno( its, 1);
    if( i1 < 1)
      i1= 1;

    ix= i1;
    if( nrpt == 0)
      k= i1-1;
    else
    {
      k= n;
      /* Reallocate tags buffer */
      mreq = n+m + (n+1-i1)*nrpt;
      mem_realloc( (void *)&itag, mreq * sizeof(int) );

      /* Reallocate wire buffers */
      mreq = (n+(n+1-i1)*nrpt) * sizeof(long double);
      mem_realloc( (void *)&x, mreq );
      mem_realloc( (void *)&y, mreq );
      mem_realloc( (void *)&z, mreq );
      mem_realloc( (void *)&x2, mreq );
      mem_realloc( (void *)&y2, mreq );
      mem_realloc( (void *)&z2, mreq );
      mem_realloc( (void *)&bi, mreq );
    }

    for( ir = 0; ir < nrp; ir++ )
    {
      for( i = i1-1; i < n; i++ )
      {
	xi= x[i];
	yi= y[i];
	zi= z[i];
	x[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= x2[i];
	yi= y2[i];
	zi= z2[i];
	x2[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y2[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z2[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	bi[k]= bi[i];
	itag[k]= itag[i];
	if( itag[i] != 0)
	  itag[k]= itag[i]+ itgi;

	k++;

      } /* for( i = i1; i < n; i++ ) */

      i1= n+1;
      n= k;

    } /* for( ir = 0; ir < nrp; ir++ ) */

  } /* if( n >= n2) */

  if( m > 0)
  {
    i1 = 0;
    if( nrpt == 0)
      k= 0;
    else
      k = m;

    /* Reallocate patch buffers */
    mreq = m * (1+nrpt) * sizeof(long double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    for( ii = 0; ii < nrp; ii++ )
    {
      for( i = i1; i < m; i++ )
      {
	xi= px[i];
	yi= py[i];
	zi= pz[i];
	px[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	py[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	pz[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= t1x[i];
	yi= t1y[i];
	zi= t1z[i];
	t1x[k]= xi* xx+ yi* xy+ zi* xz;
	t1y[k]= xi* yx+ yi* yy+ zi* yz;
	t1z[k]= xi* zx+ yi* zy+ zi* zz;
	xi= t2x[i];
	yi= t2y[i];
	zi= t2z[i];
	t2x[k]= xi* xx+ yi* xy+ zi* xz;
	t2y[k]= xi* yx+ yi* yy+ zi* yz;
	t2z[k]= xi* zx+ yi* zy+ zi* zz;
	psalp[k]= psalp[i];
	pbi[k]= pbi[i];
	k++;

      } /* for( i = i1; i < m; i++ ) */

      i1= m;
      m = k;

    } /* for( ii = 0; ii < nrp; ii++ ) */

  } /* if( m >= m2) */

  if( (nrpt == 0) && (ix == 1) )
    return;

  np= n;
  mp= m;
  ipsym=0;

  return;
}

/*-----------------------------------------------------------------------*/

/* patch generates and modifies patch geometry data */
void patch( int nx, int ny,
    long double ax1, long double ay1, long double az1,
    long double ax2, long double ay2, long double az2,
    long double ax3, long double ay3, long double az3,
    long double ax4, long double ay4, long double az4 )
{
  int mi, ntp, iy, ix, mreq;
  long double s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., xst=0.;
  long double znv, xnv, ynv, xa, xn2, yn2, zn2, salpn, xs, ys, zs, xt, yt, zt;

  /* new patches.  for nx=0, ny=1,2,3,4 patch is (respectively) */;
  /* arbitrary, rectagular, triangular, or quadrilateral. */
  /* for nx and ny  > 0 a rectangular surface is produced with */
  /* nx by ny rectangular patches. */

  m++;
  mi= m-1;

  /* Reallocate patch buffers */
  mreq = m * sizeof(long double);
  mem_realloc( (void *)&px, mreq );
  mem_realloc( (void *)&py, mreq );
  mem_realloc( (void *)&pz, mreq );
  mem_realloc( (void *)&t1x, mreq );
  mem_realloc( (void *)&t1y, mreq );
  mem_realloc( (void *)&t1z, mreq );
  mem_realloc( (void *)&t2x, mreq );
  mem_realloc( (void *)&t2y, mreq );
  mem_realloc( (void *)&t2z, mreq );
  mem_realloc( (void *)&pbi, mreq );
  mem_realloc( (void *)&psalp, mreq );

  if( nx > 0)
    ntp=2;
  else
    ntp= ny;

  if( ntp <= 1)
  {
    px[mi]= ax1;
    py[mi]= ay1;
    pz[mi]= az1;
    pbi[mi]= az2;
    znv= cosl( ax2);
    xnv= znv* cosl( ay2);
    ynv= znv* sinl( ay2);
    znv= sinl( ax2);
    xa= sqrtl( xnv* xnv+ ynv* ynv);

    if( xa >= 1.0e-6)
    {
      t1x[mi]=- ynv/ xa;
      t1y[mi]= xnv/ xa;
      t1z[mi]=0.;
    }
    else
    {
      t1x[mi]=1.;
      t1y[mi]=0.;
      t1z[mi]=0.;
    }

  } /* if( ntp <= 1) */
  else
  {
    s1x= ax2- ax1;
    s1y= ay2- ay1;
    s1z= az2- az1;
    s2x= ax3- ax2;
    s2y= ay3- ay2;
    s2z= az3- az2;

    if( nx != 0)
    {
      s1x= s1x/ nx;
      s1y= s1y/ nx;
      s1z= s1z/ nx;
      s2x= s2x/ ny;
      s2y= s2y/ ny;
      s2z= s2z/ ny;
    }

    xnv= s1y* s2z- s1z* s2y;
    ynv= s1z* s2x- s1x* s2z;
    znv= s1x* s2y- s1y* s2x;
    xa= sqrtl( xnv* xnv+ ynv* ynv+ znv* znv);
    xnv= xnv/ xa;
    ynv= ynv/ xa;
    znv= znv/ xa;
    xst= sqrtl( s1x* s1x+ s1y* s1y+ s1z* s1z);
    t1x[mi]= s1x/ xst;
    t1y[mi]= s1y/ xst;
    t1z[mi]= s1z/ xst;

    if( ntp <= 2)
    {
      px[mi]= ax1+.5*( s1x+ s2x);
      py[mi]= ay1+.5*( s1y+ s2y);
      pz[mi]= az1+.5*( s1z+ s2z);
      pbi[mi]= xa;
    }
    else
    {
      if( ntp != 4)
      {
	px[mi]=( ax1+ ax2+ ax3)/3.;
	py[mi]=( ay1+ ay2+ ay3)/3.;
	pz[mi]=( az1+ az2+ az3)/3.;
	pbi[mi]=.5* xa;
      }
      else
      {
	s1x= ax3- ax1;
	s1y= ay3- ay1;
	s1z= az3- az1;
	s2x= ax4- ax1;
	s2y= ay4- ay1;
	s2z= az4- az1;
	xn2= s1y* s2z- s1z* s2y;
	yn2= s1z* s2x- s1x* s2z;
	zn2= s1x* s2y- s1y* s2x;
	xst= sqrtl( xn2* xn2+ yn2* yn2+ zn2* zn2);
	salpn=1./(3.*( xa+ xst));
	px[mi]=( xa*( ax1+ ax2+ ax3)+ xst*( ax1+ ax3+ ax4))* salpn;
	py[mi]=( xa*( ay1+ ay2+ ay3)+ xst*( ay1+ ay3+ ay4))* salpn;
	pz[mi]=( xa*( az1+ az2+ az3)+ xst*( az1+ az3+ az4))* salpn;
	pbi[mi]=.5*( xa+ xst);
	s1x=( xnv* xn2+ ynv* yn2+ znv* zn2)/ xst;

	if( s1x <= 0.9998)
	{
	  fprintf( output_fp,
	      "\n  ERROR -- CORNERS OF QUADRILATERAL"
	      " PATCH DO NOT LIE IN A PLANE" );
	  stop(-1);
	}

      } /* if( ntp != 4) */

    } /* if( ntp <= 2) */

  } /* if( ntp <= 1) */

  t2x[mi]= ynv* t1z[mi]- znv* t1y[mi];
  t2y[mi]= znv* t1x[mi]- xnv* t1z[mi];
  t2z[mi]= xnv* t1y[mi]- ynv* t1x[mi];
  psalp[mi]=1.;

  if( nx != 0)
  {
    m += nx*ny-1;

    /* Reallocate patch buffers */
    mreq = m * sizeof(long double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    xn2= px[mi]- s1x- s2x;
    yn2= py[mi]- s1y- s2y;
    zn2= pz[mi]- s1z- s2z;
    xs= t1x[mi];
    ys= t1y[mi];
    zs= t1z[mi];
    xt= t2x[mi];
    yt= t2y[mi];
    zt= t2z[mi];

    for( iy = 0; iy < ny; iy++ )
    {
      xn2 += s2x;
      yn2 += s2y;
      zn2 += s2z;

      for( ix = 1; ix <= nx; ix++ )
      {
	xst= (long double)ix;
	px[mi]= xn2+ xst* s1x;
	py[mi]= yn2+ xst* s1y;
	pz[mi]= zn2+ xst* s1z;
	pbi[mi]= xa;
	psalp[mi]=1.;
	t1x[mi]= xs;
	t1y[mi]= ys;
	t1z[mi]= zs;
	t2x[mi]= xt;
	t2y[mi]= yt;
	t2z[mi]= zt;
	mi++;
      } /* for( ix = 0; ix < nx; ix++ ) */

    } /* for( iy = 0; iy < ny; iy++ ) */

  } /* if( nx != 0) */

  ipsym=0;
  np= n;
  mp= m;

  return;
}

/*-----------------------------------------------------------------------*/

/*** this function was an 'entry point' (part of) 'patch()' ***/
void subph( int nx, int ny )
{
  int mia, ix, iy, mi, mreq;
  long double xs, ys, zs, xa, xst, s1x, s1y, s1z, s2x, s2y, s2z, saln, xt, yt;

    /* Reallocate patch buffers */
    if( ny == 0 )
      m += 3;
    else
      m += 4;

    mreq = m * sizeof(long double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

  /* Shift patches to make room for new ones */
  if( (ny == 0) && (nx != m) )
  {
    for( iy = m-1; iy > nx+2; iy-- )
    {
      ix = iy-3;
      px[iy]= px[ix];
      py[iy]= py[ix];
      pz[iy]= pz[ix];
      pbi[iy]= pbi[ix];
      psalp[iy]= psalp[ix];
      t1x[iy]= t1x[ix];
      t1y[iy]= t1y[ix];
      t1z[iy]= t1z[ix];
      t2x[iy]= t2x[ix];
      t2y[iy]= t2y[ix];
      t2z[iy]= t2z[ix];
    }

  } /* if( (ny == 0) || (nx != m) ) */

  /* divide patch for connection */
  mi= nx-1;
  xs= px[mi];
  ys= py[mi];
  zs= pz[mi];
  xa= pbi[mi]/4.;
  xst= sqrtl( xa)/2.;
  s1x= t1x[mi];
  s1y= t1y[mi];
  s1z= t1z[mi];
  s2x= t2x[mi];
  s2y= t2y[mi];
  s2z= t2z[mi];
  saln= psalp[mi];
  xt= xst;
  yt= xst;

  if( ny == 0)
    mia= mi;
  else
  {
    mp++;
    mia= m-1;
  }

  for( ix = 1; ix <= 4; ix++ )
  {
    px[mia]= xs+ xt* s1x+ yt* s2x;
    py[mia]= ys+ xt* s1y+ yt* s2y;
    pz[mia]= zs+ xt* s1z+ yt* s2z;
    pbi[mia]= xa;
    t1x[mia]= s1x;
    t1y[mia]= s1y;
    t1z[mia]= s1z;
    t2x[mia]= s2x;
    t2y[mia]= s2y;
    t2z[mia]= s2z;
    psalp[mia]= saln;

    if( ix == 2)
      yt=- yt;

    if( (ix == 1) || (ix == 3) )
      xt=- xt;

    mia++;
  }

  if( nx <= mp)
    mp += 3;

  if( ny > 0 )
    pz[mi]=10000.;

  return;
}

/*-----------------------------------------------------------------------*/

void readgm( char *gm, int *i1, int *i2, long double *x1, long double *y1,
    long double *z1, long double *x2, long double *y2, long double *z2, long double *rad )
{
  char line_buf[134];
  int nlin, i, line_idx;
  int nint = 2, nflt = 7;
  int iarr[2] = { 0, 0 };
  long double rarr[7] = { 0., 0., 0., 0., 0., 0., 0. };


  /* read a line from input file */
  load_line( line_buf, input_fp );

  /* get line length */
  nlin= strlen( line_buf );

  /* abort if card's mnemonic too short or missing */
  if( nlin < 2 )
  {
    fprintf( output_fp,
	"\n  GEOMETRY DATA CARD ERROR:"
	"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
    stop(-1);
  }

  /* extract card's mnemonic code */
  strncpy( gm, line_buf, 2 );
  gm[2] = '\0';

  /* Exit if "XT" command read (for testing) */
  if( strcmp( gm, "XT" ) == 0 )
  {
    fprintf( stderr,
	"\nnec2c: Exiting after an \"XT\" command in readgm()\n" );
    fprintf( output_fp,
	"\n\n  nec2c: Exiting after an \"XT\" command in readgm()" );
    stop(0);
  }

  /* Return if only mnemonic on card */
  if( nlin == 2 )
  {
    *i1 = *i2 = 0;
    *x1 = *y1 = *z1 = *x2 = *y2 = *z2 = *rad = 0.;
    return;
  }

  /* read integers from line */
  line_idx = 1;
  for( i = 0; i < nint; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*x1= rarr[0];
	*y1= rarr[1];
	*z1= rarr[2];
	*x2= rarr[3];
	*y2= rarr[4];
	*z2= rarr[5];
	*rad= rarr[6];
	return;
      }

    /* read an integer from line */
    iarr[i] = atoi( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while(
	(line_buf[++line_idx] != ' ') &&
	(line_buf[  line_idx] != '	') &&
	(line_buf[  line_idx] != ',') &&
	(line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-') )
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1)  );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *x1= rarr[0];
      *y1= rarr[1];
      *z1= rarr[2];
      *x2= rarr[3];
      *y2= rarr[4];
      *z2= rarr[5];
      *rad= rarr[6];
      return;
    }

  } /* for( i = 0; i < nint; i++ ) */

  /* read long doubles from line */
  for( i = 0; i < nflt; i++ )
  {
    /* Find first numerical character */
    while( ((line_buf[++line_idx] <  '0')  ||
	    (line_buf[  line_idx] >  '9')) &&
	    (line_buf[  line_idx] != '+')  &&
	    (line_buf[  line_idx] != '-')  &&
	    (line_buf[  line_idx] != '.') )
      if( (line_buf[line_idx] == '\0') )
      {
	*i1= iarr[0];
	*i2= iarr[1];
	*x1= rarr[0];
	*y1= rarr[1];
	*z1= rarr[2];
	*x2= rarr[3];
	*y2= rarr[4];
	*z2= rarr[5];
	*rad= rarr[6];
	return;
      }

    /* read a long double from line */
    rarr[i] = atof( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' or '\0' */
    line_idx--;
    while( 
	(line_buf[++line_idx] != ' ') &&
	(line_buf[  line_idx] != '	') &&
	(line_buf[  line_idx] != ',') &&
	(line_buf[  line_idx] != '\0') )
    {
      /* test for non-numerical characters */
      if( ((line_buf[line_idx] <  '0')  ||
	   (line_buf[line_idx] >  '9')) &&
	   (line_buf[line_idx] != '.')  &&
	   (line_buf[line_idx] != '+')  &&
	   (line_buf[line_idx] != '-')  &&
	   (line_buf[line_idx] != 'E')  &&
	   (line_buf[line_idx] != 'e') )
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d.\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *x1= rarr[0];
      *y1= rarr[1];
      *z1= rarr[2];
      *x2= rarr[3];
      *y2= rarr[4];
      *z2= rarr[5];
      *rad= rarr[6];
      return;
    }

  } /* for( i = 0; i < nflt; i++ ) */

  *i1  = iarr[0];
  *i2  = iarr[1];
  *x1  = rarr[0];
  *y1  = rarr[1];
  *z1  = rarr[2];
  *x2  = rarr[3];
  *y2  = rarr[4];
  *z2  = rarr[5];
  *rad = rarr[6];

  return;
}

/*-----------------------------------------------------------------------*/

/* reflc reflects partial structure along x,y, or z axes or rotates */
/* structure to complete a symmetric structure. */
void reflc( int ix, int iy, int iz, int itx, int nop )
{
  int iti, i, nx, itagi, k, mreq;
  long double e1, e2, fnop, sam, cs, ss, xk, yk;

  np= n;
  mp= m;
  ipsym=0;
  iti= itx;

  if( ix >= 0)
  {
    if( nop == 0)
      return;

    ipsym=1;

    /* reflect along z axis */
    if( iz != 0)
    {
      ipsym=2;

      if( n > 0 )
      {
	/* Reallocate tags buffer */
	mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );

	/* Reallocate wire buffers */
	mreq = 2*n * sizeof(long double);
	mem_realloc( (void *)&x, mreq );
	mem_realloc( (void *)&y, mreq );
	mem_realloc( (void *)&z, mreq );
	mem_realloc( (void *)&x2, mreq );
	mem_realloc( (void *)&y2, mreq );
	mem_realloc( (void *)&z2, mreq );
	mem_realloc( (void *)&bi, mreq );

	for( i = 0; i < n; i++ )
	{
	  nx= i+ n;
	  e1= z[i];
	  e2= z2[i];

	  if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--SEGMENT %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  x[nx]= x[i];
	  y[nx]= y[i];
	  z[nx]=- e1;
	  x2[nx]= x2[i];
	  y2[nx]= y2[i];
	  z2[nx]=- e2;
	  itagi= itag[i];

	  if( itagi == 0)
	    itag[nx]=0;
	  if( itagi != 0)
	    itag[nx]= itagi+ iti;

	  bi[nx]= bi[i];

	} /* for( i = 0; i < n; i++ ) */

	n= n*2;
	iti= iti*2;

      } /* if( n > 0) */

      if( m > 0 )
      {
	/* Reallocate patch buffers */
	mreq = 2*m * sizeof(long double);
	mem_realloc( (void *)&px, mreq );
	mem_realloc( (void *)&py, mreq );
	mem_realloc( (void *)&pz, mreq );
	mem_realloc( (void *)&t1x, mreq );
	mem_realloc( (void *)&t1y, mreq );
	mem_realloc( (void *)&t1z, mreq );
	mem_realloc( (void *)&t2x, mreq );
	mem_realloc( (void *)&t2y, mreq );
	mem_realloc( (void *)&t2z, mreq );
	mem_realloc( (void *)&pbi, mreq );
	mem_realloc( (void *)&psalp, mreq );

	for( i = 0; i < m; i++ )
	{
	  nx = i+m;
	  if( fabsl(pz[i]) <= 1.0e-10)
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--PATCH %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  px[nx]= px[i];
	  py[nx]= py[i];
	  pz[nx]=- pz[i];
	  t1x[nx]= t1x[i];
	  t1y[nx]= t1y[i];
	  t1z[nx]=- t1z[i];
	  t2x[nx]= t2x[i];
	  t2y[nx]= t2y[i];
	  t2z[nx]=- t2z[i];
	  psalp[nx]=- psalp[i];
	  pbi[nx]= pbi[i];
	}

	m= m*2;

      } /* if( m >= m2) */

    } /* if( iz != 0) */

    /* reflect along y axis */
    if( iy != 0)
    {
      if( n > 0)
      {
	/* Reallocate tags buffer */
	mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );/*????*/

	/* Reallocate wire buffers */
	mreq = 2*n * sizeof(long double);
	mem_realloc( (void *)&x, mreq );
	mem_realloc( (void *)&y, mreq );
	mem_realloc( (void *)&z, mreq );
	mem_realloc( (void *)&x2, mreq );
	mem_realloc( (void *)&y2, mreq );
	mem_realloc( (void *)&z2, mreq );
	mem_realloc( (void *)&bi, mreq );

	for( i = 0; i < n; i++ )
	{
	  nx= i+ n;
	  e1= y[i];
	  e2= y2[i];

	  if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--SEGMENT %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  x[nx]= x[i];
	  y[nx]=- e1;
	  z[nx]= z[i];
	  x2[nx]= x2[i];
	  y2[nx]=- e2;
	  z2[nx]= z2[i];
	  itagi= itag[i];

	  if( itagi == 0)
	    itag[nx]=0;
	  if( itagi != 0)
	    itag[nx]= itagi+ iti;

	  bi[nx]= bi[i];

	} /* for( i = n2-1; i < n; i++ ) */

	n= n*2;
	iti= iti*2;

      } /* if( n >= n2) */

      if( m > 0 )
      {
	/* Reallocate patch buffers */
	mreq = 2*m * sizeof(long double);
	mem_realloc( (void *)&px, mreq );
	mem_realloc( (void *)&py, mreq );
	mem_realloc( (void *)&pz, mreq );
	mem_realloc( (void *)&t1x, mreq );
	mem_realloc( (void *)&t1y, mreq );
	mem_realloc( (void *)&t1z, mreq );
	mem_realloc( (void *)&t2x, mreq );
	mem_realloc( (void *)&t2y, mreq );
	mem_realloc( (void *)&t2z, mreq );
	mem_realloc( (void *)&pbi, mreq );
	mem_realloc( (void *)&psalp, mreq );

	for( i = 0; i < m; i++ )
	{
	  nx= i+m;
	  if( fabsl( py[i]) <= 1.0e-10)
	  {
	    fprintf( output_fp,
		"\n  GEOMETRY DATA ERROR--PATCH %d"
		" LIES IN PLANE OF SYMMETRY", i+1 );
	    stop(-1);
	  }

	  px[nx]= px[i];
	  py[nx]=- py[i];
	  pz[nx]= pz[i];
	  t1x[nx]= t1x[i];
	  t1y[nx]=- t1y[i];
	  t1z[nx]= t1z[i];
	  t2x[nx]= t2x[i];
	  t2y[nx]=- t2y[i];
	  t2z[nx]= t2z[i];
	  psalp[nx]=- psalp[i];
	  pbi[nx]= pbi[i];

	} /* for( i = m2; i <= m; i++ ) */

	m= m*2;

      } /* if( m >= m2) */

    } /* if( iy != 0) */

    /* reflect along x axis */
    if( ix == 0 )
      return;

    if( n > 0 )
    {
      /* Reallocate tags buffer */
      mem_realloc( (void *)&itag, (2*n+m) * sizeof(int) );/*????*/

      /* Reallocate wire buffers */
      mreq = 2*n * sizeof(long double);
      mem_realloc( (void *)&x, mreq );
      mem_realloc( (void *)&y, mreq );
      mem_realloc( (void *)&z, mreq );
      mem_realloc( (void *)&x2, mreq );
      mem_realloc( (void *)&y2, mreq );
      mem_realloc( (void *)&z2, mreq );
      mem_realloc( (void *)&bi, mreq );

      for( i = 0; i < n; i++ )
      {
	nx= i+ n;
	e1= x[i];
	e2= x2[i];

	if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
	{
	  fprintf( output_fp,
	      "\n  GEOMETRY DATA ERROR--SEGMENT %d"
	      " LIES IN PLANE OF SYMMETRY", i+1 );
	  stop(-1);
	}

	x[nx]=- e1;
	y[nx]= y[i];
	z[nx]= z[i];
	x2[nx]=- e2;
	y2[nx]= y2[i];
	z2[nx]= z2[i];
	itagi= itag[i];

	if( itagi == 0)
	  itag[nx]=0;
	if( itagi != 0)
	  itag[nx]= itagi+ iti;

	bi[nx]= bi[i];
      }

      n= n*2;

    } /* if( n > 0) */

    if( m == 0 )
      return;

    /* Reallocate patch buffers */
    mreq = 2*m * sizeof(long double);
    mem_realloc( (void *)&px, mreq );
    mem_realloc( (void *)&py, mreq );
    mem_realloc( (void *)&pz, mreq );
    mem_realloc( (void *)&t1x, mreq );
    mem_realloc( (void *)&t1y, mreq );
    mem_realloc( (void *)&t1z, mreq );
    mem_realloc( (void *)&t2x, mreq );
    mem_realloc( (void *)&t2y, mreq );
    mem_realloc( (void *)&t2z, mreq );
    mem_realloc( (void *)&pbi, mreq );
    mem_realloc( (void *)&psalp, mreq );

    for( i = 0; i < m; i++ )
    {
      nx= i+m;
      if( fabsl( px[i]) <= 1.0e-10)
      {
	fprintf( output_fp,
	    "\n  GEOMETRY DATA ERROR--PATCH %d"
	    " LIES IN PLANE OF SYMMETRY", i+1 );
	stop(-1);
      }

      px[nx]=- px[i];
      py[nx]= py[i];
      pz[nx]= pz[i];
      t1x[nx]=- t1x[i];
      t1y[nx]= t1y[i];
      t1z[nx]= t1z[i];
      t2x[nx]=- t2x[i];
      t2y[nx]= t2y[i];
      t2z[nx]= t2z[i];
      psalp[nx]=- psalp[i];
      pbi[nx]= pbi[i];
    }

    m= m*2;
    return;

  } /* if( ix >= 0) */

  /* reproduce structure with rotation to form cylindrical structure */
  fnop= (long double)nop;
  ipsym=-1;
  sam=TP/ fnop;
  cs= cosl( sam);
  ss= sinl( sam);

  if( n > 0)
  {
    n *= nop;
    nx= np;

    /* Reallocate tags buffer */
    mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

    /* Reallocate wire buffers */
    mreq = n * sizeof(long double);
    mem_realloc( (void *)&x, mreq );
    mem_realloc( (void *)&y, mreq );
    mem_realloc( (void *)&z, mreq );
    mem_realloc( (void *)&x2, mreq );
    mem_realloc( (void *)&y2, mreq );
    mem_realloc( (void *)&z2, mreq );
    mem_realloc( (void *)&bi, mreq );

    for( i = nx; i < n; i++ )
    {
      k= i- np;
      xk= x[k];
      yk= y[k];
      x[i]= xk* cs- yk* ss;
      y[i]= xk* ss+ yk* cs;
      z[i]= z[k];
      xk= x2[k];
      yk= y2[k];
      x2[i]= xk* cs- yk* ss;
      y2[i]= xk* ss+ yk* cs;
      z2[i]= z2[k];
      bi[i]= bi[k];
      itagi= itag[k];

      if( itagi == 0)
	itag[i]=0;
      if( itagi != 0)
	itag[i]= itagi+ iti;
    }

  } /* if( n >= n2) */

  if( m == 0 )
    return;

  m *= nop;
  nx= mp;

  /* Reallocate patch buffers */
  mreq = m * sizeof(long double);
  mem_realloc( (void *)&px, mreq  );
  mem_realloc( (void *)&py, mreq  );
  mem_realloc( (void *)&pz, mreq );
  mem_realloc( (void *)&t1x, mreq );
  mem_realloc( (void *)&t1y, mreq );
  mem_realloc( (void *)&t1z, mreq );
  mem_realloc( (void *)&t2x, mreq );
  mem_realloc( (void *)&t2y, mreq );
  mem_realloc( (void *)&t2z, mreq );
  mem_realloc( (void *)&pbi, mreq );
  mem_realloc( (void *)&psalp, mreq );

  for( i = nx; i < m; i++ )
  {
    k = i-mp;
    xk= px[k];
    yk= py[k];
    px[i]= xk* cs- yk* ss;
    py[i]= xk* ss+ yk* cs;
    pz[i]= pz[k];
    xk= t1x[k];
    yk= t1y[k];
    t1x[i]= xk* cs- yk* ss;
    t1y[i]= xk* ss+ yk* cs;
    t1z[i]= t1z[k];
    xk= t2x[k];
    yk= t2y[k];
    t2x[i]= xk* cs- yk* ss;
    t2y[i]= xk* ss+ yk* cs;
    t2z[i]= t2z[k];
    psalp[i]= psalp[k];
    pbi[i]= pbi[k];

  } /* for( i = nx; i < m; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine wire generates segment geometry */
/* data for a straight wire of ns segments. */
void wire( long double xw1, long double yw1, long double zw1,
    long double xw2, long double yw2, long double zw2, long double rad,
    long double rdel, long double rrad, int ns, int itg )
{
  int ist, i, mreq;
  long double xd, yd, zd, delz, rd, fns, radz;
  long double xs1, ys1, zs1, xs2, ys2, zs2;

  ist= n;
  n= n+ ns;
  np= n;
  mp= m;
  ipsym=0;

  if( ns < 1)
    return;

  /* Reallocate tags buffer */
  mem_realloc( (void *)&itag, (n+m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = n * sizeof(long double);
  mem_realloc( (void *)&x, mreq );
  mem_realloc( (void *)&y, mreq );
  mem_realloc( (void *)&z, mreq );
  mem_realloc( (void *)&x2, mreq );
  mem_realloc( (void *)&y2, mreq );
  mem_realloc( (void *)&z2, mreq );
  mem_realloc( (void *)&bi, mreq );

  xd= xw2- xw1;
  yd= yw2- yw1;
  zd= zw2- zw1;

  if( fabsl( rdel-1.) >= 1.0e-6)
  {
    delz= sqrtl( xd* xd+ yd* yd+ zd* zd);
    xd= xd/ delz;
    yd= yd/ delz;
    zd= zd/ delz;
    delz= delz*(1.- rdel)/(1.- powl(rdel, ns) );
    rd= rdel;
  }
  else
  {
    fns= ns;
    xd= xd/ fns;
    yd= yd/ fns;
    zd= zd/ fns;
    delz=1.;
    rd=1.;
  }

  radz= rad;
  xs1= xw1;
  ys1= yw1;
  zs1= zw1;

  for( i = ist; i < n; i++ )
  {
    itag[i]= itg;
    xs2= xs1+ xd* delz;
    ys2= ys1+ yd* delz;
    zs2= zs1+ zd* delz;
    x[i]= xs1;
    y[i]= ys1;
    z[i]= zs1;
    x2[i]= xs2;
    y2[i]= ys2;
    z2[i]= zs2;
    bi[i]= radz;
    delz= delz* rd;
    radz= radz* rrad;
    xs1= xs2;
    ys1= ys2;
    zs1= zs2;
  }

  x2[n-1]= xw2;
  y2[n-1]= yw2;
  z2[n-1]= zw2;

  return;
}

/*-----------------------------------------------------------------------*/


