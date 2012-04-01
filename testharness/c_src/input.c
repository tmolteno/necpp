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

/* common  /vsorc/ */
extern int *ivqd, *isant, *iqds, nvqd, nsant, nqds;
extern complex long double *vqd, *vqds, *vsant;

/* common  /dataj/ */
extern int iexk, ind1, indd1, ind2, indd2, ipgnd;
extern long double s, b, xj, yj, zj, cabj, sabj, salpj, rkh;
extern long double t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
extern complex long double  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;

/* common  /zload/ */
extern int nload;
extern complex long double *zarray;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/*-------------------------------------------------------------------*/

/* fill incident field array for charge discontinuity voltage source */
void qdsrc( int is, complex long double v, complex long double *e )
{
  int i, jx, j, jp1, ipr, ij, i1;
  long double xi, yi, zi, ai, cabi, sabi, salpi, tx, ty, tz;
  complex long double curd, etk, ets, etc;

  is--;
  i= icon1[is];
  icon1[is]=0;
  tbf( is+1,0);
  icon1[is]= i;
  s= si[is]*.5;
  curd= CCJ* v/(( logl(2.* s/ bi[is])-1.)*( bx[jsno-1]*
	cosl( TP* s)+ cx[jsno-1]* sinl( TP* s))* wlam);
  vqds[nqds]= v;
  iqds[nqds]= is+1;
  nqds++;

  for( jx = 0; jx < jsno; jx++ )
  {
    j= jco[jx]-1;
    jp1 = j+1;
    s= si[j];
    b= bi[j];
    xj= x[j];
    yj= y[j];
    zj= z[j];
    cabj= cab[j];
    sabj= sab[j];
    salpj= salp[j];

    if( iexk != 0)
    {
      ipr= icon1[j];

      if( ipr < 0 )
      {
	ipr=- ipr;
	ipr--;
	if( -icon1[ipr-1] != jp1 )
	  ind1=2;
	else
	{
	  xi= fabsl( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	  if( (xi < 0.999999) || (fabsl(bi[ipr]/b-1.) > 1.0e-6) )
	    ind1=2;
	  else
	    ind1=0;
	}
      }  /* if( ipr < 0 ) */
      else
	if( ipr == 0 )
	  ind1=1;
	else /* ipr > 0 */
	{
	  ipr--;
	  if( ipr != j )
	  {
	    if( icon2[ipr] != jp1)
	      ind1=2;
	    else
	    {
	      xi= fabsl( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	      if( (xi < 0.999999) || (fabsl(bi[ipr]/b-1.) > 1.0e-6) )
		ind1=2;
	      else
		ind1=0;
	    }
	  } /* if( ipr != j ) */
	  else
	  {
	    if( cabj* cabj+ sabj* sabj > 1.0e-8)
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* else */

      ipr= icon2[j];
      if( ipr < 0 )
      {
	ipr = -ipr;
	ipr--;
	if( -icon2[ipr] != jp1 )
	  ind1=2;
	else
	{
	  xi= fabsl( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	  if( (xi < 0.999999) || (fabsl(bi[ipr]/b-1.) > 1.0e-6) )
	    ind1=2;
	  else
	    ind1=0;
	}
      } /* if( ipr < 0 ) */
      else
	if( ipr == 0 )
	  ind2=1;
	else /* ipr > 0 */
	{
	  ipr--;
	  if( ipr != j )
	  {
	    if( icon1[ipr] != jp1)
	      ind2=2;
	    else
	    {
	      xi= fabsl( cabj* cab[ipr]+ sabj* sab[ipr]+ salpj* salp[ipr]);
	      if( (xi < 0.999999) || (fabsl(bi[ipr]/b-1.) > 1.0e-6) )
		ind2=2;
	      else
		ind2=0;
	    }
	  } /* if( ipr != j )*/
	  else
	  {
	    if( cabj* cabj+ sabj* sabj > 1.0e-8)
	      ind1=2;
	    else
	      ind1=0;
	  }
	} /* else */

    } /* if( iexk != 0) */

    for( i = 0; i < n; i++ )
    {
      ij= i- j;
      xi= x[i];
      yi= y[i];
      zi= z[i];
      ai= bi[i];
      efld( xi, yi, zi, ai, ij);
      cabi= cab[i];
      sabi= sab[i];
      salpi= salp[i];
      etk= exk* cabi+ eyk* sabi+ ezk* salpi;
      ets= exs* cabi+ eys* sabi+ ezs* salpi;
      etc= exc* cabi+ eyc* sabi+ ezc* salpi;
      e[i]= e[i]-( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx])* curd;
    }

    if( m != 0)
    {
      i1= n-1;
      for( i = 0; i < m; i++ )
      {
	xi= px[i];
	yi= py[i];
	zi= pz[i];
	hsfld( xi, yi, zi,0.);
	i1++;
	tx= t2x[i];
	ty= t2y[i];
	tz= t2z[i];
	etk= exk* tx+ eyk* ty+ ezk* tz;
	ets= exs* tx+ eys* ty+ ezs* tz;
	etc= exc* tx+ eyc* ty+ ezc* tz;
	e[i1] += ( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx] )* curd* psalp[i];
	i1++;
	tx= t1x[i];
	ty= t1y[i];
	tz= t1z[i];
	etk= exk* tx+ eyk* ty+ ezk* tz;
	ets= exs* tx+ eys* ty+ ezs* tz;
	etc= exc* tx+ eyc* ty+ ezc* tz;
	e[i1] += ( etk* ax[jx]+ ets* bx[jx]+ etc* cx[jx])* curd* psalp[i];
      }

    } /* if( m != 0) */

    if( nload > 0 )
      e[j] += zarray[j]* curd*(ax[jx]+ cx[jx]);

  } /* for( jx = 0; jx < jsno; jx++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

void readmn( char *gm, int *i1, int *i2, int *i3, int *i4,
    long double *f1, long double *f2, long double *f3,
    long double *f4, long double *f5, long double *f6 )
{
  char line_buf[134];
  int nlin, i, line_idx;
  int nint = 4, nflt = 6;
  int iarr[4] = { 0, 0, 0, 0 };
  long double rarr[6] = { 0., 0., 0., 0., 0., 0. };

  /* read a line from input file */
  load_line( line_buf, input_fp );

  /* get line length */
  nlin= strlen( line_buf );

  /* abort if card's mnemonic too short or missing */
  if( nlin < 2 )
  {
    fprintf( output_fp,
	"\n  COMMAND DATA CARD ERROR:"
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
    *i1 = *i2 = *i3 = *i4 = 0;
    *f1 = *f2 = *f3 = *f4 = *f5 = *f6 = 0.0;
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
	*i3= iarr[2];
	*i4= iarr[3];
	*f1= rarr[0];
	*f2= rarr[1];
	*f3= rarr[2];
	*f4= rarr[3];
	*f5= rarr[4];
	*f6= rarr[5];
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
	    "\n  COMMAND DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *i3= iarr[2];
      *i4= iarr[3];
      *f1= rarr[0];
      *f2= rarr[1];
      *f3= rarr[2];
      *f4= rarr[3];
      *f5= rarr[4];
      *f6= rarr[5];
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
	*i3= iarr[2];
	*i4= iarr[3];
	*f1= rarr[0];
	*f2= rarr[1];
	*f3= rarr[2];
	*f4= rarr[3];
	*f5= rarr[4];
	*f6= rarr[5];
	return;
      }

    /* read a long double from line */
    rarr[i] = atof( &line_buf[line_idx] );

    /* traverse numerical field to next ' ' or ',' */
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
	    "\n  COMMAND DATA CARD \"%s\" ERROR:"
	    "\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d\n",
	    gm, line_buf[line_idx], (line_idx+1) );
	stop(-1);
      }

    } /* while( (line_buff[++line_idx] ... */

    /* Return on end of line */
    if( line_buf[line_idx] == '\0' )
    {
      *i1= iarr[0];
      *i2= iarr[1];
      *i3= iarr[2];
      *i4= iarr[3];
      *f1= rarr[0];
      *f2= rarr[1];
      *f3= rarr[2];
      *f4= rarr[3];
      *f5= rarr[4];
      *f6= rarr[5];
      return;
    }

  } /* for( i = 0; i < nflt; i++ ) */

  *i1= iarr[0];
  *i2= iarr[1];
  *i3= iarr[2];
  *i4= iarr[3];
  *f1= rarr[0];
  *f2= rarr[1];
  *f3= rarr[2];
  *f4= rarr[3];
  *f5= rarr[4];
  *f6= rarr[5];

  return;
}

/*-----------------------------------------------------------------------*/

