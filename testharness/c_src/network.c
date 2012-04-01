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

/* common  /netcx/ */
extern int masym, neq, npeq, neq2, nonet, ntsol, nprint;
extern int *iseg1, *iseg2, *ntyp;
extern long double *x11r, *x11i, *x12r;
extern long double *x12i, *x22r, *x22i;
extern long double pin, pnls;
extern complex long double zped;

/* common  /vsorc/ */
extern int *ivqd, *isant, *iqds, nvqd, nsant, nqds;
extern complex long double *vqd, *vqds, *vsant;

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

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/*-------------------------------------------------------------------*/


/* subroutine netwk solves for structure currents for a given */
/* excitation including the effect of non-radiating networks if */
/* present. */
void netwk( complex long double *cm, complex long double *cmb,
    complex long double *cmc, complex long double *cmd, int *ip,
    complex long double *einc )
{
  int *ipnt = NULL, *nteqa = NULL, *ntsca = NULL;
  int jump1, jump2, nteq=0, ntsc=0, nseg2, irow2=0, j, ndimn;
  int neqz2, neqt, irow1=0, i, nseg1, isc1=0, isc2=0;
  long double asmx, asa, pwr, y11r, y11i, y12r, y12i, y22r, y22i;
  complex long double *vsrc = NULL, *rhs = NULL, *cmn = NULL;
  complex long double *rhnt = NULL, *rhnx = NULL, ymit, vlt, cux;

  neqz2= neq2;
  if( neqz2 == 0)
    neqz2=1;

  pin=0.;
  pnls=0.;
  neqt= neq+ neq2;
  ndimn = j = (2*nonet + nsant);

  /* Allocate network buffers */
  if( nonet > 0 )
  {
    mem_alloc( (void *)&rhs, np3m * sizeof(complex long double) );

    i = j * sizeof(complex long double);
    mem_alloc( (void *)&rhnt, i );
    mem_alloc( (void *)&rhnx, i );
    mem_alloc( (void *)&cmn, i * j );

    i = j * sizeof(int);
    mem_alloc( (void *)&ntsca, i );
    mem_alloc( (void *)&nteqa, i );
    mem_alloc( (void *)&ipnt, i );

    mem_alloc( (void *)&vsrc, nsant * sizeof(complex long double) );
  }

  if( ntsol == 0)
  {
    /* compute relative matrix asymmetry */
    if( masym != 0)
    {
      irow1=0;
      if( nonet != 0)
      {
	for( i = 0; i < nonet; i++ )
	{
	  nseg1= iseg1[i];
	  for( isc1 = 0; isc1 < 2; isc1++ )
	  {
	    if( irow1 == 0)
	    {
	      ipnt[irow1]= nseg1;
	      nseg1= iseg2[i];
	      irow1++;
	      continue;
	    }

	    for( j = 0; j < irow1; j++ )
	      if( nseg1 == ipnt[j])
		break;

	    if( j == irow1 )
	    {
	      ipnt[irow1]= nseg1;
	      irow1++;
	    }

	    nseg1= iseg2[i];

	  } /* for( isc1 = 0; isc1 < 2; isc1++ ) */

	} /* for( i = 0; i < nonet; i++ ) */

      } /* if( nonet != 0) */

      if( nsant != 0)
      {
	for( i = 0; i < nsant; i++ )
	{
	  nseg1= isant[i];
	  if( irow1 == 0)
	  {
	    ipnt[irow1]= nseg1;
	    irow1++;
	    continue;
	  }

	  for( j = 0; j < irow1; j++ )
	    if( nseg1 == ipnt[j])
	      break;

	  if( j == irow1 )
	  {
	    ipnt[irow1]= nseg1;
	    irow1++;
	  }

	} /* for( i = 0; i < nsant; i++ ) */

      } /* if( nsant != 0) */

      if( irow1 >= 2)
      {
	for( i = 0; i < irow1; i++ )
	{
	  isc1= ipnt[i]-1;
	  asmx= si[isc1];

	  for( j = 0; j < neqt; j++ )
	    rhs[j] = CPLX_00;

	  rhs[isc1] = CPLX_10;
	  solves( cm, ip, rhs, neq, 1, np, n, mp, m);
	  cabc( rhs);

	  for( j = 0; j < irow1; j++ )
	  {
	    isc1= ipnt[j]-1;
	    cmn[j+i*ndimn]= rhs[isc1]/ asmx;
	  }

	} /* for( i = 0; i < irow1; i++ ) */

	asmx=0.;
	asa=0.;

	for( i = 1; i < irow1; i++ )
	{
	  isc1= i;
	  for( j = 0; j < isc1; j++ )
	  {
	    cux= cmn[i+j*ndimn];
	    pwr= cabsl(( cux- cmn[j+i*ndimn])/ cux);
	    asa += pwr* pwr;

	    if( pwr < asmx)
	      continue;

	    asmx= pwr;
	    nteq= ipnt[i];
	    ntsc= ipnt[j];

	  } /* for( j = 0; j < isc1; j++ ) */

	} /* for( i = 1; i < irow1; i++ ) */

	asa= sqrtl( asa*2./ (long double)( irow1*( irow1-1)));
	fprintf( output_fp, "\n\n"
	    "   MAXIMUM RELATIVE ASYMMETRY OF THE DRIVING POINT ADMITTANCE\n"
	    "   MATRIX IS %10.3LE FOR SEGMENTS %d AND %d\n"
	    "   RMS RELATIVE ASYMMETRY IS %10.3LE",
	    asmx, nteq, ntsc, asa );

      } /* if( irow1 >= 2) */

    } /* if( masym != 0) */

    /* solution of network equations */
    if( nonet != 0)
    {
      for( i = 0; i < ndimn; i++ )
      {
	rhnx[i]=CPLX_00;
	for( j = 0; j < ndimn; j++ )
	  cmn[j+i*ndimn]=CPLX_00;
      }

      /* sort network and source data and */
      /* assign equation numbers to segments */
      nteq=0;
      ntsc=0;

      for( j = 0; j < nonet; j++ )
      {
	nseg1= iseg1[j];
	nseg2= iseg2[j];

	if( ntyp[j] <= 1)
	{
	  y11r= x11r[j];
	  y11i= x11i[j];
	  y12r= x12r[j];
	  y12i= x12i[j];
	  y22r= x22r[j];
	  y22i= x22i[j];
	}
	else
	{
	  y22r= TP* x11i[j]/ wlam;
	  y12r=0.;
	  y12i=1./( x11r[j]* sinl( y22r));
	  y11r= x12r[j];
	  y11i=- y12i* cosl( y22r);
	  y22r= x22r[j];
	  y22i= y11i+ x22i[j];
	  y11i= y11i+ x12i[j];

	  if( ntyp[j] != 2)
	  {
	    y12r=- y12r;
	    y12i=- y12i;
	  }

	} /* if( ntyp[j] <= 1) */

	jump1 = FALSE;
	if( nsant != 0)
	{
	  for( i = 0; i < nsant; i++ )
	    if( nseg1 == isant[i])
	    {
	      isc1 = i;
	      jump1 = TRUE;
	      break;
	    }
	} /* if( nsant != 0) */

	jump2 = FALSE;
	if( ! jump1 )
	{
	  isc1=-1;

	  if( nteq != 0)
	  {
	    for( i = 0; i < nteq; i++ )
	      if( nseg1 == nteqa[i])
	      {
		irow1 = i;
		jump2 = TRUE;
		break;
	      }

	  } /* if( nteq != 0) */

	  if( ! jump2 )
	  {
	    irow1= nteq;
	    nteqa[nteq]= nseg1;
	    nteq++;
	  }

	} /* if( ! jump1 ) */
	else
	{
	  if( ntsc != 0)
	  {
	    for( i = 0; i < ntsc; i++ )
	    {
	      if( nseg1 == ntsca[i])
	      {
		irow1 = ndimn- (i+1);
		jump2 = TRUE;
		break;
	      }
	    }

	  } /* if( ntsc != 0) */

	  if( ! jump2 )
	  {
	    irow1= ndimn- (ntsc+1);
	    ntsca[ntsc]= nseg1;
	    vsrc[ntsc]= vsant[isc1];
	    ntsc++;
	  }

	} /* if( ! jump1 ) */

	jump1 = FALSE;
	if( nsant != 0)
	{
	  for( i = 0; i < nsant; i++ )
	  {
	    if( nseg2 == isant[i])
	    {
	      isc2= i;
	      jump1 = TRUE;
	      break;
	    }
	  }

	} /* if( nsant != 0) */

	jump2 = FALSE;
	if( ! jump1 )
	{
	  isc2=-1;

	  if( nteq != 0)
	  {
	    for( i = 0; i < nteq; i++ )
	      if( nseg2 == nteqa[i])
	      {
		irow2= i;
		jump2 = TRUE;
		break;
	      }

	  } /* if( nteq != 0) */

	  if( ! jump2 )
	  {
	    irow2= nteq;
	    nteqa[nteq]= nseg2;
	    nteq++;
	  }

	}  /* if( ! jump1 ) */
	else
	{
	  if( ntsc != 0)
	  {
	    for( i = 0; i < ntsc; i++ )
	      if( nseg2 == ntsca[i])
	      {
		irow2 = ndimn- (i+1);
		jump2 = TRUE;
		break;
	      }

	  } /* if( ntsc != 0) */

	  if( ! jump2 )
	  {
	    irow2= ndimn- (ntsc+1);
	    ntsca[ntsc]= nseg2;
	    vsrc[ntsc]= vsant[isc2];
	    ntsc++;
	  }

	} /* if( ! jump1 ) */

	/* fill network equation matrix and right hand side vector with */
	/* network short-circuit admittance matrix coefficients. */
	if( isc1 == -1)
	{
	  cmn[irow1+irow1*ndimn] -= cmplx( y11r, y11i)* si[nseg1-1];
	  cmn[irow1+irow2*ndimn] -= cmplx( y12r, y12i)* si[nseg1-1];
	}
	else
	{
	  rhnx[irow1] += cmplx( y11r, y11i)* vsant[isc1]/wlam;
	  rhnx[irow2] += cmplx( y12r, y12i)* vsant[isc1]/wlam;
	}

	if( isc2 == -1)
	{
	  cmn[irow2+irow2*ndimn] -= cmplx( y22r, y22i)* si[nseg2-1];
	  cmn[irow2+irow1*ndimn] -= cmplx( y12r, y12i)* si[nseg2-1];
	}
	else
	{
	  rhnx[irow1] += cmplx( y12r, y12i)* vsant[isc2]/wlam;
	  rhnx[irow2] += cmplx( y22r, y22i)* vsant[isc2]/wlam;
	}

      } /* for( j = 0; j < nonet; j++ ) */

      /* add interaction matrix admittance */
      /* elements to network equation matrix */
      for( i = 0; i < nteq; i++ )
      {
	for( j = 0; j < neqt; j++ )
	  rhs[j] = CPLX_00;

	irow1= nteqa[i]-1;
	rhs[irow1]=CPLX_10;
	solves( cm, ip, rhs, neq, 1, np, n, mp, m);
	cabc( rhs);

	for( j = 0; j < nteq; j++ )
	{
	  irow1= nteqa[j]-1;
	  cmn[i+j*ndimn] += rhs[irow1];
	}

      } /* for( i = 0; i < nteq; i++ ) */

      /* factor network equation matrix */
      factr( nteq, cmn, ipnt, ndimn);

    } /* if( nonet != 0) */

  } /* if( ntsol != 0) */

  if( nonet != 0)
  {
    /* add to network equation right hand side */
    /* the terms due to element interactions */
    for( i = 0; i < neqt; i++ )
      rhs[i]= einc[i];

    solves( cm, ip, rhs, neq, 1, np, n, mp, m);
    cabc( rhs);

    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      rhnt[i]= rhnx[i]+ rhs[irow1];
    }

    /* solve network equations */
    solve( nteq, cmn, ipnt, rhnt, ndimn);

    /* add fields due to network voltages to electric fields */
    /* applied to structure and solve for induced current */
    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      einc[irow1] -= rhnt[i];
    }

    solves( cm, ip, einc, neq, 1, np, n, mp, m);
    cabc( einc);

    if( nprint == 0)
    {
      fprintf( output_fp, "\n\n\n"
	  "                          "
	  "--------- STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS --------" );

      fprintf( output_fp, "\n"
	  "  TAG   SEG       VOLTAGE (VOLTS)          CURRENT (AMPS)        "
	  " IMPEDANCE (OHMS)       ADMITTANCE (MHOS)     POWER\n"
	  "  No:   No:     REAL      IMAGINARY     REAL      IMAGINARY    "
	  " REAL      IMAGINARY     REAL      IMAGINARY   (WATTS)" );
    }

    for( i = 0; i < nteq; i++ )
    {
      irow1= nteqa[i]-1;
      vlt= rhnt[i]* si[irow1]* wlam;
      cux= einc[irow1]* wlam;
      ymit= cux/ vlt;
      zped= vlt/ cux;
      irow2= itag[irow1];
      pwr=.5* creall( vlt* conjl( cux));
      pnls= pnls- pwr;

      if( nprint == 0)
	fprintf( output_fp, "\n"
	    " %4d %5d %11.4LE %11.4LE %11.4LE %11.4LE"
	    " %11.4LE %11.4LE %11.4LE %11.4LE %11.4LE",
	    irow2, irow1+1, creall(vlt), cimagl(vlt), creall(cux), cimagl(cux),
	    creall(zped), cimagl(zped), creall(ymit), cimagl(ymit), pwr );
    }

    if( ntsc != 0)
    {
      for( i = 0; i < ntsc; i++ )
      {
	irow1= ntsca[i]-1;
	vlt= vsrc[i];
	cux= einc[irow1]* wlam;
	ymit= cux/ vlt;
	zped= vlt/ cux;
	irow2= itag[irow1];
	pwr=.5* creall( vlt* conjl( cux));
	pnls= pnls- pwr;

	if( nprint == 0)
	  fprintf( output_fp, "\n"
	      " %4d %5d %11.4LE %11.4LE %11.4LE %11.4LE"
	      " %11.4LE %11.4LE %11.4LE %11.4LE %11.4LE",
	      irow2, irow1+1, creall(vlt), cimagl(vlt), creall(cux), cimagl(cux),
	      creall(zped), cimagl(zped), creall(ymit), cimagl(ymit), pwr );

      } /* for( i = 0; i < ntsc; i++ ) */

    } /* if( ntsc != 0) */

  } /* if( nonet != 0) */
  else
  {
    /* solve for currents when no networks are present */
    solves( cm, ip, einc, neq, 1, np, n, mp, m);
    cabc( einc);
    ntsc=0;
  }

  if( (nsant+nvqd) == 0)
    return;

  fprintf( output_fp, "\n\n\n"
      "                        "
      "--------- ANTENNA INPUT PARAMETERS ---------" );

  fprintf( output_fp, "\n"
      "  TAG   SEG       VOLTAGE (VOLTS)         "
      "CURRENT (AMPS)         IMPEDANCE (OHMS)    "
      "    ADMITTANCE (MHOS)     POWER\n"
      "  No:   No:     REAL      IMAGINARY"
      "     REAL      IMAGINARY     REAL      "
      "IMAGINARY    REAL       IMAGINARY   (WATTS)" );

  if( nsant != 0)
  {
    for( i = 0; i < nsant; i++ )
    {
      isc1= isant[i]-1;
      vlt= vsant[i];

      if( ntsc == 0)
      {
	cux= einc[isc1]* wlam;
	irow1=0;
      }
      else
      {
	for( j = 0; j < ntsc; j++ )
	  if( ntsca[j] == isc1+1)
	    break;

	irow1= ndimn- (j+1);
	cux= rhnx[irow1];
	for( j = 0; j < nteq; j++ )
	  cux -= cmn[j+irow1*ndimn]*rhnt[j];
	cux=(einc[isc1]+ cux)* wlam;
	irow1++;

      } /* if( ntsc == 0) */

      ymit= cux/ vlt;
      zped= vlt/ cux;
      pwr=.5* creall( vlt* conjl( cux));
      pin= pin+ pwr;

      if( irow1 != 0)
	pnls= pnls+ pwr;

      irow2= itag[isc1];
      fprintf( output_fp, "\n"
	  " %4d %5d %11.4LE %11.4LE %11.4LE %11.4LE"
	  " %11.4LE %11.4LE %11.4LE %11.4LE %11.4LE",
	  irow2, isc1+1, creall(vlt), cimagl(vlt), creall(cux), cimagl(cux),
	  creall(zped), cimagl(zped), creall(ymit), cimagl(ymit), pwr );

    } /* for( i = 0; i < nsant; i++ ) */

  } /* if( nsant != 0) */

  if( nvqd == 0)
    return;

  for( i = 0; i < nvqd; i++ )
  {
    isc1= ivqd[i]-1;
    vlt= vqd[i];
    cux= cmplx( air[isc1], aii[isc1]);
    ymit= cmplx( bir[isc1], bii[isc1]);
    zped= cmplx( cir[isc1], cii[isc1]);
    pwr= si[isc1]* TP*.5;
    cux=( cux- ymit* sinl( pwr)+ zped* cosl( pwr))* wlam;
    ymit= cux/ vlt;
    zped= vlt/ cux;
    pwr=.5* creall( vlt* conjl( cux));
    pin= pin+ pwr;
    irow2= itag[isc1];

    fprintf( output_fp,	"\n"
	" %4d %5d %11.4LE %11.4LE %11.4LE %11.4LE"
	" %11.4LE %11.4LE %11.4LE %11.4LE %11.4LE",
	irow2, isc1+1, creall(vlt), cimagl(vlt), creall(cux), cimagl(cux),
	creall(zped), cimagl(zped), creall(ymit), cimagl(ymit), pwr );

  } /* for( i = 0; i < nvqd; i++ ) */

  /* Free network buffers */
  free_ptr( (void *)&ipnt );
  free_ptr( (void *)&nteqa );
  free_ptr( (void *)&ntsca );
  free_ptr( (void *)&vsrc );
  free_ptr( (void *)&rhs );
  free_ptr( (void *)&cmn );
  free_ptr( (void *)&rhnt );
  free_ptr( (void *)&rhnx );

  return;
}

/*-----------------------------------------------------------------------*/

