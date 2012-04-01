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

/*** common data are implemented as global variables ***/
/* common  /data/ */
int n, np, m, mp, ipsym, npm, np2m, np3m; /* n+m,n+2m,n+3m */
int *icon1, *icon2, *itag;
long double *x, *y, *z, *si, *bi;
long double *x2, *y2, *z2, *cab, *sab, *salp;
long double *t1x, *t1y, *t1z, *t2x, *t2y, *t2z;
long double *px, *py, *pz, *pbi, *psalp;
long double wlam;

/* common  /cmb/ */
complex long double *cm;

/* common  /matpar/ */
int icase, npblk, nlast;
int imat, nbbx, npbx, nlbx, nbbl, npbl, nlbl;

/* common  /save/ */
int *ip;
long double epsr, sig, scrwlt, scrwrt, fmhz;

/* common  /crnt/ */
long double *air, *aii, *bir, *bii, *cir, *cii;
complex long double *cur;

/* common  /gnd/ */
int ksymp, ifar, iperf, nradl;
long double t2, cl, ch, scrwl, scrwr;
complex long double zrati, zrati2, t1, frati;

/* common  /zload/ */
int nload;
complex long double *zarray;

/* common  /yparm/ */
int ncoup, icoup, *nctag, *ncseg;
complex long double *y11a, *y12a;

/* common  /segj/ */
int *jco, jsno, nscon, maxcon; /* Max. no. connections */
long double *ax, *bx, *cx;

/* common  /vsorc/ */
int *ivqd, *isant, *iqds, nvqd, nsant, nqds;
complex long double *vqd, *vqds, *vsant;

/* common  /netcx/ */
int masym, neq, npeq, neq2, nonet, ntsol, nprint;
int *iseg1, *iseg2, *ntyp;
long double *x11r, *x11i, *x12r;
long double *x12i, *x22r, *x22i;
long double pin, pnls;
complex long double zped;

/* common  /fpat/ */
int near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp, inor, iax, ixtyp;
long double thets, phis, dth, dph, rfld, gnor, clt, cht, epsr2, sig2;
long double xpr6, pinr, pnlr, ploss, xnr, ynr, znr, dxnr, dynr, dznr;

/*common  /ggrid/ */
extern int nxa[3], nya[3];
extern long double dxa[3], dya[3], xsa[3], ysa[3];
extern complex long double epscf, *ar1, *ar2, *ar3;

/* common  /gwav/ */
long double r1, r2, zmh, zph;
complex long double u, u2, xx1, xx2;

/* common  /plot/ */
int iplp1, iplp2, iplp3, iplp4;

/* common  /dataj/ */
int iexk, ind1, indd1, ind2, indd2, ipgnd;
long double s, b, xj, yj, zj, cabj, sabj, salpj, rkh;
long double t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
complex long double  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;

/* common  /smat/ */
int nop; /* My addition */
complex long double *ssx;

/* common  /incom/ */
int isnor;
long double xo, yo, zo, sn, xsn, ysn;

/* common  /tmi/ */
int ija; /* changed to ija to avoid conflict */
long double zpk, rkb2;

/*common  /tmh/ */
long double zpka, rhks;

/* pointers to input/output files */
FILE *input_fp=NULL, *output_fp=NULL, *plot_fp=NULL;

/* signal handler */
static void sig_handler( int signal );

/*-------------------------------------------------------------------*/

int main( int argc, char **argv )
{
  char infile[81] = "", otfile[81] = "";
  char ain[3], line_buf[81];

  /* input card mnemonic list */
#define NUM_CMNDS  20
  char *atst[NUM_CMNDS] =
  {
    "FR", "LD", "GN", "EX", "NT", "TL", \
    "XQ", "GD", "RP", "NX", "PT", "KH", \
    "NE", "NH", "PQ", "EK", "CP", "PL", \
    "EN", "WG"
  };

  char *hpol[3] = { "LINEAR", "RIGHT", "LEFT" };
  char *pnet[3] = { "        ", "STRAIGHT", " CROSSED" };

  int *ldtyp, *ldtag, *ldtagf, *ldtagt;
  int ifrtmw, ifrtmp, mpcnt, ib11=0, ic11=0, id11=0, ix11, igo, nfrq;
  int iexk, jump, iptflg, iptflq, iped, iflow, itmp1, iresrv;
  int itmp3, itmp2, itmp4, nthi=0, nphi=0, iptag=0, iptagf=0, iptagt=0;
  int iptaq=0, iptaqf=0, iptaqt=0, nphic=0, inc=0;
  int i, j, itmp5, nthic=0, mhz=0, ifrq=0, isave=0;

  int
    igox,       /* used in place of "igo" in freq loop */
    next_job,   /* start next job (next sructure) flag */
    idx,        /* general purpose index    */
    ain_num,    /* ain mnemonic as a number */
    jmp_iloop,  /* jump to input loop flag  */
    jmp_floop=0,/* jump to freq. loop flag  */
    mreq;       /* Size req. for malloc's   */

  long double *zlr, *zli, *zlc, *fnorm;
  long double *xtemp, *ytemp, *ztemp, *sitemp, *bitemp;
  long double fmhz1, rkh, tmp1, delfrq=0., tmp2, tmp3, tmp4, tmp5, tmp6;
  long double xpr1=0., xpr2=0., xpr3=0., xpr4=0., xpr5=0.;
  long double zpnorm=0., thetis=0., phiss=0., extim;
  long double tim1, tim, tim2, etha, fr, fr2, cmag, ph, ethm, ephm, epha;
  complex long double eth, eph, curi, ex, ey, ez, epsc;

  /* getopt() variables */
  extern char *optarg;
  extern int optind, opterr, optopt;
  int option;

  /*** signal handler related code ***/
  /* new and old actions for sigaction() */
  struct sigaction sa_new, sa_old;


  /* initialize new actions */
  sa_new.sa_handler = sig_handler;
  sigemptyset( &sa_new.sa_mask );
  sa_new.sa_flags = 0;

  /* register function to handle signals */
  sigaction( SIGINT,  &sa_new, &sa_old );
  sigaction( SIGSEGV, &sa_new, 0 );
  sigaction( SIGFPE,  &sa_new, 0 );
  sigaction( SIGTERM, &sa_new, 0 );
  sigaction( SIGABRT, &sa_new, 0 );

  /*** command line arguments handler ***/
  if( argc == 1 )
  {
    usage();
    exit(-1);
  }

  /* process command line options */
  while( (option = getopt(argc, argv, "i:o:hv") ) != -1 )
  {
    switch( option )
    {
      case 'i' : /* specify input file name */
	if( strlen(optarg) > 75 )
	  abort_on_error(-1);
	strcpy( infile, optarg );
	break;

      case 'o' : /* specify output file name */
	if( strlen(optarg) > 75 )
	  abort_on_error(-2);
	strcpy( otfile, optarg );
	break;

      case 'h' : /* print usage and exit */
	usage();
	exit(0);

      case 'v' : /* print nec2c version */
	puts( version );
	exit(0);

      default: /* print usage and exit */
	usage();
	exit(-1);

    } /* end of switch( option ) */

  } /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */

  /*** open input file ***/
  if( (input_fp = fopen(infile, "r")) == NULL )
  {
    char mesg[88] = "nec2c: ";

    strcat( mesg, infile );
    perror( mesg );
    exit(-1);
  }

  /* make an output file name if not */
  /* specified by user on invocation */
  if( strlen( otfile ) == 0 )
  {
    /* strip file name extension if there is one */
    idx = 0;
    while( (infile[++idx] != '.') && (infile[idx] != '\0') );
    infile[idx] = '\0';

    /* make the output file name from input file */
    strcpy( otfile, infile );
    strcat( otfile, ".out" ); /* add extension */
  }

  /* open output file */
  if( (output_fp = fopen(otfile, "w")) == NULL )
  {
    char mesg[88] = "nec2c: ";

    strcat( mesg, otfile );
    perror( mesg );
    exit(-1);
  }

  /*** here we had code to read interactively input/output ***/
  /*** file names. this is done non-interactively above.   ***/

  secnds( &extim );

  /* Null buffer pointers */
  /* type int */
  icon1 = icon2 = ncseg = nctag = ivqd = isant = iqds = NULL;
  itag = ip = ldtyp = ldtag = ldtagf = ldtagt = jco = NULL;
  /* type long double */
  air = aii = bir = bii = cir = cii = zlr = zli = zlc = fnorm = NULL;
  ax = bx = cx = xtemp = ytemp = ztemp = sitemp = bitemp = NULL;
  x = y = z = si = bi = x2 = y2 = z2 = cab = sab = salp = NULL;
  t1x = t1y = t1z = t2x = t2y = t2z = px = py = pz = pbi = psalp = NULL;
  /* type complex long double */
  ar1 = ar2 = ar3 = cur = cm = zarray = NULL;
  y11a = y12a = vqd = vqds = vsant = ssx = NULL;

  /* Allocate some buffers */
  mem_alloc( (void *)&ar1, sizeof(complex long double)*11*10*4 );
  mem_alloc( (void *)&ar2, sizeof(complex long double)*17*5*4 );
  mem_alloc( (void *)&ar3, sizeof(complex long double)*9*8*4 );


  /* l_1: */
  /* main execution loop, exits at various points */
  /* depending on error conditions or end of jobs */
  while( TRUE )
  {
    ifrtmw=0;
    ifrtmp=0;

    /* print the nec2c header to output file */
    fprintf( output_fp,	"\n\n\n"
	"                              "
	" __________________________________________\n"
	"                              "
	"|                                          |\n"
	"                              "
	"|  NUMERICAL ELECTROMAGNETICS CODE (nec2c) |\n"
	"                              "
	"|   Translated to 'C' in Double Precision  |\n"
	"                              "
	"|__________________________________________|\n" );

    /* read a line from input file */
    if( load_line(line_buf, input_fp) == EOF )
      abort_on_error(-3);

    /* separate card's id mnemonic */
    strncpy( ain, line_buf, 2 );
    ain[2] = '\0';

    /* if its a "cm" or "ce" card start reading comments */
    if( (strcmp(ain, "CM") == 0) ||
	(strcmp(ain, "CE") == 0) )
    {
      fprintf( output_fp, "\n\n\n"
	  "                               "
	  "---------------- COMMENTS ----------------\n" );

      /* write comment to output file */
      fprintf( output_fp,
	  "                              %s\n",
	  &line_buf[2] );

      /* Keep reading till a non "CM" card */
      while( strcmp(ain, "CM") == 0 )
      {
	/* read a line from input file */
	if( load_line(line_buf, input_fp) == EOF )
	  abort_on_error(-3);

	/* separate card's id mnemonic */
	strncpy( ain, line_buf, 2 );
	ain[2] = '\0';

	/* write comment to output file */
	fprintf( output_fp,
	    "                              %s\n",
	    &line_buf[2] );

      } /* while( strcmp(ain, "CM") == 0 ) */

      /* no "ce" card at end of comments */
      if( strcmp(ain, "CE") != 0 )
      {
	fprintf( output_fp,
	    "\n\n  ERROR: INCORRECT LABEL FOR A COMMENT CARD" );
	abort_on_error(-4);
      }

    } /* if( strcmp(ain, "CM") == 0 ... */
    else
      rewind( input_fp );

    /* Free some buffer pointers.
     * These are allocated by realloc()
     * so they need to be free()'d
     * before reallocation for a new job
     */
    free_ptr( (void *)&itag );
    free_ptr( (void *)&fnorm );
    free_ptr( (void *)&ldtyp );
    free_ptr( (void *)&ldtag );
    free_ptr( (void *)&ldtagf );
    free_ptr( (void *)&ldtagt );
    free_ptr( (void *)&zlr );
    free_ptr( (void *)&zli );
    free_ptr( (void *)&zlc );
    free_ptr( (void *)&jco );
    free_ptr( (void *)&ax );
    free_ptr( (void *)&bx );
    free_ptr( (void *)&cx );
    free_ptr( (void *)&ivqd );
    free_ptr( (void *)&iqds );
    free_ptr( (void *)&vqd );
    free_ptr( (void *)&vqds );
    free_ptr( (void *)&isant );
    free_ptr( (void *)&vsant );
    free_ptr( (void *)&x );
    free_ptr( (void *)&y );
    free_ptr( (void *)&z );
    free_ptr( (void *)&x2 );
    free_ptr( (void *)&y2 );
    free_ptr( (void *)&z2 );
    free_ptr( (void *)&px );
    free_ptr( (void *)&py );
    free_ptr( (void *)&pz );
    free_ptr( (void *)&t1x );
    free_ptr( (void *)&t1y );
    free_ptr( (void *)&t1z );
    free_ptr( (void *)&t2x );
    free_ptr( (void *)&t2y );
    free_ptr( (void *)&t2z );
    free_ptr( (void *)&si );
    free_ptr( (void *)&bi );
    free_ptr( (void *)&cab );
    free_ptr( (void *)&sab );
    free_ptr( (void *)&salp );
    free_ptr( (void *)&pbi );
    free_ptr( (void *)&psalp );

    /* initializations etc from original fortran code */
    mpcnt=0;
    imat=0;

    /* set up geometry data in subroutine datagn */
    datagn();
    iflow=1;

    /* Allocate some buffers */
    mreq = npm * sizeof(long double);
    mem_alloc( (void *)&air, mreq );
    mem_alloc( (void *)&aii, mreq );
    mem_alloc( (void *)&bir, mreq );
    mem_alloc( (void *)&bii, mreq );
    mem_alloc( (void *)&cir, mreq );
    mem_alloc( (void *)&cii, mreq );
    mem_alloc( (void *)&xtemp,  mreq );
    mem_alloc( (void *)&ytemp,  mreq );
    mem_alloc( (void *)&ztemp,  mreq );
    mem_alloc( (void *)&sitemp, mreq );
    mem_alloc( (void *)&bitemp, mreq );

    mreq = np2m * sizeof(int);
    mem_alloc( (void *)&ip, mreq );

    mreq = np3m * sizeof( complex long double);
    mem_alloc( (void *)&cur, mreq );

    /* Matrix parameters */
    if( imat == 0)
    {
      neq= n+2*m;
      neq2=0;
      ib11=0;
      ic11=0;
      id11=0;
      ix11=0;
    }

    fprintf( output_fp, "\n\n\n" );

    /* default values for input parameters and flags */
    npeq= np+2*mp;
    iplp1=0;
    iplp2=0;
    iplp3=0;
    iplp4=0;
    igo=1;
    nfrq=1;
    rkh=1.;
    iexk=0;
    ixtyp=0;
    nload=0;
    nonet=0;
    near=-1;
    iptflg=-2;
    iptflq=-1;
    ifar=-1;
    zrati=CPLX_10;
    iped=0;
    ncoup=0;
    icoup=0;
    fmhz= CVEL;
    ksymp=1;
    nradl=0;
    iperf=0;

    /* l_14: */

    /* main input section, exits at various points */
    /* depending on error conditions or end of job */
    next_job = FALSE;
    while( ! next_job )
    {
      jmp_iloop = FALSE;

      /* main input section - standard read statement - jumps */
      /* to appropriate section for specific parameter set up */
      readmn( ain, &itmp1, &itmp2, &itmp3, &itmp4,
	  &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6 );

      /* If its an "XT" card, exit */
      if( strcmp(ain, "XT" ) == 0 )
      {
	fprintf( stderr,
	    "\nnec2c: Exiting after an \"XT\" command in main()\n" );
	fprintf( output_fp,
	    "\n\n  nec2c: Exiting after an \"XT\" command in main()" );
	stop(0);
      }

      mpcnt++;
      fprintf( output_fp,
	  "\n  DATA CARD No: %3d "
	  "%s %3d %5d %5d %5d %12.5LE %12.5LE %12.5LE %12.5LE %12.5LE %12.5LE",
	  mpcnt, ain, itmp1, itmp2, itmp3, itmp4,
	  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );

      /* identify card id mnemonic (except "ce" and "cm") */
      for( ain_num = 0; ain_num < NUM_CMNDS; ain_num++ )
	if( strncmp( ain, atst[ain_num], 2) == 0 )
	  break;

      /* take action according to card id mnemonic */
      switch( ain_num )
      {
	case 0: /* "fr" card, frequency parameters */

	  ifrq= itmp1;
	  nfrq= itmp2;
	  if( nfrq == 0)
	    nfrq=1;
	  fmhz= tmp1;
	  delfrq= tmp2;
	  if( iped == 1)
	    zpnorm=0.;
	  igo=1;
	  iflow=1;

	  continue; /* continue card input loop */

	case 1: /* "ld" card, loading parameters */
	  {
	    int idx;

	    if( iflow != 3 )
	    {
	      iflow=3;
	      /* Free loading buffers */
	      nload=0;
	      free_ptr( (void *)&ldtyp );
	      free_ptr( (void *)&ldtag );
	      free_ptr( (void *)&ldtagf );
	      free_ptr( (void *)&ldtagt );
	      free_ptr( (void *)&zlr );
	      free_ptr( (void *)&zli );
	      free_ptr( (void *)&zlc );

	      if( igo > 2 )
		igo=2;
	      if( itmp1 == -1 )
		continue; /* continue card input loop */
	    }

	    /* Reallocate loading buffers */
	    nload++;
	    idx = nload * sizeof(int);
	    mem_realloc( (void *)&ldtyp,  idx );
	    mem_realloc( (void *)&ldtag,  idx );
	    mem_realloc( (void *)&ldtagf, idx );
	    mem_realloc( (void *)&ldtagt, idx );
	    idx = nload * sizeof(long double);
	    mem_realloc( (void *)&zlr, idx );
	    mem_realloc( (void *)&zli, idx );
	    mem_realloc( (void *)&zlc, idx );

	    idx = nload-1;
	    ldtyp[idx]= itmp1;
	    ldtag[idx]= itmp2;
	    if( itmp4 == 0)
	      itmp4= itmp3;
	    ldtagf[idx]= itmp3;
	    ldtagt[idx]= itmp4;

	    if( itmp4 < itmp3 )
	    {
	      fprintf( output_fp,
		  "\n\n  DATA FAULT ON LOADING CARD No: %d: ITAG "
		  "STEP1: %d IS GREATER THAN ITAG STEP2: %d",
		  nload, itmp3, itmp4 );
	      stop(-1);
	    }

	    zlr[idx]= tmp1;
	    zli[idx]= tmp2;
	    zlc[idx]= tmp3;
	  }

	  continue; /* continue card input loop */

	case 2: /* "gn" card, ground parameters under the antenna */

	  iflow=4;

	  if( igo > 2)
	    igo=2;

	  if( itmp1 == -1 )
	  {
	    ksymp=1;
	    nradl=0;
	    iperf=0;
	    continue; /* continue card input loop */
	  }

	  iperf= itmp1;
	  nradl= itmp2;
	  ksymp=2;
	  epsr= tmp1;
	  sig= tmp2;

	  if( nradl != 0)
	  {
	    if( iperf == 2)
	    {
	      fprintf( output_fp,
		  "\n\n  RADIAL WIRE G.S. APPROXIMATION MAY "
		  "NOT BE USED WITH SOMMERFELD GROUND OPTION" );
	      stop(-1);
	    }

	    scrwlt= tmp3;
	    scrwrt= tmp4;
	    continue; /* continue card input loop */
	  }

	  epsr2= tmp3;
	  sig2= tmp4;
	  clt= tmp5;
	  cht= tmp6;

	  continue; /* continue card input loop */

	case 3: /* "ex" card, excitation parameters */

	  if( iflow != 5)
	  {
	    /* Free vsource buffers */
	    free_ptr( (void *)&ivqd );
	    free_ptr( (void *)&iqds );
	    free_ptr( (void *)&vqd );
	    free_ptr( (void *)&vqds );
	    free_ptr( (void *)&isant );
	    free_ptr( (void *)&vsant );

	    nsant=0;
	    nvqd=0;
	    iped=0;
	    iflow=5;
	    if( igo > 3)
	      igo=3;
	  }

	  masym= itmp4/10;
	  if( (itmp1 == 0) || (itmp1 == 5) )
	  {
	    ixtyp= itmp1;
	    ntsol=0;

	    if( ixtyp != 0)
	    {
	      nvqd++;
	      mem_realloc( (void *)&ivqd, nvqd * sizeof(int) );
	      mem_realloc( (void *)&iqds, nvqd * sizeof(int) );
	      mem_realloc( (void *)&vqd,  nvqd * sizeof(complex long double) );
	      mem_realloc( (void *)&vqds, nvqd * sizeof(complex long double) );

	      {
		int indx = nvqd-1;

		ivqd[indx]= isegno( itmp2, itmp3);
		vqd[indx]= cmplx( tmp1, tmp2);
		if( cabsl( vqd[indx]) < 1.e-20)
		  vqd[indx] = CPLX_10;

		iped= itmp4- masym*10;
		zpnorm= tmp3;
		if( (iped == 1) && (zpnorm > 0.0) )
		  iped=2;
		continue; /* continue card input loop */
	      }

	    } /* if( ixtyp != 0) */

	    nsant++;
	    mem_realloc( (void *)&isant, nsant * sizeof(int) );
	    mem_realloc( (void *)&vsant, nsant * sizeof(complex long double) );

	    {
	      int indx = nsant-1;

	      isant[indx]= isegno( itmp2, itmp3);
	      vsant[indx]= cmplx( tmp1, tmp2);
	      if( cabsl( vsant[indx]) < 1.e-20)
		vsant[indx] = CPLX_10;

	      iped= itmp4- masym*10;
	      zpnorm= tmp3;
	      if( (iped == 1) && (zpnorm > 0.0) )
		iped=2;
	      continue; /* continue card input loop */
	    }

	  } /* if( (itmp1 <= 0) || (itmp1 == 5) ) */

	  if( (ixtyp == 0) || (ixtyp == 5) )
	    ntsol=0;

	  ixtyp= itmp1;
	  nthi= itmp2;
	  nphi= itmp3;
	  xpr1= tmp1;
	  xpr2= tmp2;
	  xpr3= tmp3;
	  xpr4= tmp4;
	  xpr5= tmp5;
	  xpr6= tmp6;
	  nsant=0;
	  nvqd=0;
	  thetis= xpr1;
	  phiss= xpr2;

	  continue; /* continue card input loop */

	case 4: case 5: /* "nt" & "tl" cards, network parameters */
	  {
	    int idx;

	    if( iflow != 6)
	    {
	      nonet=0;
	      ntsol=0;
	      iflow=6;

	      /* Free network buffers */
	      free_ptr( (void *)&ntyp );
	      free_ptr( (void *)&iseg1 );
	      free_ptr( (void *)&iseg2 );
	      free_ptr( (void *)&x11r );
	      free_ptr( (void *)&x11i );
	      free_ptr( (void *)&x12r );
	      free_ptr( (void *)&x12i );
	      free_ptr( (void *)&x22r );
	      free_ptr( (void *)&x22i );

	      if( igo > 3)
		igo=3;

	      if( itmp2 == -1 )
		continue; /* continue card input loop */
	    }

	    /* Re-allocate network buffers */
	    nonet++;
	    idx = nonet * sizeof(int);
	    mem_realloc( (void *)&ntyp, idx );
	    mem_realloc( (void *)&iseg1, idx );
	    mem_realloc( (void *)&iseg2, idx );
	    idx = nonet * sizeof(long double);
	    mem_realloc( (void *)&x11r, idx );
	    mem_realloc( (void *)&x11i, idx );
	    mem_realloc( (void *)&x12r, idx );
	    mem_realloc( (void *)&x12i, idx );
	    mem_realloc( (void *)&x22r, idx );
	    mem_realloc( (void *)&x22i, idx );

	    idx = nonet-1;
	    if( ain_num == 4 )
	      ntyp[idx]=1;
	    else
	      ntyp[idx]=2;

	    iseg1[idx]= isegno( itmp1, itmp2);
	    iseg2[idx]= isegno( itmp3, itmp4);
	    x11r[idx]= tmp1;
	    x11i[idx]= tmp2;
	    x12r[idx]= tmp3;
	    x12i[idx]= tmp4;
	    x22r[idx]= tmp5;
	    x22i[idx]= tmp6;

	    if( (ntyp[idx] == 1) || (tmp1 > 0.) )
	      continue; /* continue card input loop */

	    ntyp[idx]=3;
	    x11r[idx]=- tmp1;

	    continue; /* continue card input loop */
	  }

	case 6: /* "xq" execute card - calc. including radiated fields */

	  if( ((iflow == 10) && (itmp1 == 0)) ||
	      ((nfrq  ==  1) && (itmp1 == 0) && (iflow > 7)) )
	    continue; /* continue card input loop */

	  if( itmp1 == 0)
	  {
	    if( iflow > 7)
	      iflow=11;
	    else
	      iflow=7;
	  }
	  else
	  {
	    ifar=0;
	    rfld=0.;
	    ipd=0;
	    iavp=0;
	    inor=0;
	    iax=0;
	    nth=91;
	    nph=1;
	    thets=0.;
	    phis=0.;
	    dth=1.0;
	    dph=0.;

	    if( itmp1 == 2)
	      phis=90.;

	    if( itmp1 == 3)
	    {
	      nph=2;
	      dph=90.;
	    }

	  } /* if( itmp1 == 0) */

	  break;

	case 7: /* "gd" card, ground representation */

	  epsr2= tmp1;
	  sig2= tmp2;
	  clt= tmp3;
	  cht= tmp4;
	  iflow=9;

	  continue; /* continue card input loop */

	case 8: /* "rp" card, standard observation angle parameters */

	  ifar= itmp1;
	  nth= itmp2;
	  nph= itmp3;

	  if( nth == 0)
	    nth=1;
	  if( nph == 0)
	    nph=1;

	  ipd= itmp4/10;
	  iavp= itmp4- ipd*10;
	  inor= ipd/10;
	  ipd= ipd- inor*10;
	  iax= inor/10;
	  inor= inor- iax*10;

	  if( iax != 0)
	    iax=1;
	  if( ipd != 0)
	    ipd=1;
	  if( (nth < 2) || (nph < 2) || (ifar == 1) )
	    iavp=0;

	  thets= tmp1;
	  phis= tmp2;
	  dth= tmp3;
	  dph= tmp4;
	  rfld= tmp5;
	  gnor= tmp6;
	  iflow=10;

	  break;

	case 9: /* "nx" card, do next job */
	  next_job = TRUE;
	  continue; /* continue card input loop */

	case 10: /* "pt" card, print control for current */

	  iptflg= itmp1;
	  iptag= itmp2;
	  iptagf= itmp3;
	  iptagt= itmp4;

	  if( (itmp3 == 0) && (iptflg != -1) )
	    iptflg=-2;
	  if( itmp4 == 0)
	    iptagt= iptagf;

	  continue; /* continue card input loop */

	case 11: /* "kh" card, matrix integration limit */

	  rkh= tmp1;
	  if( igo > 2)
	    igo=2;
	  iflow=1;

	  continue; /* continue card input loop */

	case 12: case 13:  /* "ne"/"nh" cards, near field calculation parameters */

	  if( ain_num == 13 )
	    nfeh=1;
	  else
	    nfeh=0;

	  if( (iflow == 8) && (nfrq != 1) )
	  {
	    fprintf( output_fp,
		"\n\n  WHEN MULTIPLE FREQUENCIES ARE REQUESTED, "
		"ONLY ONE NEAR FIELD CARD CAN BE USED -"
		"\n  LAST CARD READ WILL BE USED" );
	  }

	  near= itmp1;
	  nrx= itmp2;
	  nry= itmp3;
	  nrz= itmp4;
	  xnr= tmp1;
	  ynr= tmp2;
	  znr= tmp3;
	  dxnr= tmp4;
	  dynr= tmp5;
	  dznr= tmp6;
	  iflow=8;

	  if( nfrq != 1)
	    continue; /* continue card input loop */

	  break;

	case 14: /* "pq" card, write control for charge */

	  iptflq= itmp1;
	  iptaq= itmp2;
	  iptaqf= itmp3;
	  iptaqt= itmp4;

	  if( (itmp3 == 0) && (iptflq != -1) )
	    iptflq=-2;
	  if( itmp4 == 0)
	    iptaqt= iptaqf;

	  continue; /* continue card input loop */

	case 15: /* "ek" card,  extended thin wire kernel option */

	  iexk=1;
	  if( itmp1 == -1)
	    iexk=0;
	  if( igo > 2)
	    igo=2;
	  iflow=1;

	  continue; /* continue card input loop */

	case 16: /* "cp" card, maximum coupling between antennas */

	  if( iflow != 2)
	  {
	    ncoup=0;
	    free_ptr( (void *)&nctag );
	    free_ptr( (void *)&ncseg );
	    free_ptr( (void *)&y11a );
	    free_ptr( (void *)&y12a );
	  }

	  icoup=0;
	  iflow=2;

	  if( itmp2 == 0)
	    continue; /* continue card input loop */

	  ncoup++;
	  mem_realloc( (void *)&nctag, (ncoup) * sizeof(int) );
	  mem_realloc( (void *)&ncseg, (ncoup) * sizeof(int) );
	  nctag[ncoup-1]= itmp1;
	  ncseg[ncoup-1]= itmp2;

	  if( itmp4 == 0)
	    continue; /* continue card input loop */

	  ncoup++;
	  mem_realloc( (void *)&nctag, (ncoup) * sizeof(int) );
	  mem_realloc( (void *)&ncseg, (ncoup) * sizeof(int) );
	  nctag[ncoup-1]= itmp3;
	  ncseg[ncoup-1]= itmp4;

	  continue; /* continue card input loop */

	case 17: /* "pl" card, plot flags */

	  iplp1= itmp1;
	  iplp2= itmp2;
	  iplp3= itmp3;
	  iplp4= itmp4;

	  if( plot_fp == NULL )
	  {
	    char plotfile[81];

	    /* Make a plot file name */
	    strcpy( plotfile, infile );
	    strcat( plotfile, ".plt" );

	    /* Open plot file */
	    if( (plot_fp = fopen(plotfile, "w")) == NULL )
	    {
	      char mesg[88] = "nec2c: ";

	      strcat( mesg, plotfile );
	      perror( mesg );
	      exit(-1);
	    }
	  }

	  continue; /* continue card input loop */

	case 19: /* "wg" card, not supported */
	  abort_on_error(-5);

	default:
	  if( ain_num != 18 )
	  {
	    fprintf( output_fp,
		"\n\n  FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION" );
	    stop(-1);
	  }

	  /******************************************************
	   *** normal exit of nec2c when all jobs complete ok ***
	   ******************************************************/

	  /* time the process */
	  secnds( &tmp1 );
	  tmp1 -= extim;
	  fprintf( output_fp, "\n\n  TOTAL RUN TIME: %d msec", (int)tmp1 );
	  stop(0);

      } /* switch( ain_num ) */

      /**************************************
       *** end of the main input section. ***
       *** beginning of frequency do loop ***
       **************************************/

      /* Allocate to normalization buffer */
      {
	int mreq1, mreq2;

	mreq1 = mreq2 = 0;
	if( iped )
	  mreq1 = 4*nfrq * sizeof(long double);
	if( iptflg >= 2 )
	  mreq2 = nthi*nphi * sizeof(long double);

	if( (mreq1 > 0) || (mreq2 > 0) )
	{
	  if( mreq1 > mreq2 )
	    mem_alloc( (void *)&fnorm, mreq1 );
	  else
	    mem_alloc( (void *)&fnorm, mreq2 );
	}
      }

      /* igox is used in place of "igo" in the   */
      /* freq loop. below is a special igox case */
      if( ((ain_num == 6) || (ain_num == 8)) && (igo == 5) )
	igox = 6;
      else
	igox = igo;

      switch( igox )
      {
	case 1: /* label 41 */
	  /* Memory allocation for primary interacton matrix. */
	  iresrv = np2m * (np+2*mp);
	  mem_alloc( (void *)&cm, iresrv * sizeof(complex long double) );

	  /* Memory allocation for symmetry array */
	  nop = neq/npeq;
	  mem_alloc( (void *)&ssx, nop*nop * sizeof( complex long double) );

	  mhz=1;

	  if( (n != 0) && (ifrtmw != 1) )
	  {
	    ifrtmw=1;
	    for( i = 0; i < n; i++ )
	    {
	      xtemp[i]= x[i];
	      ytemp[i]= y[i];
	      ztemp[i]= z[i];
	      sitemp[i]= si[i];
	      bitemp[i]= bi[i];
	    }
	  }

	  if( (m != 0) && (ifrtmp != 1) )
	  {
	    ifrtmp=1;
	    for( i = 0; i < m; i++ )
	    {
	      j = i+n;
	      xtemp[j]= px[i];
	      ytemp[j]= py[i];
	      ztemp[j]= pz[i];
	      bitemp[j]= pbi[i];
	    }
	  }

	  fmhz1= fmhz;

	  /* irngf is not used (NGF function not implemented) */
	  if( imat == 0)
	    fblock( npeq, neq, iresrv, ipsym);

	  /* label 42 */
	  /* frequency do loop */
	  do
	  {
	    jmp_floop = FALSE;

	    if( mhz != 1)
	    {
	      if( ifrq == 1)
		fmhz *= delfrq;
	      else
		fmhz += delfrq;
	    }

	    fr= fmhz/ CVEL;
	    wlam= CVEL/ fmhz;
	    fprintf( output_fp, "\n\n\n"
		"                               "
		"--------- FREQUENCY --------\n"
		"                                "
		"FREQUENCY :%11.4LE MHz\n"
		"                                "
		"WAVELENGTH:%11.4LE Mtr", fmhz, wlam );

	    fprintf( output_fp, "\n\n"
		"                        "
		"APPROXIMATE INTEGRATION EMPLOYED FOR SEGMENTS \n"
		"                        "
		"THAT ARE MORE THAN %.3LF WAVELENGTHS APART", rkh );

	    if( iexk == 1)
	      fprintf( output_fp, "\n"
		  "                        "
		  "THE EXTENDED THIN WIRE KERNEL WILL BE USED" );

	    /* frequency scaling of geometric parameters */
	    if( n != 0)
	    {
	      for( i = 0; i < n; i++ )
	      {
		x[i]= xtemp[i]* fr;
		y[i]= ytemp[i]* fr;
		z[i]= ztemp[i]* fr;
		si[i]= sitemp[i]* fr;
		bi[i]= bitemp[i]* fr;
	      }
	    }

	    if( m != 0)
	    {
	      fr2= fr* fr;
	      for( i = 0; i < m; i++ )
	      {
		j = i+n;
		px[i]= xtemp[j]* fr;
		py[i]= ytemp[j]* fr;
		pz[i]= ztemp[j]* fr;
		pbi[i]= bitemp[j]* fr2;
	      }
	    }

	    igo = 2;

	    /* label 46 */
	    case 2: /* structure segment loading */
	    fprintf( output_fp, "\n\n\n"
		"                          "
		"------ STRUCTURE IMPEDANCE LOADING ------" );

	    if( nload != 0)
	      load( ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc );

	    if( nload == 0 )
	      fprintf( output_fp, "\n"
		  "                                 "
		  "THIS STRUCTURE IS NOT LOADED" );

	    fprintf( output_fp, "\n\n\n"
		"                            "
		"-------- ANTENNA ENVIRONMENT --------" );

	    if( ksymp != 1)
	    {
	      frati=CPLX_10;

	      if( iperf != 1)
	      {
		if( sig < 0.)
		  sig=- sig/(59.96*wlam);

		epsc= cmplx( epsr, -sig*wlam*59.96);
		zrati=1./ csqrtl( epsc);
		u= zrati;
		u2= u* u;

		if( nradl != 0)
		{
		  scrwl= scrwlt/ wlam;
		  scrwr= scrwrt/ wlam;
		  t1= CPLX_01*2367.067/ (long double)nradl;
		  t2= scrwr* (long double)nradl;

		  fprintf( output_fp, "\n"
		      "                            "
		      "RADIAL WIRE GROUND SCREEN\n"
		      "                            "
		      "%d WIRES\n"
		      "                            "
		      "WIRE LENGTH: %8.2LF METERS\n"
		      "                            "
		      "WIRE RADIUS: %10.3LE METERS",
		      nradl, scrwlt, scrwrt );

		  fprintf( output_fp, "\n"
		      "                            "
		      "MEDIUM UNDER SCREEN -" );

		} /* if( nradl != 0) */

		if( iperf != 2)
		  fprintf( output_fp, "\n"
		      "                            "
		      "FINITE GROUND - REFLECTION COEFFICIENT APPROXIMATION" );
		else
		{
		  somnec( epsr, sig, fmhz );
		  frati=( epsc-1.)/( epsc+1.);
		  if( cabsl(( epscf- epsc)/ epsc) >= 1.0e-3 )
		  {
		    fprintf( output_fp,
			"\n ERROR IN GROUND PARAMETERS -"
			"\n COMPLEX DIELECTRIC CONSTANT FROM FILE IS: %12.5LE%+12.5LEj"
			"\n                                REQUESTED: %12.5LE%+12.5LEj",
			creall(epscf), cimagl(epscf), creall(epsc), cimagl(epsc) );
		    stop(-1);
		  }

		  fprintf( output_fp, "\n"
		      "                            "
		      "FINITE GROUND - SOMMERFELD SOLUTION" );

		} /* if( iperf != 2) */

		fprintf( output_fp, "\n"
		    "                            "
		    "RELATIVE DIELECTRIC CONST: %.3LF\n"
		    "                            "
		    "CONDUCTIVITY: %10.3LE MHOS/METER\n"
		    "                            "
		    "COMPLEX DIELECTRIC CONSTANT: %11.4LE%+11.4LEj",
		    epsr, sig, creall(epsc), cimagl(epsc) );

	      } /* if( iperf != 1) */
	      else
		fprintf( output_fp, "\n"
		    "                            "
		    "PERFECT GROUND" );

	    } /* if( ksymp != 1) */
	    else
	      fprintf( output_fp, "\n"
		  "                            "
		  "FREE SPACE" );

	    /* label 50 */
	    /* fill and factor primary interaction matrix */
	    secnds( &tim1 );
	    cmset( neq, cm, rkh, iexk );
	    secnds( &tim2 );
	    tim= tim2- tim1;
	    factrs( npeq, neq, cm, ip );
	    secnds( &tim1 );
	    tim2= tim1- tim2;
	    fprintf( output_fp, "\n\n\n"
		"                             "
		"---------- MATRIX TIMING ----------\n"
		"                               "
		"FILL: %d msec  FACTOR: %d msec",
		(int)tim, (int)tim2 );

	    igo=3;
	    ntsol=0;

	    /* label 53 */
	    case 3: /* excitation set up (right hand side, -e inc.) */

	    nthic=1;
	    nphic=1;
	    inc=1;
	    nprint=0;

	    /* l_54 */
	    do
	    {
	      if( (ixtyp != 0) && (ixtyp != 5) )
	      {
		if( (iptflg <= 0) || (ixtyp == 4) )
		  fprintf( output_fp, "\n\n\n"
		      "                             "
		      "---------- EXCITATION ----------" );

		tmp5= TA* xpr5;
		tmp4= TA* xpr4;

		if( ixtyp == 4)
		{
		  tmp1= xpr1/ wlam;
		  tmp2= xpr2/ wlam;
		  tmp3= xpr3/ wlam;
		  tmp6= xpr6/( wlam* wlam);

		  fprintf( output_fp, "\n"
		      "                                  "
		      "    CURRENT SOURCE\n"
		      "                     -- POSITION (METERS) -- "
		      "      ORIENTATION (DEG)\n"
		      "                     X          Y          Z "
		      "      ALPHA        BETA   DIPOLE MOMENT\n"
		      "               %10.5LF %10.5LF %10.5LF "
		      " %7.2LF     %7.2LF    %8.3LF",
		      xpr1, xpr2, xpr3, xpr4, xpr5, xpr6 );
		}
		else
		{
		  tmp1= TA* xpr1;
		  tmp2= TA* xpr2;
		  tmp3= TA* xpr3;
		  tmp6= xpr6;

		  if( iptflg <= 0)
		    fprintf( output_fp,
			"\n  PLANE WAVE - THETA: %7.2LF deg, PHI: %7.2LF deg,"
			" ETA=%7.2LF DEG, TYPE - %s  AXIAL RATIO: %6.3LF",
			xpr1, xpr2, xpr3, hpol[ixtyp-1], xpr6 );

		} /* if( ixtyp == 4) */

	      } /* if( (ixtyp  != 0) && (ixtyp <= 4) ) */

	      /* fills e field right-hand matrix */
	      etmns( tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, ixtyp, cur);

	      /* matrix solving  (netwk calls solves) */
	      if( (nonet != 0) && (inc <= 1) )
	      {
		fprintf( output_fp, "\n\n\n"
		    "                                            "
		    "---------- NETWORK DATA ----------" );

		itmp3=0;
		itmp1= ntyp[0];

		for( i = 0; i < 2; i++ )
		{
		  if( itmp1 == 3)
		    itmp1=2;

		  if( itmp1 == 2)
		    fprintf( output_fp, "\n"
			"  -- FROM -  --- TO --      TRANSMISSION LINE       "
			" --------- SHUNT ADMITTANCES (MHOS) ---------   LINE\n"
			"  TAG   SEG  TAG   SEG    IMPEDANCE      LENGTH    "
			" ----- END ONE -----      ----- END TWO -----   TYPE\n"
			"  No:   No:  No:   No:         OHMS      METERS     "
			" REAL      IMAGINARY      REAL      IMAGINARY" );
		  else
		    if( itmp1 == 1)
		      fprintf( output_fp, "\n"
			  "  -- FROM -  --- TO --            --------"
			  " ADMITTANCE MATRIX ELEMENTS (MHOS) ---------\n"
			  "  TAG   SEG  TAG   SEG   ----- (ONE,ONE) ------  "
			  " ----- (ONE,TWO) -----   ----- (TWO,TWO) -------\n"
			  "  No:   No:  No:   No:      REAL      IMAGINARY     "
			  " REAL     IMAGINARY       REAL      IMAGINARY" );

		  for( j = 0; j < nonet; j++)
		  {
		    itmp2= ntyp[j];

		    if( (itmp2/itmp1) != 1 )
		      itmp3 = itmp2;
		    else
		    {
		      int idx4, idx5;

		      itmp4= iseg1[j];
		      itmp5= iseg2[j];
		      idx4 = itmp4-1;
		      idx5 = itmp5-1;

		      if( (itmp2 >= 2) && (x11i[j] <= 0.) )
		      {
			long double xx, yy, zz;

			xx = x[idx5]- x[idx4];
			yy = y[idx5]- y[idx4];
			zz = z[idx5]- z[idx4];
			x11i[j]= wlam* sqrtl( xx*xx + yy*yy + zz*zz );
		      }

		      fprintf( output_fp, "\n"
			  " %4d %5d %4d %5d  %11.4LE %11.4LE  "
			  "%11.4LE %11.4LE  %11.4LE %11.4LE %s",
			  itag[idx4], itmp4, itag[idx5], itmp5,
			  x11r[j], x11i[j], x12r[j], x12i[j],
			  x22r[j], x22i[j], pnet[itmp2-1] );

		    } /* if(( itmp2/ itmp1) == 1) */

		  } /* for( j = 0; j < nonet; j++) */

		  if( itmp3 == 0)
		    break;

		  itmp1= itmp3;

		} /* for( j = 0; j < nonet; j++) */

	      } /* if( (nonet != 0) && (inc <= 1) ) */

	      if( (inc > 1) && (iptflg > 0) )
		nprint=1;

	      netwk( cm, &cm[ib11], &cm[ic11], &cm[id11], ip, cur );
	      ntsol=1;

	      if( iped != 0)
	      {
		itmp1= 4*( mhz-1);

		fnorm[itmp1  ]= creall( zped);
		fnorm[itmp1+1]= cimagl( zped);
		fnorm[itmp1+2]= cabsl( zped);
		fnorm[itmp1+3]= cang( zped);

		if( iped != 2 )
		{
		  if( fnorm[itmp1+2] > zpnorm)
		    zpnorm= fnorm[itmp1+2];
		}

	      } /* if( iped != 0) */

	      /* printing structure currents */
	      if( n != 0)
	      {
		if( iptflg != -1)
		{
		  if( iptflg <= 0)
		  {
		    fprintf( output_fp, "\n\n\n"
			"                           "
			"-------- CURRENTS AND LOCATION --------\n"
			"                                  "
			"DISTANCES IN WAVELENGTHS" );

		    fprintf( output_fp,	"\n\n"
			"   SEG  TAG    COORDINATES OF SEGM CENTER     SEGM"
			"    ------------- CURRENT (AMPS) -------------\n"
			"   No:  No:       X         Y         Z      LENGTH"
			"     REAL      IMAGINARY    MAGN        PHASE" );
		  }
		  else
		  {
		    if( (iptflg != 3) && (inc <= 1) )
		      fprintf( output_fp, "\n\n\n"
			  "             -------- RECEIVING PATTERN PARAMETERS --------\n"
			  "                      ETA: %7.2LF DEGREES\n"
			  "                      TYPE: %s\n"
			  "                      AXIAL RATIO: %6.3LF\n\n"
			  "            THETA     PHI      ----- CURRENT ----    SEG\n"
			  "            (DEG)    (DEG)     MAGNITUDE    PHASE    No:",
			  xpr3, hpol[ixtyp-1], xpr6 );

		  } /* if( iptflg <= 0) */

		} /* if( iptflg != -1) */

		ploss=0.;
		itmp1=0;
		jump= iptflg+1;

		for( i = 0; i < n; i++ )
		{
		  curi= cur[i]* wlam;
		  cmag= cabsl( curi);
		  ph= cang( curi);

		  if( (nload != 0) && (fabsl(creall(zarray[i])) >= 1.e-20) )
		    ploss += 0.5* cmag* cmag* creall( zarray[i])* si[i];

		  if( jump == 0)
		    continue;

		  if( jump > 0 )
		  {
		    if( (iptag != 0) && (itag[i] != iptag) )
		      continue;

		    itmp1++;
		    if( (itmp1 < iptagf) || (itmp1 > iptagt) )
		      continue;

		    if( iptflg != 0)
		    {
		      if( iptflg >= 2 )
		      {
			fnorm[inc-1]= cmag;
			isave= (i+1);
		      }

		      if( iptflg != 3)
		      {
			fprintf( output_fp, "\n"
			    "          %7.2LF  %7.2LF   %11.4LE  %7.2LF  %5d",
			    xpr1, xpr2, cmag, ph, i+1 );

			continue;
		      }

		    } /* if( iptflg != 0) */

		  } /* if( jump > 0 ) */
		  else
		  {
		    fprintf( output_fp, "\n"
			" %5d %4d %9.4LF %9.4LF %9.4LF %9.5LF"
			" %11.4LE %11.4LE %11.4LE %8.3LF",
			i+1, itag[i], x[i], y[i], z[i], si[i],
			creall(curi), cimagl(curi), cmag, ph );

		    if( iplp1 != 1 )
		      continue;

		    if( iplp2 == 1)
		      fprintf( plot_fp, "%12.4LE %12.4LE\n", creall(curi), cimagl(curi) );
		    else
		      if( iplp2 == 2)
			fprintf( plot_fp, "%12.4LE %12.4LE\n", cmag, ph );
		  }

		} /* for( i = 0; i < n; i++ ) */

		if( iptflq != -1)
		{
		  fprintf( output_fp, "\n\n\n"
		      "                                  "
		      "------ CHARGE DENSITIES ------\n"
		      "                                  "
		      "   DISTANCES IN WAVELENGTHS\n\n"
		      "   SEG   TAG    COORDINATES OF SEG CENTER     SEG        "
		      "  CHARGE DENSITY (COULOMBS/METER)\n"
		      "   No:   No:     X         Y         Z       LENGTH   "
		      "  REAL      IMAGINARY     MAGN        PHASE" );

		  itmp1 = 0;
		  fr = 1.e-6/fmhz;

		  for( i = 0; i < n; i++ )
		  {
		    if( iptflq != -2 )
		    {
		      if( (iptaq != 0) && (itag[i] != iptaq) )
			continue;

		      itmp1++;
		      if( (itmp1 < iptaqf) || (itmp1 > iptaqt) )
			continue;

		    } /* if( iptflq == -2) */

		    curi= fr* cmplx(- bii[i], bir[i]);
		    cmag= cabsl( curi);
		    ph= cang( curi);

		    fprintf( output_fp, "\n"
			" %5d %4d %9.4LF %9.4LF %9.4LF %9.5LF"
			" %11.4LE %11.4LE %11.4LE %9.3LF",
			i+1, itag[i], x[i], y[i], z[i], si[i],
			creall(curi), cimagl(curi), cmag, ph );

		  } /* for( i = 0; i < n; i++ ) */

		} /* if( iptflq != -1) */

	      } /* if( n != 0) */

	      if( m != 0)
	      {
		fprintf( output_fp, "\n\n\n"
		    "                                      "
		    " --------- SURFACE PATCH CURRENTS ---------\n"
		    "                                                "
		    " DISTANCE IN WAVELENGTHS\n"
		    "                                                "
		    " CURRENT IN AMPS/METER\n\n"
		    "                                 ---------"
		    " SURFACE COMPONENTS --------   "
		    " ---------------- RECTANGULAR COMPONENTS ----------------\n"
		    "  PCH   --- PATCH CENTER ---     TANGENT VECTOR 1    "
		    " TANGENT VECTOR 2    ------- X ------    ------- Y ------   "
		    " ------- Z ------\n  No:    X       Y       Z       MAG.    "
		    "   PHASE     MAG.       PHASE    REAL   IMAGINARY    REAL  "
		    " IMAGINARY    REAL   IMAGINARY" );

		j= n-3;
		itmp1= -1;

		for( i = 0; i < m; i++ )
		{
		  j += 3;
		  itmp1++;
		  ex= cur[j];
		  ey= cur[j+1];
		  ez= cur[j+2];
		  eth= ex* t1x[itmp1]+ ey* t1y[itmp1]+ ez* t1z[itmp1];
		  eph= ex* t2x[itmp1]+ ey* t2y[itmp1]+ ez* t2z[itmp1];
		  ethm= cabsl( eth);
		  etha= cang( eth);
		  ephm= cabsl( eph);
		  epha= cang( eph);

		  fprintf( output_fp, "\n"
		      " %4d %7.3LF %7.3LF %7.3LF %11.4LE %8.2LF %11.4LE %8.2LF"
		      " %9.2LE %9.2LE %9.2LE %9.2LE %9.2LE %9.2LE",
		      i+1, px[itmp1], py[itmp1], pz[itmp1],
		      ethm, etha, ephm, epha, creall(ex), cimagl(ex),
		      creall(ey), cimagl(ey), creall(ez), cimagl(ez) );

		  if( iplp1 != 1)
		    continue;

		  if( iplp3 == 1)
		    fprintf( plot_fp, "%12.4LE %12.4LE\n", creall(ex), cimagl(ex) );
		  if( iplp3 == 2)
		    fprintf( plot_fp, "%12.4LE %12.4LE\n", creall(ey), cimagl(ey) );
		  if( iplp3 == 3)
		    fprintf( plot_fp, "%12.4LE %12.4LE\n", creall(ez), cimagl(ez) );
		  if( iplp3 == 4)
		    fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE %12.4LE %12.4LE %12.4LE\n",
			creall(ex),cimagl(ex),creall(ey),cimagl(ey),creall(ez),cimagl(ez) );

		} /* for( i=0; i<m; i++ ) */

	      } /* if( m != 0) */

	      if( (ixtyp == 0) || (ixtyp == 5) )
	      {
		tmp1= pin- pnls- ploss;
		tmp2= 100.* tmp1/ pin;

		fprintf( output_fp, "\n\n\n"
		    "                               "
		    "---------- POWER BUDGET ---------\n"
		    "                               "
		    "INPUT POWER   = %11.4LE Watts\n"
		    "                               "
		    "RADIATED POWER= %11.4LE Watts\n"
		    "                               "
		    "STRUCTURE LOSS= %11.4LE Watts\n"
		    "                               "
		    "NETWORK LOSS  = %11.4LE Watts\n"
		    "                               "
		    "EFFICIENCY    = %7.2LF Percent",
		    pin, tmp1, ploss, pnls, tmp2 );

	      } /* if( (ixtyp == 0) || (ixtyp == 5) ) */

	      igo = 4;

	      if( ncoup > 0)
		couple( cur, wlam );

	      if( iflow == 7)
	      {
		if( (ixtyp > 0) && (ixtyp < 4) )
		{
		  nthic++;
		  inc++;
		  xpr1 += xpr4;

		  if( nthic <= nthi )
		    continue; /* continue excitation loop */

		  nthic=1;
		  xpr1= thetis;
		  xpr2= xpr2+ xpr5;
		  nphic++;

		  if( nphic <= nphi )
		    continue; /* continue excitation loop */

		  break;

		} /* if( (ixtyp >= 1) && (ixtyp <= 3) ) */

		if( nfrq != 1)
		{
		  jmp_floop = TRUE;
		  break; /* continue the freq loop */
		}

		fprintf( output_fp, "\n\n\n" );
		jmp_iloop = TRUE;

		break; /* continue card input loop */

	      } /*if( iflow == 7) */


	      case 4: /* label_71 */
	      igo = 5;

	      /* label_72 */
	      case 5: /* near field calculation */

	      if( near != -1)
	      {
		nfpat();

		if( mhz == nfrq)
		  near=-1;

		if( nfrq == 1)
		{
		  fprintf( output_fp, "\n\n\n" );
		  jmp_iloop = TRUE;
		  break; /* continue card input loop */
		}

	      } /* if( near != -1) */


	      /* label_78 */
	      case 6: /* standard far field calculation */

	      if( ifar != -1)
	      {
		pinr= pin;
		pnlr= pnls;
		rdpat();
	      }

	      if( (ixtyp == 0) || (ixtyp >= 4) )
	      {
		if( mhz == nfrq )
		  ifar=-1;

		if( nfrq != 1)
		{
		  jmp_floop = TRUE;
		  break;
		}

		fprintf( output_fp, "\n\n\n" );
		jmp_iloop = TRUE;
		break;

	      } /* if( (ixtyp == 0) || (ixtyp >= 4) ) */

	      nthic++;
	      inc++;
	      xpr1 += xpr4;

	      if( nthic <= nthi )
		continue; /* continue excitation loop */

	      nthic = 1;
	      xpr1  = thetis;
	      xpr2 += xpr5;
	      nphic++;

	      if( nphic > nphi )
		break;

	    } /* do (l_54) */
	    while( TRUE );

	    /* jump to freq. or input loop */
	    if( jmp_iloop )
	      break;

	    if( jmp_floop )
	      continue;

	    nphic = 1;
	    xpr2  = phiss;

	    /* normalized receiving pattern printed */
	    if( iptflg >= 2)
	    {
	      itmp1= nthi* nphi;

	      tmp1= fnorm[0];
	      for( j = 1; j < itmp1; j++ )
		if( fnorm[j] > tmp1)
		  tmp1= fnorm[j];

	      fprintf( output_fp, "\n\n\n"
		  "                     "
		  "---- NORMALIZED RECEIVING PATTERN ----\n"
		  "                      "
		  "NORMALIZATION FACTOR: %11.4LE\n"
		  "                      "
		  "ETA: %7.2LF DEGREES\n"
		  "                      "
		  "TYPE: %s\n"
		  "                      AXIAL RATIO: %6.3LF\n"
		  "                      SEGMENT No: %d\n\n"
		  "                      "
		  "THETA     PHI       ---- PATTERN ----\n"
		  "                      "
		  "(DEG)    (DEG)       DB     MAGNITUDE",
		  tmp1, xpr3, hpol[ixtyp-1], xpr6, isave );

	      for( j = 0; j < nphi; j++ )
	      {
		itmp2= nthi*j;

		for( i = 0; i < nthi; i++ )
		{
		  itmp3= i + itmp2;

		  if( itmp3 < itmp1)
		  {
		    tmp2= fnorm[itmp3]/ tmp1;
		    tmp3= db20( tmp2);

		    fprintf( output_fp, "\n"
			"                    %7.2LF  %7.2LF   %7.2LF  %11.4LE",
			xpr1, xpr2, tmp3, tmp2 );

		    xpr1 += xpr4;
		  }

		} /* for( i = 0; i < nthi; i++ ) */

		xpr1= thetis;
		xpr2 += xpr5;

	      } /* for( j = 0; j < nphi; j++ ) */

	      xpr2= phiss;

	    } /* if( iptflg >= 2) */

	    if( mhz == nfrq)
	      ifar=-1;

	    if( nfrq == 1)
	    {
	      fprintf( output_fp, "\n\n\n" );
	      jmp_iloop = TRUE;
	      break; /* continue card input loop */
	    }

	  } /*** do (frequency loop) (l_42) ***/
	  while( (++mhz <= nfrq) );

	  /* Jump to card input loop */
	  if( jmp_iloop )
	    break;

	  if( iped != 0)
	  {
	    int iss;

	    if( nvqd > 0)
	      iss = ivqd[nvqd-1];
	    else
	      iss = isant[nsant-1];

	    fprintf( output_fp, "\n\n\n"
		"                            "
		" -------- INPUT IMPEDANCE DATA --------\n"
		"                                     "
		" SOURCE SEGMENT No: %d\n"
		"                                  "
		" NORMALIZATION FACTOR:%12.5LE\n\n"
		"              ----------- UNNORMALIZED IMPEDANCE ----------  "
		"  ------------ NORMALIZED IMPEDANCE -----------\n"
		"      FREQ    RESISTANCE    REACTANCE    MAGNITUDE    PHASE  "
		"  RESISTANCE    REACTANCE    MAGNITUDE    PHASE\n"
		"       MHz       OHMS         OHMS         OHMS     DEGREES  "
		"     OHMS         OHMS         OHMS     DEGREES",
		iss, zpnorm );

	    itmp1= nfrq;
	    if( ifrq == 0)
	      tmp1= fmhz-( nfrq-1)* delfrq;
	    else
	      if( ifrq == 1)
		tmp1= fmhz/( powl(delfrq, (nfrq-1)) );

	    for( i = 0; i < itmp1; i++ )
	    {
	      itmp2= 4*i;
	      tmp2= fnorm[itmp2  ]/ zpnorm;
	      tmp3= fnorm[itmp2+1]/ zpnorm;
	      tmp4= fnorm[itmp2+2]/ zpnorm;
	      tmp5= fnorm[itmp2+3];

	      fprintf( output_fp, "\n"
		  " %9.3LF   %11.4LE  %11.4LE  %11.4LE  %7.2LF  "
		  " %11.4LE  %11.4LE  %11.4LE  %7.2LF",
		  tmp1, fnorm[itmp2], fnorm[itmp2+1], fnorm[itmp2+2],
		  fnorm[itmp2+3], tmp2, tmp3, tmp4, tmp5 );

	      if( ifrq == 0)
		tmp1 += delfrq;
	      else
		if( ifrq == 1)
		  tmp1 *= delfrq;

	    } /* for( i = 0; i < itmp1; i++ ) */

	    fprintf( output_fp, "\n\n\n" );

	  } /* if( iped != 0) */

	  nfrq=1;
	  mhz=1;

      } /* switch( igox ) */

    } /* while( ! next_job ): Main input section (l_14) */

  } /* while(TRUE): Main execution loop (l_1) */

  return(0);

} /* end of main() */

/*-----------------------------------------------------------------------*/

/* prnt sets up the print formats for impedance loading */
void prnt( int in1, int in2, int in3, long double fl1, long double fl2,
    long double fl3, long double fl4, long double fl5, long double fl6, char *ia, int ichar )
{
  /* record to be output and buffer used to make it */
  char record[101+ichar*4], buf[15];
  int in[3], i1, i;
  long double fl[6];

  in[0]= in1;
  in[1]= in2;
  in[2]= in3;
  fl[0]= fl1;
  fl[1]= fl2;
  fl[2]= fl3;
  fl[3]= fl4;
  fl[4]= fl5;
  fl[5]= fl6;

  /* integer format */
  i1=0;
  strcpy( record, "\n " );

  if( (in1 == 0) && (in2 == 0) && (in3 == 0) )
  {
    strcat( record, " ALL" );
    i1=1;
  }

  for( i = i1; i < 3; i++ )
  {
    if( in[i] == 0)
      strcat( record, "     " );
    else
    {
      sprintf( buf, "%5d", in[i] );
      strcat( record, buf );
    }
  }

  /* floating point format */
  for( i = 0; i < 6; i++ )
  {
    if( fabsl( fl[i]) >= 1.0e-20 )
    {
      sprintf( buf, " %11.4LE", fl[i] );
      strcat( record, buf );
    }
    else
      strcat( record, "            " );
  }

  strcat( record, "   " );
  strcat( record, ia );
  fprintf( output_fp, "%s", record );

  return;
}

/*-----------------------------------------------------------------------*/

static void sig_handler( int signal )
{
  switch( signal )
  {
    case SIGINT :
      fprintf( stderr, "\n%s\n", "nec2c: exiting via user interrupt" );
      exit( signal );

    case SIGSEGV :
      fprintf( stderr, "\n%s\n", "nec2c: segmentation fault" );
      exit( signal );

    case SIGFPE :
      fprintf( stderr, "\n%s\n", "nec2c: floating point exception" );
      exit( signal );

    case SIGABRT :
      fprintf( stderr, "\n%s\n", "nec2c: abort signal received" );
      exit( signal );

    case SIGTERM :
      fprintf( stderr, "\n%s\n", "nec2c: termination request received" );

      stop( signal );
  }

} /* end of sig_handler() */

/*------------------------------------------------------------------------*/

