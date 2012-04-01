#ifndef	NEC2C_H
#define	NEC2C_H 1

#include <complex.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>

#ifndef	TRUE
#define	TRUE	1
#endif

#ifndef	FALSE
#define	FALSE	0
#endif

/* commonly used complex constants */
#define	CPLX_00	(0.0+0.0fj)
#define	CPLX_01	(0.0+1.0fj)
#define	CPLX_10	(1.0+0.0fj)
#define	CPLX_11	(1.0+1.0fj)

/* common constants */
#define	PI	3.141592654
#define	POT	1.570796327
#define	TP	6.283185308
#define	PTP	.6283185308
#define	TPJ	(0.0+6.283185308fj)
#define PI8	25.13274123
#define PI10	31.41592654
#define	TA	1.745329252E-02
#define	TD	57.29577951
#define	ETA	376.73
#define	CVEL	299.8
#define	RETA	2.654420938E-3
#define	TOSP	1.128379167
#define ACCS	1.E-12
#define	SP	1.772453851
#define	FPI	12.56637062
#define	CCJ	(0.0-0.01666666667fj)
#define	CONST1	(0.0+4.771341189fj)
#define	CONST2	4.771341188
#define	CONST3	(0.0-29.97922085fj)
#define	CONST4	(0.0+188.365fj)
#define	GAMMA	.5772156649
#define C1	-.02457850915
#define C2	.3674669052
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
#define MAXH	20
#define CRIT	1.0E-4
#define NM	131072
#define NTS	4
#define	SMIN	1.e-3

/* Replaces the "10000" limit used to */
/* identify segment/patch connections */
#define	PCHCON  100000

/* carriage return and line feed */
#define	CR	0x0d
#define	LF	0x0a

/* max length of a line read from input file */
#define	LINE_LEN	132/* version of fortran source for the -v option */
#define		version "nec2c v0.1"

/* Function prototypes */
int 	main(int argc, char **argv);
void 	arc(int itg, int ns, long double rada, long double ang1,
	long double ang2, long double rad);
void 	blckot(complex long double *ar, int nunit, int ix1, int ix2, int nblks, int neof);
void 	blckin(complex long double *ar, int nunit, int ix1, int ix2, int nblks, int neof);
void 	cabc(complex long double *curx);
void 	cmset(int nrow, complex long double *cm, long double rkhx, int iexkx);
void 	cmss(int j1, int j2, int im1, int im2,
	complex long double *cm, int nrow, int itrp);
void 	cmsw(int j1, int j2, int i1, int i2, complex long double *cm,
	complex long double *cw, int ncw, int nrow, int itrp);
void 	cmws(int j, int i1, int i2, complex long double *cm, int nr,
	complex long double *cw, int nw, int itrp);
void 	cmww(int j, int i1, int i2, complex long double *cm, int nr,
	complex long double *cw, int nw, int itrp);
void 	conect(int ignd);
void 	couple(complex long double *cur, long double wlam);
void 	datagn(void);
long double db10(long double x);
long double db20(long double x);
void 	efld(long double xi, long double yi, long double zi, long double ai, int ij);
void 	eksc(long double s, long double z, long double rh, long double xk, int ij,
	complex long double *ezs, complex long double *ers, complex long double *ezc,
	complex long double *erc, complex long double *ezk, complex long double *erk);
void 	ekscx(long double bx, long double s, long double z, long double rhx, long double xk,
	int ij, int inx1, int inx2, complex long double *ezs,
	complex long double *ers, complex long double *ezc, complex long double *erc,
	complex long double *ezk, complex long double *erk);
void 	etmns(long double p1, long double p2, long double p3, long double p4, long double p5,
	long double p6, int ipr, complex long double *e);
void 	factr(int n, complex long double *a, int *ip, int ndim);
void 	factrs(int np, int nrow, complex long double *a, int *ip);
complex long double fbar(complex long double p);
void 	fblock(int nrow, int ncol, int imax, int ipsym);
void 	ffld(long double thet, long double phi,
        complex long double *eth, complex long double *eph);
void 	fflds(long double rox, long double roy, long double roz, complex long double *scur,
	complex long double *ex, complex long double *ey, complex long double *ez);
void 	gf(long double zk, long double *co, long double *si);
void 	gfld(long double rho, long double phi, long double rz, complex long double *eth,
	complex long double *epi, complex long double *erd, complex long double ux, int ksymp);
void 	gh(long double zk, long double *hr, long double *hi);
void 	gwave(complex long double *erv, complex long double *ezv,
	complex long double *erh, complex long double *ezh, complex long double *eph);
void 	gx(long double zz, long double rh, long double xk,
	complex long double *gz, complex long double *gzp);
void 	gxx(long double zz, long double rh, long double a, long double a2, long double xk,
	int ira, complex long double *g1, complex long double *g1p, complex long double *g2,
	complex long double *g2p, complex long double *g3, complex long double *gzp);
void 	helix(long double s, long double hl, long double a1, long double b1,
	long double a2,long double b2, long double rad, int ns, int itg);
void 	hfk(long double el1, long double el2, long double rhk,
	long double zpkx, long double *sgr, long double *sgi);
void 	hintg(long double xi, long double yi, long double zi);
void 	hsfld(long double xi, long double yi, long double zi, long double ai);
void 	hsflx(long double s, long double rh, long double zpx, complex long double *hpk,
	complex long double *hps, complex long double *hpc);
void 	intrp(long double x, long double y, complex long double *f1,
	complex long double *f2, complex long double *f3, complex long double *f4);
void 	intx(long double el1, long double el2, long double b, int ij,
	long double *sgr, long double *sgi);
int 	isegno(int itagi, int mx);
void 	lfactr(complex long double *a, int nrow, int ix1, int ix2, int *ip);
void 	load(int *ldtyp, int *ldtag, int *ldtagf, int *ldtagt,
	long double *zlr, long double *zli, long double *zlc);
void 	lunscr(complex long double *a, int nrow, int nop,
	int *ix, int *ip, int iu2, int iu3, int iu4);
void 	move(long double rox, long double roy, long double roz, long double xs,
	long double ys, long double zs, int its, int nrpt, int itgi);
void 	nefld(long double xob, long double yob, long double zob, complex long double *ex,
	complex long double *ey, complex long double *ez);
void 	netwk(complex long double *cm, complex long double *cmb, complex long double *cmc,
	complex long double *cmd, int *ip, complex long double *einc);
void 	nfpat(void);
void 	nhfld(long double xob, long double yob, long double zob, complex long double *hx,
	complex long double *hy, complex long double *hz);
void 	patch(int nx, int ny, long double ax1, long double ay1, long double az1,
	long double ax2, long double ay2, long double az2, long double ax3, long double ay3,
	long double az3, long double ax4, long double ay4, long double az4);
void 	subph(int nx, int ny);
void 	pcint(long double xi, long double yi, long double zi, long double cabi,
	long double sabi, long double salpi, complex long double *e);
void 	prnt(int in1, int in2, int in3, long double fl1, long double fl2,
	long double fl3, long double fl4, long double fl5,
	long double fl6, char *ia, int ichar);
void 	qdsrc(int is, complex long double v, complex long double *e);
void 	rdpat(void);
void 	readgm(char *gm, int *i1, int *i2, long double *x1, long double *y1,
	long double *z1, long double *x2, long double *y2, long double *z2, long double *rad);
void 	readmn(char *gm, int *i1, int *i2, int *i3, int *i4, long double *f1,
	long double *f2, long double *f3, long double *f4, long double *f5, long double *f6);
void 	reflc(int ix, int iy, int iz, int itx, int nop);
void 	rom2(long double a, long double b, complex long double *sum, long double dmin);
void 	sbf(int i, int is, long double *aa, long double *bb, long double *cc);
void 	sflds(long double t, complex long double *e);
void 	solgf(complex long double *a, complex long double *b, complex long double *c,
	complex long double *d, complex long double *xy, int *ip, int np, int n1,
	int n, int mp, int m1, int m, int n1c, int n2c, int n2cz);
void 	solve(int n, complex long double *a, int *ip, complex long double *b, int ndim);
void 	solves(complex long double *a, int *ip, complex long double *b, int neq,
	int nrh, int np, int n, int mp, int m);
void 	tbf(int i, int icap);
void 	test(long double f1r, long double f2r, long double *tr, long double f1i,
	long double f2i, long double *ti, long double dmin);
void 	trio(int j);
void 	unere(long double xob, long double yob, long double zob);
void 	wire(long double xw1, long double yw1, long double zw1,
	long double xw2, long double yw2,
	long double zw2, long double rad, long double rdel,
	long double rrad, int ns, int itg);
complex long double zint(long double sigl, long double rolam);
int 	min(int a, int b);
void 	usage(void);
void 	abort_on_error(int why);
complex long double cmplx(long double a, long double j);
void 	secnds(long double *x);
int 	stop(int flag);
int 	load_line(char *buff, FILE *pfile);
void	mem_alloc( void **ptr, int req );
void	mem_realloc( void **ptr, int req );
void	free_ptr( void **ptr );
void 	somnec(long double epr, long double sig, long double fmhz);
void 	bessel(complex long double z, complex long double *j0, complex long double *j0p);
void 	evlua(complex long double *erv, complex long double *ezv,
	complex long double *erh, complex long double *eph);
void 	gshank(complex long double start, complex long double dela,
	complex long double *sum, int nans, complex long double *seed,
	int ibk, complex long double bk, complex long double delb);
void 	hankel(complex long double z, complex long double *h0, complex long double *h0p);
void 	lambda(long double t, complex long double *xlam, complex long double *dxlam);
void 	rom1(int n, complex long double *sum, int nx);
void 	saoa( long double t, complex long double *ans);
long double cang(complex long double z);

#endif
