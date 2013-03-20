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





From: <BURKE_at_email.domain.hidden>
Date: Wed Apr 10 1996 - 22:17:10 EDT

Re: NEC benchmarks, and a way to speedup your code (maybe)

In running benchmarks on different platforms with NEC codes of different
origins, anomalous timings could result from the transposed (or not)
storage of the matrix. This was discussed in the July 1995 ACES
Newsletter. Since many NEC users probably did not see that, there will
be versions of NEC around with and without this do-it-yourself change
that can greatly speedup the solution. All codes sent out from LLNL
since May 1995 have the matrix un-transposed for in-core solution.

In NEC-2, 3 and 4 the matrix is stored in transposed form (column,row)
for convenience in the out-of-core solution. Transposed storage is
unnecessary when the matrix fits into RAM, and will slow down the matrix
solution since inner loops in subroutines FACTR and SOLVE are changing
the second index in the array. The difference between transposed and
normal storage was not very large on older computers, but can be very
significant on PCs and workstations that depend on cache storage. On a
VAX 6330 the matrix-factor time with matrix transposed was slower than
with normal storage by a factor of 1.2. However, on an Alpha 3000/400
the factor was 2.4 and on a Mac 8100 it was 3.5 to 4.

Fortunately it is easy to fix this problem if you have a FORTRAN
compiler. Rather than changing the entire matrix fill and NGF, you can
simply re-transpose the matrix at the beginning of subroutine FACTR.
Instructions for doing this were given in the ACES Newsletter, but
FACTR and SOLVE routines with the changes made are included below.
These should just replace the same routines in NEC-2, 3 or 4 after the
parameter MAXSEG has been set correctly.

(A confusing factor if you are not used to archaic FORTRAN: 2*MAXSEG
sets the size of the complex array D in COMMON/SCRATM/ in FACTR and
SOLVE. COMMON/SCRATM/ is contained in several other routines in NEC,
and the lengths must match. However, in subroutine RDPAT the array is
real, so the dimension would be 4*MAXSEG. If the array in your FACTR
and SOLVE has a number for a dimension, then MAXSEG should be half that
number. If your code gets the dimension from an INCLUDE file, you
should use that INCLUDE.)

Re-transposing the matrix takes a little time, proportional to N**2,
but saves time proportional to N**3. The amount of time wasted is
negligible, and the risk of breaking something if you try to undo the
transposed storage throughout the code is great. I do not know whether
this change has gotten into the codes on ftp.netcom.com. At least for
benchmarks it should be checked that any codes compared are the same,
re-transposed for speed or not for "classic" NEC.

Jerry Burke
LLNL

FACTR and SOLVE routines with matrix re-transposed:

      SUBROUTINE FACTR (N,A,IP,NDIM)
C
C Purpose:
C FACTR computes the LU decomposition of the matrix in A. The
C algorithm is described on pages 411-416 OF A. Ralston--A First
C Course in Numerical Analysis. Comments below refer to comments in
C Ralston's text. (MATRIX TRANSPOSED.)
C
      PARAMETER (MAXSEG=300)
      COMPLEX A,D,ARJ
      DIMENSION A(NDIM,NDIM), IP(NDIM)
      COMMON /SCRATM/ D(2*MAXSEG)
      INTEGER R,RM1,RP1,PJ,PR
C
C Un-transpose the matrix for Gauss elimination
C
      DO 12 I=2,N
         DO 11 J=1,I-1
            ARJ=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=ARJ
11 CONTINUE
12 CONTINUE
      IFLG=0
      DO 9 R=1,N
C
C STEP 1
C
      DO 1 K=1,N
      D(K)=A(K,R)
1 CONTINUE
C
C STEPS 2 AND 3
C
      RM1=R-1
      IF (RM1.LT.1) GO TO 4
      DO 3 J=1,RM1
      PJ=IP(J)
      ARJ=D(PJ)
      A(J,R)=ARJ
      D(PJ)=D(J)
      JP1=J+1
      DO 2 I=JP1,N
      D(I)=D(I)-A(I,J)*ARJ
2 CONTINUE
3 CONTINUE
4 CONTINUE
C
C STEP 4
C
      DMAX=REAL(D(R)*CONJG(D(R)))
      IP(R)=R
      RP1=R+1
      IF (RP1.GT.N) GO TO 6
      DO 5 I=RP1,N
      ELMAG=REAL(D(I)*CONJG(D(I)))
      IF (ELMAG.LT.DMAX) GO TO 5
      DMAX=ELMAG
      IP(R)=I
5 CONTINUE
6 CONTINUE
      IF (DMAX.LT.1.E-10) IFLG=1
      PR=IP(R)
      A(R,R)=D(PR)
      D(PR)=D(R)
C
C STEP 5
C
      IF (RP1.GT.N) GO TO 8
      ARJ=1./A(R,R)
      DO 7 I=RP1,N
      A(I,R)=D(I)*ARJ
7 CONTINUE
8 CONTINUE
      IF (IFLG.EQ.0) GO TO 9
      WRITE(3,10) R,DMAX
      IFLG=0
9 CONTINUE
      RETURN
C
10 FORMAT (' FACTR: PIVOT(',I3,')=',1PE16.8)
      END
      SUBROUTINE SOLVE (N,A,IP,B,NDIM)
C=======================================================================
C (C) Copyright 1992
C The Regents of the University of California. All rights reserved.
C=======================================================================
C
C SOLVE solves the matrix equation LU*X=B where L is a unit
C lower triangular matrix and U is an upper triangular matrix both
C of which are stored in A. The RHS vector B is input and the
C solution is returned through vector B.
C
      PARAMETER (MAXSEG=300)
      COMPLEX A,B,Y,SUM
      INTEGER PI
      COMMON /SCRATM/ Y(2*MAXSEG)
      DIMENSION A(NDIM,NDIM), IP(NDIM), B(NDIM)
C
C FORWARD SUBSTITUTION
C
      DO 3 I=1,N
      PI=IP(I)
      Y(I)=B(PI)
      B(PI)=B(I)
      IP1=I+1
      IF (IP1.GT.N) GO TO 2
      DO 1 J=IP1,N
      B(J)=B(J)-A(J,I)*Y(I)
1 CONTINUE
2 CONTINUE
3 CONTINUE
C
C BACKWARD SUBSTITUTION
C
      DO 6 K=1,N
      I=N-K+1
      SUM=(0.,0.)
      IP1=I+1
      IF (IP1.GT.N) GO TO 5
      DO 4 J=IP1,N
      SUM=SUM+A(I,J)*B(J)
4 CONTINUE
5 CONTINUE
      B(I)=(Y(I)-SUM)/A(I,I)
6 CONTINUE
      RETURN
      END
*/
#include "math_util.h"
#include <iostream>
using namespace std;

#include "nec_exception.h"
#include "matrix_algebra.h"
#include "nec_output.h"

#include "config.h"

#ifdef NEC_MATRIX_CHECK
//#define NEC_MATRIX_CHECK

void to_octave(nec_complex& x);
void to_octave(nec_complex& x)
{
	cout << real(x) << " + " << imag(x) << "I";
}

void to_octave(int& x);
void to_octave(int& x)
{
	cout << x;
}

void to_octave(nec_complex* a, int n, int ndim);
void to_octave(nec_complex* a, int n, int ndim)
{
	cout << "[";
	for (int row = 0; row < n; row++ )
	{
		int row_offset = row*ndim;
		for (int i = 0; i < n; i++ )
		{
			to_octave(a[i+row_offset]);
			if (i < n-1)
				cout << ", ";
		}
		if (row < n-1)
			cout << "; ";
	}
	cout << "];" << endl;
}

void to_octave(complex_array& a, int n, int ndim);
void to_octave(complex_array& a, int n, int ndim)
{
	to_octave(a.get_ptr(),n,ndim);
}

void to_octave(int* a, int n);
void to_octave(int* a, int n)
{
	cout << "[";
	for (int i = 0; i < n; i++ )
	{
		to_octave(a[i]);
		if (i < n-1)
			cout << ", ";
	}
	cout << "];" << endl;
}

void to_octave(int_array& a, int n);
void to_octave(int_array& a, int n)
{
	to_octave(a.get_ptr(),n);
}
#endif



#ifndef LAPACK
/*! \brief Solve The system of equations using Gaussian Elimination.
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using the Gauss-Doolittle algorithm 
	presented on pages 411-416 of A. Ralston -- a first course in 
	numerical analysis.
	
	Comments below refer to comments in Ralstons text.
	
	(matrix transposed.)
*/
void lu_decompose(nec_output_file& s_output, int n, complex_array& a, int_array& ip, int ndim)
{
	DEBUG_TRACE("lu_decompose_ge(" << n << "," << ndim << ")");
	
#ifdef NEC_MATRIX_CHECK
	// Debug output to try and figure out the LAPACK stuff
	cout << "a = ";
	to_octave(a,n,ndim);
#endif
	
	/* Allocate scratch memory */
	complex_array scm;
	scm.resize(n);
	
	/* Un-transpose the matrix for Gauss elimination */
	for (int i = 1; i < n; i++ )
	{
		int i_offset = i * ndim;
		int j_offset = 0;
		for (int j = 0; j < i; j++ )
		{
			nec_complex aij = a[i+j_offset];
			a[i+j_offset] = a[j+i_offset];
			a[j+i_offset] = aij;
			
			j_offset += ndim;
		}
	}
	
	bool iflg=false;
	/* step 1 */
	for (int r = 0; r < n; r++ )
	{
		int r_offset = r*ndim;
		
		for (int k = 0; k < n; k++ )
			scm[k]= a[k+r_offset];
		
		/* steps 2 and 3 */
		int rm1 = r;
		for (int j = 0; j < rm1; j++ )
		{
			int pj= ip[j]-1;
			nec_complex arj = scm[pj];
			a[j+r_offset]= arj;
			scm[pj]= scm[j];
			int jp1 = j+1;
		
			int j_offset = j*ndim;
			for (int i = jp1; i < n; i++ )
				scm[i] -= a[i+j_offset]* arj;
		
		}
		
		/* step 4 */
		nec_float dmax = norm(scm[r]);
		
		int rp1 = r+1;
		ip[r]= rp1;
		for (int i = rp1; i < n; i++ )
		{
			nec_float elmag = norm(scm[i]);
			if ( elmag >= dmax)
			{
				dmax = elmag;
				ip[r] = i+1;	// set the permute array element
			}
		}
		
		if ( dmax < 1.e-10)
			iflg=true;
		
		int pr = ip[r]-1;
		a[r+r_offset] = scm[pr];
		scm[pr] = scm[r];
		
		/* step 5 */
		if ( rp1 < n)
		{
			nec_complex arr = cplx_10() / a[r+r_offset];
			
			for (int i = rp1; i < n; i++ )
				a[i+r_offset]= scm[i]* arr;
		}
		
		if ( true == iflg )
		{
			s_output.string("\n  PIVOT(");
			s_output.integer(r);
			s_output.string(")= ");
			s_output.real(dmax);
			iflg=false;
		}	
	} /* for( r=0; r < n; r++ ) */
	

#ifdef NEC_MATRIX_CHECK
	cout << "solved = ";
	to_octave(a,n,ndim);

	cout << "ip = ";
	to_octave(ip,n);
#endif
	
}


/*! \brief
	Subroutine to solve the matrix equation lu*x=b where l is a unit
	lower triangular matrix and u is an upper triangular matrix both
	of which are stored in a.  the rhs vector b is input and the
	solution is returned through vector b.   (matrix transposed)
      
	  COMPLEX*16 A,B,Y,SUM
      INTEGER PI
      COMMON /SCRATM/ Y(2*MAXSEG)
      DIMENSION A(NDIM,NDIM), IP(NDIM), B(NDIM)
C
C     FORWARD SUBSTITUTION
C
      DO 3 I=1,N
      PI=IP(I)
      Y(I)=B(PI)
      B(PI)=B(I)
      IP1=I+1
      IF (IP1.GT.N) GO TO 2
      DO 1 J=IP1,N
      B(J)=B(J)-A(J,I)*Y(I)
1     CONTINUE
2     CONTINUE
3     CONTINUE
C
C     BACKWARD SUBSTITUTION
C
      DO 6 K=1,N
      I=N-K+1
      SUM=(0.,0.)
      IP1=I+1
      IF (IP1.GT.N) GO TO 5
      DO 4 J=IP1,N
      SUM=SUM+A(I,J)*B(J)
4     CONTINUE
5     CONTINUE
      B(I)=(Y(I)-SUM)/A(I,I)
6     CONTINUE
      RETURN
      END
	
*/
void solve( int n, complex_array& a, int_array& ip,
    complex_array& b, int ndim )
{
	DEBUG_TRACE("solve(" << n << "," << ndim << ")");
/*	
	We should use zgetrs from LAPACK to solve this.
*/	
	complex_array y(n);
	
	/* forward substitution */
	for (int i = 0; i < n; i++ )
	{
		int pivot_index = ip[i]-1;
		y[i]= b[pivot_index];
		b[pivot_index]= b[i];
		int ip1= i+1;
		
		int i_offset = i*ndim;
		for (int j = ip1; j < n; j++ )
			b[j] -= a[j+i_offset] * y[i];
	}
	
	/* backward substitution */
	for (int k = 0; k < n; k++ )
	{
		int i= n-k-1;
		nec_complex sum(cplx_00());
		int ip1= i+1;
		
		for (int j = ip1; j < n; j++ )
			sum += a[i+j*ndim]* b[j];
		
		b[i] = (y[i]- sum) / a[i+i*ndim];
	}
}


#else /*  LAPACK */
#warning("Using lapack")

extern "C"
{
//#include <atlas/atlas_enum.h>
#include <atlas/clapack.h>
}

// extern "C" void zgetrf_( int *M, int *N, void *A, int *LDA, int *IPIV, int *INFOp);
// //SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
// int lapack_zgetrf(int M, int N, void *A, int LDA, int *IPIV)
// {
// 	int info;
// 	zgetrf_(&M, &N, (void*)A, &LDA, IPIV, &info);
// 	return info;
// }
// 
// 
// #warning "Using LAPACK for Matrix Operations"

/*! Use lapack to perform LU decomposition
*/
void lu_decompose(nec_output_file& s_output,  int n, complex_array& a_in, int_array& ip, int ndim)
{
	DEBUG_TRACE("factor_lapack(" << n << "," << ndim << ")");
	ASSERT(n <= ndim);

#ifdef NEC_MATRIX_CHECK
	cout << "atlas_a = ";
	to_octave(a_in,n,ndim);
#endif

	/* Un-transpose the matrix for Gauss elimination */
	for (int i = 1; i < n; i++ )
	{
		int i_offset = i * ndim;
		int j_offset = 0;
		for (int j = 0; j < i; j++ )
		{
			nec_complex aij = a_in[i+j_offset];
			a_in[i+j_offset] = a_in[j+i_offset];
			a_in[j+i_offset] = aij;
			
			j_offset += ndim;
		}
	}
	
		
	/* Now call the LAPACK LU-Decomposition
	ZGETRF computes an LU factorization of a general M-by-N matrix A
	*  using partial pivoting with row interchanges.
	*
	*  The factorization has the form
	*     A = P * L * U
	*  where P is a permutation matrix, L is lower triangular with unit
	*  diagonal elements (lower trapezoidal if m > n), and U is upper
	*  triangular (upper trapezoidal if m < n).

	Arguments
	*  =========
	*
	*  M       (input) INTEGER
	*          The number of rows of the matrix A.  M >= 0.
	*
	*  N       (input) INTEGER
	*          The number of columns of the matrix A.  N >= 0.
	*
	*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
	*          On entry, the M-by-N matrix to be factored.
	*          On exit, the factors L and U from the factorization
	*          A = P*L*U; the unit diagonal elements of L are not stored.
	*
	*  LDA     (input) INTEGER
	*          The leading dimension of the array A.  LDA >= max(1,M).
	*
	*  IPIV    (output) INTEGER array, dimension (min(M,N))
	*          The pivot indices; for 1 <= i <= min(M,N), row i of the
	*          matrix was interchanged with row IPIV(i).
	*
	*  INFO    (output) INTEGER
	*          = 0:  successful exit
	*          < 0:  if INFO = -i, the i-th argument had an illegal value
	*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
	*                has been completed, but the factor U is exactly
	*                singular, and division by zero will occur if it is used
	*                to solve a system of equations.
	*/
	int info = clapack_zgetrf (CblasColMajor, n, n, 
		(void*) a_in.data(), ndim, ip.data());
	
	if (0 != info)
	{
		/*
			The factorization has been completed, but the factor U is exactly singular,
			and division by zero will occur if it is used to solve a system of equations. 
		*/
		throw new nec_exception("nec++:  LU Decomposition Failed: ",info);
	}
	
	
#ifdef NEC_MATRIX_CHECK
	cout << "atlas_solved = ";
	to_octave(a_in,n,ndim);

	cout << "atlas_ip = ";
	to_octave(ip,n);
#endif
} 

// extern "C" void zgetrs_(char* TRANS, int* N, int* NRHS, void *A, int* LDA, int *IPIV, void* B, int* LDB, int *INFOp);
// // SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
// int lapack_zgetrs(int N, int NRHS, void *A, int LDA, int *IPIV, void* B, int LDB)
// {
// 	int info;
// 	char transp = 'N';
// 	zgetrs_(&transp, &N, &NRHS, (void*)A, &LDA, IPIV, B, &LDB, &info);
// 	return info;
// }

/*! \brief
	Subroutine to solve the matrix equation lu*x=b where l is a unit
	lower triangular matrix and u is an upper triangular matrix both
	of which are stored in a.  the rhs vector b is input and the
	solution is returned through vector b.   (matrix transposed)
*/
void solve( int n, complex_array& a, int_array& ip,
    complex_array& b, int ndim )
{
	DEBUG_TRACE("solve(" << n << "," << ndim << ")");

	int info = clapack_zgetrs (CblasColMajor, CblasNoTrans, 
		n, 1, (void*) a.data(), ndim, ip.data(), b.data(), n);
	
	if (0 != info)
	{
		/*
			The factorization has been completed, but the factor U is exactly singular,
			and division by zero will occur if it is used to solve a system of equations. 
		*/
		throw new nec_exception("nec++: Solving Failed: ",info);
	}
	
}

#endif /*  LAPACK */


/*-----------------------------------------------------------------------*/

/**	factrs

	For symmetric structure, transforms submatricies to form
	matricies of the symmetric modes and calls routine to LU decompose
	matricies.
	
	If no symmetry [nrow = np], the routine is called to LU decompose the
	complete matrix.
*/
void factrs(nec_output_file& s_output,  int np, int nrow, complex_array& a, int_array& ip )
{
	DEBUG_TRACE("factrs(" << np << "," << nrow << ")");
	if (nrow == np) // no symmetry
	{
		lu_decompose(s_output,  np, a, ip, nrow );
		return;
	}
	
	int num_symmetric_modes = nrow / np;
	DEBUG_TRACE("\tnum_symmetric_modes = " << num_symmetric_modes);
	
	for (int mode = 0; mode < num_symmetric_modes; mode++ )
	{
		int mode_offset = mode * np;
		
		complex_array a_temp = a.segment(mode_offset, a.size()-mode_offset);
		int_array ip_temp = ip.segment(mode_offset, ip.size()-mode_offset);
		
		lu_decompose(s_output,  np, a_temp, ip_temp, nrow );
	}
}

/*-----------------------------------------------------------------------*/



/**
	Subroutine solves, for symmetric structures, handles the
	transformation of the right hand side vector and solution
	of the matrix eq.
	\param neq number of equations?
	\param nrh dimension of right hand  vector?
*/
void solves(complex_array& a, int_array& ip, complex_array& b, int neq,
	int nrh, int np, int n, int mp, int m, int nop, 
	complex_array& symmetry_array
)
{
	DEBUG_TRACE("solves(" << neq << "," << nrh << "," << np << "," << n << ")");
	DEBUG_TRACE("      ( nop=" << nop << ")");

	
	/* Allocate some scratch memory */
	complex_array scm;
	scm.resize(n + 2*m);
	
	int npeq= np+ 2*mp;
	nec_float fnop = nop;
	nec_float fnorm = 1.0/ fnop;
	int nrow= neq;
	
	if ( nop != 1)
	{
		for (int ic = 0; ic < nrh; ic++ )
		{
			int column_offset = ic*neq;
			if ( (n != 0) && (m != 0) )
			{
				for (int i = 0; i < neq; i++ )
					scm[i]= b[i+column_offset];
			
				int j= np-1;
			
				for (int k = 0; k < nop; k++ )
				{
					if ( k != 0 )
					{
						int ia= np-1;
						for (int i = 0; i < np; i++ )
						{
							ia++;
							j++;
							b[j+column_offset]= scm[ia];
						}
				
						if ( k == (nop-1) )
							continue;
					} /* if ( k != 0 ) */
				
					int mp2 = 2*mp;
					int ib= n-1;
					for (int i = 0; i < mp2; i++ )
					{
						ib++;
						j++;
						b[j+column_offset]= scm[ib];
					}
				} /* for( k = 0; k < nop; k++ ) */
			
			} /* if ( (n != 0) && (m != 0) ) */
		
			/* transform matrix eq. rhs vector according to symmetry modes */
			for (int i = 0; i < npeq; i++ )
			{
				for (int k = 0; k < nop; k++ )
				{
					int ia= i+ k* npeq;
					scm[k]= b[ia+column_offset];
				}
			
				nec_complex sum_normal(scm[0]);
				for (int k = 1; k < nop; k++ )
					sum_normal += scm[k];
			
				b[i+column_offset]= sum_normal * fnorm;
			
				for (int k = 1; k < nop; k++ )
				{
					int ia= i+ k* npeq;
					nec_complex sum(scm[0]);
				
					for (int j = 1; j < nop; j++ )
						sum += scm[j]* conj( symmetry_array[k+j*nop]);
				
					b[ia+column_offset]= sum* fnorm;
				}
			} /* for( i = 0; i < npeq; i++ ) */
		
		} /* for( ic = 0; ic < nrh; ic++ ) */
	
	} /* if ( nop != 1) */
	
	/* solve each mode equation */
	for (int kk = 0; kk < nop; kk++ )
	{
		int ia= kk* npeq;
	
		for (int ic = 0; ic < nrh; ic++ )
		{
			int column_offset = ic*neq;
			complex_array a_sub = a.segment(ia, a.size()-ia);
			complex_array b_sub = b.segment(ia+column_offset, b.size() - (ia+column_offset) );
			int_array ip_sub = ip.segment(ia, ip.size()-ia);
			solve( npeq, a_sub, ip_sub, b_sub, nrow );
		}
	
	} /* for( kk = 0; kk < nop; kk++ ) */
	
	if ( nop == 1)
	{
		return;
	}
	
	/* inverse transform the mode solutions */
	for (int ic = 0; ic < nrh; ic++ )
	{
		int column_offset = ic*neq;
		for (int i = 0; i < npeq; i++ )
		{
			for (int k = 0; k < nop; k++ )
			{
				int ia= i+ k* npeq;
				scm[k]= b[ia+column_offset];
			}
		
			nec_complex sum_normal(scm[0]);
			for (int k = 1; k < nop; k++ )
				sum_normal += scm[k];
		
			b[i+column_offset]= sum_normal;
			
			for (int k = 1; k < nop; k++ )
			{
				int ia= i+ k* npeq;
				
				nec_complex sum(scm[0]);
			
				for (int j = 1; j < nop; j++ )
					sum += scm[j]* symmetry_array[k+j*nop];
			
				b[ia+column_offset]= sum;
			}
		
		} /* for( i = 0; i < npeq; i++ ) */
	
		if ( (n == 0) || (m == 0) )
			continue;
	
		for (int i = 0; i < neq; i++ )
			scm[i]= b[i+column_offset];
	
		int j = np-1;
	
		for (int k = 0; k < nop; k++ )
		{
			if ( k != 0 )
			{
				int ia = np-1;
				for (int i = 0; i < np; i++ )
				{
					ia++;
					j++;
					b[ia+column_offset]= scm[j];
				}
			
				if ( k == nop)
					continue;
			
			} /* if ( k != 0 ) */
		
			int ib = n-1;
			int mp2 = 2* mp;
			for (int i = 0; i < mp2; i++ )
			{
				ib++;
				j++;
				b[ib+column_offset]= scm[j];
			}
		} /* for( k = 0; k < nop; k++ ) */
	
	} /* for( ic = 0; ic < nrh; ic++ ) */
}

/*-----------------------------------------------------------------------*/
/** \brief test for convergence in numerical integration */
void test(
	nec_float f1r, nec_float f2r, nec_float *tr,
 	nec_float f1i, nec_float f2i, nec_float *ti,
	nec_float dmin )
{
	static nec_float _min_val =  1.0e-37;
	
/*
{
  double den;

  den= fabs( f2r);
  *tr= fabs( f2i);

  if( den < *tr)
    den= *tr;
  if( den < dmin)
    den= dmin;

  if( den < 1.0e-37)
  {
    *tr=0.;
    *ti=0.;
    return;
  }

  *tr= fabs(( f1r- f2r)/ den);
  *ti= fabs(( f1i- f2i)/ den);
}

*/	
	nec_float den = fabs( f2r);
	nec_float temp_tr = fabs( f2i);
	
	if( den < temp_tr)
		den = temp_tr;
	if( den < dmin)
		den = dmin;
	
	if( den < _min_val)
	{
		*tr = 0.0;
		*ti = 0.0;
		return;
	}
	
	*tr= fabs((f1r - f2r)/ den);
	*ti= fabs((f1i - f2i)/ den); 
}

/**
	Simpler test for convergence in numerical integration.
	This tests only one number. It is a special case of the
	test() function above.
*/
nec_float test_simple( nec_float f1r, nec_float f2r, nec_float dmin )
{
	static nec_float _min_val =  1.0e-37;
	
	nec_float den = fabs(f2r);
	
	if( den < dmin)
		den = dmin;
	if (den < _min_val)
	{
		return 0.0;
	}
	
	return fabs((f1r - f2r) / den);	
}
