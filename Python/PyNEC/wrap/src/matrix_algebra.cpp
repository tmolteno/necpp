/*
	Copyright (C) 2004-2005  Timothy C.A. Molteno
	
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
#include "math_util.h"
#include <iostream>
using namespace std;

#include "nec_exception.h"
#include "matrix_algebra.h"
#include "nec_output.h"

#ifdef NEC_ERROR_CHECK
//#define NEC_MATRIX_CHECK
#endif

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

/*
	Using Octave to debug this stuff.
	
	[l,u,p] = lu(a);
	
	orig = p*l*u;
	
	diff = orig .- a;
	
	sum(sum(diff))
	
	
	
a = [-5.47894 + 2965.96I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -2.21543 + -0.962361I, -0.795364 + -0.113862I; -2.49805 + -1034.62I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -1.01 + -0.527357I; -1.67535 + -27.5893I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -1.21709 + -1.25704I; -1.56046 + -6.58465I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -1.40441 + -2.71197I; -1.40441 + -2.71197I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -1.56046 + -6.58465I; -1.21709 + -1.25704I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -1.67535 + -27.5892I; -1.01 + -0.527357I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -2.49805 + -1034.62I; -0.795364 + -0.113862I, -2.21543 + -0.962361I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -5.47894 + 2965.96I];

solved = [-5.47894 + 2965.96I, -0.34883 + 0.00148662I, -0.00930089 + 0.00058204I, -0.00221909 + 0.000530222I, -0.000913486 + 0.000475198I, -0.000423063 + 0.000411133I, -0.000177173 + 0.000340856I, -3.7894e-05 + 0.000268234I; -5.16779 + -1625.83I, -15.1883 + 2508.51I, -0.654113 + 0.00641696I, -0.0209388 + 0.00199762I, -0.00563214 + 0.00175056I, -0.00239586 + 0.0015366I, -0.00108368 + 0.00130452I, -0.000401231 + 0.00105952I; -3.8195 + -48.9503I, -6.57291 + -1642.9I, -25.8745 + 2000.59I, -0.829723 + 0.0150408I, -0.0290603 + 0.00375449I, -0.00826549 + 0.00315121I, -0.00352744 + 0.00270246I, -0.00152297 + 0.00224776I; -3.52843 + -12.6722I, -5.06917 + -53.3655I, -8.86623 + -1660.82I, -43.5325 + 1696.62I, -0.986053 + 0.0322509I, -0.0368356 + 0.00638373I, -0.0108087 + 0.00506807I, -0.00453205 + 0.00421639I; -3.14757 + -5.34769I, -4.63435 + -14.5329I, -6.97654 + -58.4747I, -11.9718 + -1674.55I, -16.5775 + -1687.93I, -0.842565 + -0.0540462I, 0.0398738 + -0.00709961I, 0.0120791 + -0.00628218I; -2.70084 + -2.45195I, -4.09335 + -6.19899I, -6.27227 + -16.722I, -9.38065 + -62.8562I, -11.8415 + 3073.26I, -27.0157 + -1749.13I, -0.513164 + -0.118154I, 0.049591 + -0.0125818I; -2.21543 + -0.962361I, -3.47508 + -2.78436I, -5.4597 + -7.15434I, -8.24983 + -18.5785I, -5.67188 + -1626.5I, 0.575028 + 3140.29I, -34.565 + -1761.99I, -0.095149 + -0.172344I; -0.795364 + -0.113862I, -1.28761 + -0.565893I, -2.07043 + -1.61953I, -3.17657 + -4.0337I, -1.84467 + -27.7235I, -2.29644 + -1033.56I, 7.796 + 3017.51I, -403.323 + -276.157I];
	
[l,u,p] = lu(a);
ans = (l .- eye(size(a))) .+ u;
solved .- ans
	
	
	
	

atlas_a = [-5.47894 + 2965.96I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -2.21543 + -0.962361I, -0.795364 + -0.113862I; -2.49805 + -1034.62I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -1.01 + -0.527357I; -1.67535 + -27.5893I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -1.21709 + -1.25704I; -1.56046 + -6.58465I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -1.40441 + -2.71197I; -1.40441 + -2.71197I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -1.56046 + -6.58465I; -1.21709 + -1.25704I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -1.67535 + -27.5892I; -1.01 + -0.527357I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -2.49805 + -1034.62I; -0.795364 + -0.113862I, -2.21543 + -0.962361I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -5.47894 + 2965.96I];

atlas_solved = [-5.47894 + 2965.96I, -0.548157 + 0.00275496I, -0.0165016 + 0.00131826I, -0.00427032 + 0.00119753I, -0.000825012 + 0.000912136I, -0.000323087 + 0.00074755I, -3.7894e-05 + 0.000268234I, -0.00180105 + 0.00106456I; -2.49805 + -1034.62I, -15.1883 + 2508.51I, -0.65489 + 0.00658541I, -0.0212608 + 0.00214951I, -0.00246121 + 0.00164669I, -0.00110154 + 0.00139199I, -0.000222473 + 0.000514644I, -0.00578205 + 0.00188246I; -1.67535 + -27.5893I, -6.16215 + -1640.95I, -25.8745 + 2000.59I, -0.82997 + 0.0151662I, -0.00831659 + 0.00324278I, -0.00354023 + 0.00277484I, -0.00079601 + 0.0010452I, -0.0291788 + 0.00386463I; -1.56046 + -6.58465I, -4.69302 + -52.5555I, -8.62173 + -1660.32I, -43.5325 + 1696.62I, -0.0368818 + 0.00647535I, -0.0108184 + 0.00514009I, -0.00232792 + 0.00193202I, -0.98616 + 0.0323595I; -1.40441 + -2.71197I, -4.30575 + -14.1549I, -6.75926 + -58.2348I, -11.7922 + -1674.36I, -16.7274 + -1688I, 0.0398818 + -0.00716753I, 0.00624829 + -0.0028398I, -0.842524 + -0.0541183I; -1.21709 + -1.25704I, -3.81819 + -6.03339I, -6.0904 + -16.6174I, -9.22723 + -62.7739I, -11.8415 + 3073.26I, -27.2273 + -1749.16I, 0.026934 + -0.00558256I, -0.513109 + -0.118254I; -1.01 + -0.527357I, -3.25593 + -2.73825I, -5.31524 + -7.12687I, -8.12807 + -18.5589I, -5.66889 + -1626.5I, 0.685432 + 3140.3I, -15.4618 + -1109.08I, -0.152662 + -0.272982I; -0.795364 + -0.113862I, -2.65173 + -1.02258I, -4.45744 + -3.105I, -6.95632 + -7.87272I, -4.18315 + -49.2046I, -4.79308 + -1623.94I, 3.84241 + 3009.97I, -643.167 + -435.236I];
[l,u,p] = lu(atlas_a');
atlas_ans = (l .- eye(size(atlas_a))) .+ u;
atlas_solved' .- ans
*/
	
#ifndef ATLAS
/*
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using the Gauss-Doolittle algorithm 
	presented on pages 411-416 of A. Ralston -- a first course in 
	numerical analysis.
	
	Comments below refer to comments in Ralstons text.
	
	(matrix transposed.)
*/
void lu_decompose(nec_output_file& s_output, int n, complex_array& a, int_array& ip, int ndim)
{
	DEBUG_TRACE("lu_decompose(" << n << "," << ndim << ")");
	
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

#else /*  ATLAS */
/* 
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using ATLAS (matrix transposed.)

    zgetrf computes an LU  factorization  of  a  general  M-by-N
     matrix A using partial pivoting with row interchanges.

     The factorization has the form
        A = P * L * U
     where P is a permutation matrix, L is lower triangular  with
     unit  diagonal  elements (lower trapezoidal if m > n), and U
     is upper triangular (upper trapezoidal if m < n).

     This is the right-looking Level 3 BLAS version of the  algo-
     rithm.
	
	int clapack_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   void *A, const int lda, int *ipiv);
		   
		CblasRowMajor 	assumes we are storing A in C style
		CblasColMajor	assume we are storing A in FORTRAN style
		   
    A (input/output)
		On entry, the M-by-N matrix to  be  factored.   On
		exit, the factors L and U from the factorization A
		= P*L*U; the unit diagonal elements of L  are  not
		stored.
			   	
	ipiv	(output) INTEGER array, dimension (min(M,N))   
		The pivot indices; for 1 <= i <= min(M,N), row i of the   
		matrix was interchanged with row IPIV(i).
	
	input parameters
		ndim		actual full matrix dimension
		n		dimension of submatrix to be factored
		a_in[ndim,ndim]	The full input matrix
		ip[ndim]	The pivot points
		

 */
extern "C"
{
#include <atlas_enum.h>
#include <clapack.h>
}

void lapack_test();
void lapack_test()
{
/*	test_a = [1 + 0I, 2 + 0I; 3 + 0I, 4 + 0I];
	test_ans = [2 + 0I, 0.5 + 0I; 4 + 0I, 1 + 0I];
	piv = [1,1];
	
	We need to use the transpose operator here since we are storing the matrix as column major
	
	[l,u,p] = lu(test_a');
	ans = (l .- eye(size(test_a))) .+ u;
	
	We also compare to the transpose as the result matrix is in column major order as well
	test_ans' .- ans
*/
	complex_array A(4);
	int_array piv(2);
	
	{
		A[0] = 3;
		A[1] = 1;
		A[2] = 4;
		A[3] = 2;
		
		cout << "test_a = ";
		to_octave(A,2,2);
		
		// Now call the LAPACK LU-Decomposition
		int info = clapack_zgetrf (CblasColMajor, 2, 2, (void*) A.get_ptr(), 2, piv.get_ptr());
		
		std::cout << "CblasColMajor: " << info << " : " << endl;
		cout << "test_ans = ";
		to_octave(A,2,2);
		std::cout << "piv = ";
		to_octave(piv,2);
	}
}

/*! Use lapack to perform LU decomposition
*/
void lu_decompose(nec_output_file& s_output,  int n, complex_array& a_in, int_array& ip, int ndim)
{
	DEBUG_TRACE("factor_lapack(" << n << "," << ndim << ")");
	ASSERT(n <= ndim);

#ifdef NEC_MATRIX_CHECK
	lapack_test();
	
	cout << "atlas_a = ";
	to_octave(a_in,n,ndim);
#endif

	// copy the input matrix a_in into a temporary array.
	// transposing as we go... Should use cblas_zgemm...
	complex_array A(n,n);
	int_array piv(n);
	
	/* Un-transpose the matrix for Gauss elimination */
	for (int row = 0; row < n; row++ )
	{
		int row_start = row * ndim;
		for (int col = 0; col < n; col++ )
		{
			A.set_col_major(n,col,row,a_in[col + row_start]);
		}
	}
	
	int lead_dim = std::max(1, n);
		
	// Now call the LAPACK LU-Decomposition
	int info = clapack_zgetrf (CblasColMajor, n, n, 
		(void*) A.get_ptr(), lead_dim, piv.get_ptr());
	
	if (0 != info)
	{
		/*
			The factorization has been completed, but the factor U is exactly singular,
			and division by zero will occur if it is used to solve a system of equations. 
		*/
		nec_exception* nex = new nec_exception("nec++: LU Decomposition Failed: ");
		nex->append(info);
		throw nex;
	}
	
	
	/*
	IPIV	(output) INTEGER array, dimension (min(M,N))   
		The pivot indices; for 1 <= i <= min(M,N), row i of the   
		matrix was interchanged with row IPIV(i).
	*/  
	for (int j = 0; j < n; j++ )
	{
		ip[j] = piv[j] + 1;
	}
		
	// copy the output back into the a_in array.
	for (int row = 0; row < n; row++ )
	{
		int row_start = row*ndim;
		
		for (int col = 0; col < n; col++ )
		{
			a_in[row_start + col] = A.get_col_major(n,col,row);
		}
	}
	
#ifdef NEC_MATRIX_CHECK
	cout << "atlas_solved = ";
	to_octave(a_in,n,ndim);

	cout << "atlas_ip = ";
	to_octave(ip,n);
#endif
} 
#endif /*  ATLAS */


/*-----------------------------------------------------------------------*/

/*	factrs

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
		
		complex_array a_temp = a.sub_array(mode_offset);
		int_array ip_temp = ip.sub_array(mode_offset);
		
		lu_decompose(s_output,  np, a_temp, ip_temp, nrow );
	}
}

/*-----------------------------------------------------------------------*/

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



/*
	Subroutine solves, for symmetric structures, handles the
	transformation of the right hand side vector and solution
	of the matrix eq.
*/
void solves(complex_array& a, int_array& ip, complex_array& b, int neq,
	int nrh, int np, int n, int mp, int m, int nop, 
	complex_array& symmetry_array
)
{
	DEBUG_TRACE("solves(" << neq << "," << nrh << "," << np << "," << n << ")");
	int  ic, i, kk, ia, ib, j, k;

	nec_complex  sum;
	
	/* Allocate to scratch memory */
	complex_array scm;
	scm.resize(n + 2*m);
	
	int npeq= np+ 2*mp;
	nec_float fnop = nop;
	nec_float fnorm = 1.0/ fnop;
	int nrow= neq;
	
	if ( nop != 1)
	{
		for( ic = 0; ic < nrh; ic++ )
		{
			if ( (n != 0) && (m != 0) )
			{
				for( i = 0; i < neq; i++ )
					scm[i]= b[i+ic*neq];
			
				kk=2* mp;
				ia= np-1;
				ib= n-1;
				j= np-1;
			
				for( k = 0; k < nop; k++ )
				{
					if ( k != 0 )
					{
						for( i = 0; i < np; i++ )
						{
							ia++;
							j++;
							b[j+ic*neq]= scm[ia];
						}
				
						if ( k == (nop-1) )
							continue;
					} /* if ( k != 0 ) */
				
					for( i = 0; i < kk; i++ )
					{
						ib++;
						j++;
						b[j+ic*neq]= scm[ib];
					}
				} /* for( k = 0; k < nop; k++ ) */
			
			} /* if ( (n != 0) && (m != 0) ) */
		
			/* transform matrix eq. rhs vector according to symmetry modes */
			for( i = 0; i < npeq; i++ )
			{
				for( k = 0; k < nop; k++ )
				{
					ia= i+ k* npeq;
					scm[k]= b[ia+ic*neq];
				}
			
				sum= scm[0];
				for( k = 1; k < nop; k++ )
					sum += scm[k];
			
				b[i+ic*neq]= sum* fnorm;
			
				for( k = 1; k < nop; k++ )
				{
					ia= i+ k* npeq;
					sum= scm[0];
				
					for( j = 1; j < nop; j++ )
						sum += scm[j]* conj( symmetry_array[k+j*nop]);
				
					b[ia+ic*neq]= sum* fnorm;
				}
			} /* for( i = 0; i < npeq; i++ ) */
		
		} /* for( ic = 0; ic < nrh; ic++ ) */
	
	} /* if ( nop != 1) */
	
	/* solve each mode equation */
	for( kk = 0; kk < nop; kk++ )
	{
		ia= kk* npeq;
		ib= ia;
	
		for( ic = 0; ic < nrh; ic++ )
		{
			complex_array a_sub = a.sub_array(ib);
			complex_array b_sub = b.sub_array(ia+ic*neq);
			int_array ip_sub = ip.sub_array(ia);
			solve( npeq, a_sub, ip_sub, b_sub, nrow );
		}
	
	} /* for( kk = 0; kk < nop; kk++ ) */
	
	if ( nop == 1)
	{
		return;
	}
	
	/* inverse transform the mode solutions */
	for( ic = 0; ic < nrh; ic++ )
	{
		for( i = 0; i < npeq; i++ )
		{
			for( k = 0; k < nop; k++ )
			{
				ia= i+ k* npeq;
				scm[k]= b[ia+ic*neq];
			}
		
			sum= scm[0];
			for( k = 1; k < nop; k++ )
				sum += scm[k];
		
			b[i+ic*neq]= sum;
			for( k = 1; k < nop; k++ )
			{
				ia= i+ k* npeq;
				sum= scm[0];
			
				for( j = 1; j < nop; j++ )
					sum += scm[j]* symmetry_array[k+j*nop];
			
				b[ia+ic*neq]= sum;
			}
		
		} /* for( i = 0; i < npeq; i++ ) */
	
		if ( (n == 0) || (m == 0) )
			continue;
	
		for( i = 0; i < neq; i++ )
			scm[i]= b[i+ic*neq];
	
		kk=2* mp;
		ia= np-1;
		ib= n-1;
		j= np-1;
	
		for( k = 0; k < nop; k++ )
		{
			if ( k != 0 )
			{
				for( i = 0; i < np; i++ )
				{
					ia++;
					j++;
					b[ia+ic*neq]= scm[j];
				}
			
				if ( k == nop)
					continue;
			
			} /* if ( k != 0 ) */
		
			for( i = 0; i < kk; i++ )
			{
				ib++;
				j++;
				b[ib+ic*neq]= scm[j];
			}
		} /* for( k = 0; k < nop; k++ ) */
	
	} /* for( ic = 0; ic < nrh; ic++ ) */
}

/*-----------------------------------------------------------------------*/
/* test for convergence in numerical integration */
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

/*
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
