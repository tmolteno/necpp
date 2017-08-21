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



	Routines for checking the LU decomposition in NEC. Goal is to get it going
	using LAPACK.

a = {	{3 + 0I, 3 + 0I, 2 + 0I, -2 + 0I},	{1 + 0I, 1 + 0I, 13 + 0I, 3 + 0I},	{-4 + 0I, 0 + 0I, -1 + 0I, -1 + 0I},	{ 2 + 0I, 2 + 0I, 0 + 0I, 4 + 0I}}


FORTRAN output

    3.00000000000000
    1.00000000000000
   0.666666666666667
   0.666666666666667


    1.00000000000000
    12.3333333333333
    0.00000000000000
   0.297297297297297


    0.00000000000000
    1.00000000000000
    4.00000000000000
   0.175675675675676


    2.00000000000000
    1.33333333333333
    0.00000000000000
    5.72972972972973
            2
            3
            3
            4


test_a = [3 + 0I, 1 + 0I, -4 + 0I, 2 + 0I; 3 + 0I, 1 + 0I, 0 + 0I, 2 + 0I; 2 + 0I, 13 + 0I, -1 + 0I, 0 + 0I; -2 + 0I, 3 + 0I, -1 + 0I, 4 + 0I];
test_ans = [-4 + 0I, -0.75 + -0I, -0.25 + -0I, -0.5 + -0I; 0 + 0I, 3 + 0I, 0.333333 + 0I, 0.666667 + 0I; -1 + 0I, 1.25 + 0I, 12.3333 + -0I, -0.108108 + -0I; -1 + 0I, -2.75 + 0I, 3.66667 + 0I, 5.72973 + -0I];
piv = [2, 2, 2, 3];
lu_dcpse_old = [3 + 0I, 1 + 0I, 0.666667 + 0I, -0.666667 + 0I; 1 + 0I, 12.3333 + 0I, 0 + 0I, 0.297297 + 0I; 0 + 0I, -1 + 0I, -4 + 0I, 0.175676 + 0I; 2 + 0I, -1.33333 + 0I, 0 + 0I, 5.72973 + 0I];
piv = [2, 3, 3, 4];
lu_dcpse_bke = [3 + 0I, 1 + 0I, 0.666667 + 0I, -0.666667 + 0I; 1 + 0I, 12.3333 + 0I, 0 + 0I, 0.297297 + 0I; 0 + 0I, -1 + 0I, -4 + 0I, 0.175676 + 0I; 2 + 0I, -1.33333 + 0I, 0 + 0I, 5.72973 + 0I];
piv = [2, 3, 3, 4];
lu_decompose = [3 + 0I, 1 + 0I, -4 + 0I, 2 + 0I; 0.666667 + 0I, 12.3333 + -0I, 1.66667 + 0I, -1.33333 + 0I; -0.666667 + 0I, 0.297297 + -0I, -4.16216 + -0I, 5.72973 + 0I; 1 + 0I, -0 + -0I, -0.961039 + 0I, 5.50649 + -0I];
piv = [1, 3, 4, 4];



function [l u p]=aprt(a)
	l = eye(size(a)) .+ tril(a,-1);
	u = triu(a);
	p = 1:length(a);
endfunction

test_a

[l u p] = lu(test_a)
l * u / test_a > 0.001
[l,u,pp] = aprt(lu_dcpse_bke)
l * u / test_a > 0.001
[l,u,pp] = aprt(lu_decompose)
l * u / test_a > 0.001
[l,u,pp] = aprt(lu_dcpse_old)
l * u / test_a > 0.001


*/
#if NEC_ERROR_CHECK
#include "math_util.h"

#include <iostream>
using namespace std;


/* Various Routines for sending matrices out to the console in a form that can
be pasted into OCTAVE (http://www.octave.org)
*/
void to_octave(nec_complex& x);
void to_octave(nec_complex& x)
{
	cout << real(x);// << " + " << imag(x) << "I";
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
	
	
	
a = [-5.47894 + 2965.96I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -2.21543 + -0.962361I, -0.795364 + -0.113862I; -2.49805 + -1034.62I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -1.01 + -0.527357I; -1.67535 + -27.5893I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -1.21709 + -1.25704I; -1.56046 + -6.58465I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -1.40441 + -2.71197I; -1.40441 + -2.71197I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -1.56046 + -6.58465I; -1.21709 + -1.25704I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, 
-5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -1.67535 + -27.5892I; -1.01 + -0.527357I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -2.49805 + -1034.62I; -0.795364 + -0.113862I, -2.21543 + -0.962361I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -5.47894 + 2965.96I];

solved = [-5.47894 + 2965.96I, -0.34883 + 0.00148662I, -0.00930089 + 0.00058204I, -0.00221909 + 0.000530222I, -0.000913486 + 0.000475198I, -0.000423063 + 0.000411133I, -0.000177173 + 0.000340856I, -3.7894e-05 + 0.000268234I; -5.16779 + -1625.83I, -15.1883 + 2508.51I, -0.654113 + 0.00641696I, -0.0209388 + 0.00199762I, -0.00563214 + 0.00175056I, -0.00239586 + 0.0015366I, -0.00108368 + 0.00130452I, -0.000401231 + 0.00105952I; -3.8195 + -48.9503I, -6.57291 + -1642.9I, -25.8745 + 2000.59I, -0.829723 + 0.0150408I, -0.0290603 + 0.00375449I, -0.00826549 + 0.00315121I, -0.00352744 + 0.00270246I, -0.00152297 + 0.00224776I; -3.52843 + -12.6722I, -5.06917 + -53.3655I, -8.86623 + -1660.82I, -43.5325 + 1696.62I, -0.986053 + 0.0322509I, -0.0368356 + 0.00638373I, -0.0108087 + 0.00506807I, -0.00453205 + 0.00421639I; -3.14757 + -5.34769I, -4.63435 + -14.5329I, -6.97654 + -58.4747I, -11.9718 + -1674.55I, -16.5775 + -1687.93I, -0.842565 + -0.0540462I, 
0.0398738 + -0.00709961I, 0.0120791 + -0.00628218I; -2.70084 + -2.45195I, -4.09335 + -6.19899I, -6.27227 + -16.722I, -9.38065 + -62.8562I, -11.8415 + 3073.26I, -27.0157 + -1749.13I, -0.513164 + -0.118154I, 0.049591 + -0.0125818I; -2.21543 + -0.962361I, -3.47508 + -2.78436I, -5.4597 + -7.15434I, -8.24983 + -18.5785I, -5.67188 + -1626.5I, 0.575028 + 3140.29I, -34.565 + -1761.99I, -0.095149 + -0.172344I; -0.795364 + -0.113862I, -1.28761 + -0.565893I, -2.07043 + -1.61953I, -3.17657 + -4.0337I, -1.84467 + -27.7235I, -2.29644 + -1033.56I, 7.796 + 3017.51I, -403.323 + -276.157I];
	
[l,u,p] = lu(a);
ans = (l .- eye(size(a))) .+ u;
solved .- ans
	
	
	

atlas_a = [-5.47894 + 2965.96I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -2.21543 + -0.962361I, -0.795364 + -0.113862I; -2.49805 + -1034.62I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -2.70084 + -2.45195I, -1.01 + -0.527357I; -1.67535 + -27.5893I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -3.14757 + -5.34769I, -1.21709 + -1.25704I; -1.56046 + -6.58465I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -3.52843 + -12.6722I, -1.40441 + -2.71197I; -1.40441 + -2.71197I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -3.8195 + -48.9503I, -1.56046 + -6.58465I; -1.21709 + -1.25704I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, 
-5.16779 + -1625.83I, -10.9686 + 3075.64I, -5.16779 + -1625.83I, -1.67535 + -27.5892I; -1.01 + -0.527357I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -10.9686 + 3075.64I, -2.49805 + -1034.62I; -0.795364 + -0.113862I, -2.21543 + -0.962361I, -2.70084 + -2.45195I, -3.14757 + -5.34769I, -3.52843 + -12.6722I, -3.8195 + -48.9504I, -5.16779 + -1625.83I, -5.47894 + 2965.96I];

atlas_solved = [-5.47894 + 2965.96I, -0.548157 + 0.00275496I, -0.0165016 + 0.00131826I, -0.00427032 + 0.00119753I, -0.000825012 + 0.000912136I, -0.000323087 + 0.00074755I, -3.7894e-05 + 0.000268234I, -0.00180105 + 0.00106456I; -2.49805 + -1034.62I, -15.1883 + 2508.51I, -0.65489 + 0.00658541I, -0.0212608 + 0.00214951I, -0.00246121 + 0.00164669I, -0.00110154 + 0.00139199I, -0.000222473 + 0.000514644I, -0.00578205 + 0.00188246I; -1.67535 + -27.5893I, -6.16215 + -1640.95I, -25.8745 + 2000.59I, -0.82997 + 0.0151662I, -0.00831659 + 0.00324278I, -0.00354023 + 0.00277484I, -0.00079601 + 0.0010452I, -0.0291788 + 0.00386463I; -1.56046 + -6.58465I, -4.69302 + -52.5555I, -8.62173 + -1660.32I, -43.5325 + 1696.62I, -0.0368818 + 0.00647535I, -0.0108184 + 0.00514009I, -0.00232792 + 0.00193202I, -0.98616 + 0.0323595I; -1.40441 + -2.71197I, -4.30575 + -14.1549I, -6.75926 + -58.2348I, -11.7922 + -1674.36I, -16.7274 + -1688I, 0.0398818 + -0.00716753I, 
0.00624829 + -0.0028398I, -0.842524 + -0.0541183I; -1.21709 + -1.25704I, -3.81819 + -6.03339I, -6.0904 + -16.6174I, -9.22723 + -62.7739I, -11.8415 + 3073.26I, -27.2273 + -1749.16I, 0.026934 + -0.00558256I, -0.513109 + -0.118254I; -1.01 + -0.527357I, -3.25593 + -2.73825I, -5.31524 + -7.12687I, -8.12807 + -18.5589I, -5.66889 + -1626.5I, 0.685432 + 3140.3I, -15.4618 + -1109.08I, -0.152662 + -0.272982I; -0.795364 + -0.113862I, -2.65173 + -1.02258I, -4.45744 + -3.105I, -6.95632 + -7.87272I, -4.18315 + -49.2046I, -4.79308 + -1623.94I, 3.84241 + 3009.97I, -643.167 + -435.236I];
[l,u,p] = lu(atlas_a');
atlas_ans = (l .- eye(size(atlas_a))) .+ u;
atlas_solved' .- ans
*/
	
/*! \brief Solve The system of equations using Gaussian Elimination.
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using the Gauss-Doolittle algorithm 
	presented on pages 411-416 of A. Ralston -- a first course in 
	numerical analysis.
	
	Comments below refer to comments in Ralstons text.
	
	(matrix transposed.)
*/
void lu_decompose_old(int n, complex_array& a, int_array& ip, int ndim);
void lu_decompose_old(int n, complex_array& a, int_array& ip, int ndim)
{	
	// Debug output to try and figure out the LAPACK stuff
	/* Allocate scratch memory */
	complex_array scm;
	scm.resize(n);
	
	/* Un-transpose the matrix for Gauss elimination */
/*	for (int i = 1; i < n; i++ )
	{
		for (int j = 0; j < i; j++ )
			std::swap(a.get(i,j),a.get(j,i));
	}*/
	
	bool iflg=false;
	for (int r = 0; r < n; r++ )
	{
		int r_offset = r*ndim;
		
		/* step 1 */
		for (int k = 0; k < n; k++ )
			scm[k]= a[k+r_offset];
		
		/* steps 2 and 3 */
		for (int j = 0; j < r; j++ )
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
			{
				a[i+r_offset]= scm[i]* arr;
			}
		}
		
		if ( true == iflg )
		{
			cout << "\n  PIVOT(" << r << ")= " << dmax;
			iflg=false;
		}	
	} /* for( r=0; r < n; r++ ) */
	
/*
	cout << "solved = ";
	to_octave(a,n,ndim);

	cout << "ip = ";
	to_octave(ip,n);
	*/
}



void lu_decompose_burke(int n, complex_array& a, int_array& ip, int ndim);
void lu_decompose_burke(int n, complex_array& a, int_array& ip, int ndim)
{	
/*C
C Un-transpose the matrix for Gauss elimination
C
      DO 12 I=2,N
         DO 11 J=1,I-1
            ARJ=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=ARJ
11 CONTINUE
12 CONTINUE
*/
// 	for (int i = 1; i < n; i++ )
// 	{
// 		for (int j = 0; j < i; j++ )
// 			std::swap(a.get(i,j),a.get(j,i));
// 	}

 /*
     IFLG=0
      DO 9 R=1,N
*/
	complex_array d(n);
	for (int r = 0; r < n; r++ )
	{
		bool iflg=false;
/*C
C STEP 1
C
      DO 1 K=1,N
      D(K)=A(K,R)
1 CONTINUE*/
		for (int k = 0; k < n; k++ )
			d[k]= a.get(k,r);

/*
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
*/
		int rm1 = r - 1;
		if (rm1 >= 0)
		{
			for (int j=0; j < r; j++)
			{
				int pj = ip[j];
				nec_complex arj = d[pj];
				a.set(j,r,arj);
				d[pj] = d[j];
				int jp1 = j + 1;
				for (int i=jp1;i<n;i++)
					d[i] -= a.get(i,j)*arj;
			}
		}
/*
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
*/
		nec_float dmax = norm(d[r]);;
		ip[r] = r;
		int rp1 = r + 1;
		if (rp1 < n)
		{
			for (int i=rp1; i<n; i++)
			{
				nec_float elmag = norm(d[i]);
				if (elmag >= dmax)
				{
					dmax = elmag;
					ip[r] = i;
				}
			}
		}
		if (dmax < 1e-10)
			iflg = true;
		int pr = ip[r];
		a.set(r,r,d[pr]);
		d[pr] = d[r];

/*
C
C STEP 5
C
      IF (RP1.GT.N) GO TO 8
      ARJ=1./A(R,R)
      DO 7 I=RP1,N
      A(I,R)=D(I)*ARJ
7 CONTINUE
8 CONTINUE
*/
		if (rp1 < n)
		{
			nec_complex arj = cplx_10() / a.get(r,r);
			for (int i=rp1; i<n; i++)
			{
				a.set(i,r,d[i] * arj);
			}
		}
/*
      IF (IFLG.EQ.0) GO TO 9
      WRITE(3,10) R,DMAX
      IFLG=0
*/
		if (iflg == true)
		{
			cout << "FACTR: PIVOT(" << r << ")=" << dmax;
			iflg = false;
		}
	} 

	// increment ip array
	for (int i=0;i<n;i++)
		ip[i] += 1;
/*
	cout << "solved = ";
	to_octave(a,n,ndim);

	cout << "ip = ";
	to_octave(ip,n);
	*/
}

extern "C"{
#include <atlas_enum.h>
#include <clapack.h>
}


/*! Use lapack to perform LU decomposition
*/
void lu_decompose(int n, complex_array& a_in, int_array& ip, int ndim);
void lu_decompose(int n, complex_array& a_in, int_array& ip, int ndim)
{
/*	
	cout << "atlas_a = ";
	to_octave(a_in,n,ndim);

*/
	// copy the input matrix a_in into a temporary array (transposing as we go)
//	complex_array A(n,n);
	int_array piv(n);
	
/*	for (int row = 0; row < n; row++ )
	{
		int col_index = row * ndim;
		for (int col = 0; col < n; col++ )
		{
			A.set(row,col,a_in[col_index++]);
		}
	}*/
	
	int info = clapack_zgetrf (CblasColMajor, n, n, 
		(void*) a_in.get_ptr(), ndim, piv.get_ptr());
	
	if (0 != info)
	{
		/*
			The factorization has been completed, but the factor U is exactly singular,
			and division by zero will occur if it is used to solve a system of equations. 
		*/
		cout << "nec++: LU Decomposition Failed: " << info;
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
/*	for (int row = 0; row < n; row++ )
	{
		int col_index = row*ndim;
		
		for (int col = 0; col < n; col++ )
		{
			a_in[col_index++] = A.get(row,col);
		}
	}*/
/*	
	cout << "atlas_solved = ";
	to_octave(a_in,n,ndim);

	cout << "atlas_ip = ";
	to_octave(ip,n);*/
} 

void matrix_setup(complex_array& A);
void matrix_setup(complex_array& A)
{	int n = 4;
	A.get(0,0) = 3.0;
	A.get(0,1) =1.0;
	A.get(0,2) = -4.0;
	A.get(0,3) =2.0;
	A.get(1,0) =3.0;
	A.get(1,1) =1.0;
	A.get(1,2) =0.0;
	A.get(1,3) =2.0;
	A.get(2,0) =2.0;
	A.get(2,1) =13.0;
	A.get(2,2) = -1.0;
	A.get(2,3) =0.0;
	A.get(3,0) = -2.0;
	A.get(3,1) =3.0;
	A.get(3,2) = -1.0;
	A.get(3,3) =4.0;
	for (int i = 1; i < n; i++ )
	{
		for (int j = 0; j < i; j++ )
			std::swap(A.get(i,j),A.get(j,i));
	}
}

int main()
{
/*	test_a = [1 + 0I, 2 + 0I; 3 + 0I, 4 + 0I];
	test_ans = [2 + 0I, 0.5 + 0I; 4 + 0I, 1 + 0I];
	piv = [1,1];
	
	We need to use the transpose operator here since we are storing the matrix as column major
	
	[l,u,p] = lu(test_a');
	ret = (l .- eye(size(test_a))) .+ u;
	
	We also compare to the transpose as the result matrix is in column major order as well
	test_ans' .- ret
*/
	int N = 4;
	{
		complex_array A(N,N);
		int_array piv(N);
		matrix_setup(A);
		
		cout << "test_a = ";
		to_octave(A,N,N);
		
		// Now call the LAPACK LU-Decomposition
		int info = clapack_zgetrf (CblasColMajor, N, N, (void*) A.get_ptr(), N, piv.get_ptr());

//		std::cout << "CblasColMajor: " << info << " : " << endl;
		cout << "test_ans = ";
		to_octave(A,N,N);
		std::cout << "piv = ";
		to_octave(piv,N);

//		info = clapack_zgetrs (CblasColMajor, 2, 2, (void*) A.get_ptr(), 2, piv.get_ptr());
//		cout << "test_ans = ";
//		to_octave(A,2,2);

	}

	{
		complex_array A(N,N);
		int_array piv(N);
		matrix_setup(A);
		
		lu_decompose_old(N, A, piv, N);
		cout << "lu_dcpse_old = ";
		to_octave(A,N,N);
		std::cout << "piv = ";
		to_octave(piv,N);
	}
	{
		complex_array A(N,N);
		int_array piv(N);
		matrix_setup(A);
		
		lu_decompose_burke(N, A, piv, N);
		cout << "lu_dcpse_bke = ";
		to_octave(A,N,N);
		std::cout << "piv = ";
		to_octave(piv,N);
	}
	{
		complex_array A(N,N);
		int_array piv(N);
		matrix_setup(A);
		
		lu_decompose(N, A, piv, N);
		cout << "lu_decompose = ";
		to_octave(A,N,N);
		std::cout << "piv = ";
		to_octave(piv,N);
	}


	return 0;
}

#endif /* NEC_ERROR_CHECK */
