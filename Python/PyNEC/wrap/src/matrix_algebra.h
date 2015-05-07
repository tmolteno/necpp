#ifndef __matrix_algebra__
#define __matrix_algebra__

/*
  Copyright (C) 2004, 2015  Timothy C.A. Molteno
  
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
#include "nec_output.h"

/** \brief LU that uses either lapack or built in
 * */
void lu_decompose(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim);
void solve( int n, complex_array& a, int_array& ip, complex_array& b, int64_t ndim );

/** \brief LU that uses built in functions (no library dependency)
 * */
void lu_decompose_ge(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim);
void solve_ge( int n, complex_array& a, int_array& ip, complex_array& b, int64_t ndim );

#if LAPACK
/** \brief LU that uses LAPACK. These are not built unless a LAPACK is chosen at
 * compile time.
 * */
void lu_decompose_lapack(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim);
void solve_lapack( int n, complex_array& a, int_array& ip, complex_array& b, int64_t ndim );
#endif


void factrs(nec_output_file& s_output,  int64_t np, int64_t nrow, complex_array& a, int_array& ip );
void solves(complex_array& a, int_array& ip, complex_array& b, int64_t neq,
  int64_t nrh, int64_t np, int64_t n, int64_t mp, int64_t m, int64_t nop, 
  complex_array& symmetry_array);


/* Do some simple tests for integration convergence */

void test(nec_float f1r, nec_float f2r, nec_float *tr, nec_float f1i,
  nec_float f2i, nec_float *ti, nec_float dmin);

nec_float test_simple( nec_float f1r, nec_float f2r, nec_float dmin );

#ifdef NEC_ERROR_CHECK
void to_octave(nec_complex& x);
void to_octave(int& x);
void to_octave(nec_complex* a, int n, int ndim);
void to_octave(complex_array& a, int n, int ndim);
void to_octave(int* a, int n);
void to_octave(int_array& a, int n);
#endif /*  NEC_ERROR_CHECK */

#endif /* __matrix_algebra__ */
