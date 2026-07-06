#pragma once

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

/** \brief LU decomposition and solve using Eigen PartialPivLU. */
void lu_decompose(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim);
void solve( int64_t n, complex_array& a, int_array& ip, complex_array& b, int64_t ndim );

/** \brief Gauss-Doolittle LU (reference implementation). */
void lu_decompose_ge(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim);
void solve_ge( int64_t n, complex_array& a, int_array& ip, complex_array& b, int64_t ndim );


void factrs(nec_output_file& s_output,  int64_t np, int64_t nrow, complex_array& a, int_array& ip );
void solves(complex_array& a, int_array& ip, complex_array& b, int64_t neq,
  int64_t nrh, int64_t np, int64_t n, int64_t mp, int64_t m, int64_t nop, 
  complex_array& symmetry_array);

#ifdef NEC_ERROR_CHECK
void to_octave(nec_complex& x);
void to_octave(int& x);
void to_octave(nec_complex* a, int n, int ndim);
void to_octave(complex_array& a, int n, int ndim);
void to_octave(int* a, int n);
void to_octave(int_array& a, int n);
#endif
