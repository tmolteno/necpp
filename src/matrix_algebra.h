#ifndef __matrix_algebra__
#define __matrix_algebra__

/*
	Copyright (C) 2004  Timothy C.A. Molteno
	
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

void lu_decompose(nec_output_file& s_output, int n, complex_array& a, int_array& ip, int ndim);
void factrs(nec_output_file& s_output,  int np, int nrow, complex_array& a, int_array& ip );
void solve( int n, complex_array& a, int_array& ip, complex_array& b, int ndim );

void solves(complex_array& a, int_array& ip, complex_array& b, int neq,
	int nrh, int np, int n, int mp, int m, int nop, 
	complex_array& symmetry_array);


/* Do some simple tests for integration convergence */

void 	test(nec_float f1r, nec_float f2r, nec_float *tr, nec_float f1i,
	nec_float f2i, nec_float *ti, nec_float dmin);

nec_float test_simple( nec_float f1r, nec_float f2r, nec_float dmin );

#endif /* __matrix_algebra__ */
