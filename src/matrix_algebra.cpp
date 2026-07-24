/*
  NEC-2 matrix algebra: LU decomposition (Eigen PartialPivLU),
  forward/back substitution, symmetry handling (factrs/solves),
  matrix assembly helpers.

  Matrix storage: column-major transposed (NEC convention).
  LU uses Eigen::Map with OuterStride for zero-copy factorization.
  Gauss-Doolittle reference implementation retained for testing.
*/
#include "math_util.h"
#include <iostream>
using namespace std;

#include "nec_exception.h"
#include "matrix_algebra.h"
#include "nec_output.h"

#include "common.h"
#include <vector>

#ifdef NEC_ERROR_CHECK

void to_octave(nec_complex& x) {
    cout << real(x) << " + " << imag(x) << "I";
}

void to_octave(int& x) {
    cout << x;
}

void to_octave(nec_complex* a, int n, int ndim) {
    cout << "[";
    for (int row = 0; row < n; row++ ) {
        int64_t row_offset = row*ndim;
        for (int i = 0; i < n; i++ ) {
            to_octave(a[i+row_offset]);
            if (i < n-1)
                cout << ", ";
        }
        if (row < n-1)
            cout << "; ";
    }
    cout << "];" << endl;
}

void to_octave(complex_array& a, int n, int ndim) {
    to_octave(a.data(),n,ndim);
}

void to_octave(int* a, int n) {
    cout << "[";
    for (int i = 0; i < n; i++ ) {
        to_octave(a[i]);
        if (i < n-1)
            cout << ", ";
    }
    cout << "];" << endl;
}

void to_octave(int_array& a, int n) {
    to_octave(a.data(),n);
}
#endif /* NEC_ERROR_CHECK */



/*! \brief Solve The system of equations using Gaussian Elimination.
    Subroutine to factor a matrix into a unit lower triangular matrix 
    and an upper triangular matrix using the Gauss-Doolittle algorithm 
    presented on pages 411-416 of A. Ralston -- a first course in 
    numerical analysis.
    
    Comments below refer to comments in Ralstons text.
    
    (matrix transposed.)
*/
void lu_decompose_ge(nec_output_file& s_output, int64_t n, complex_array& a, int_array& ip, int64_t ndim)
{
    DEBUG_TRACE("lu_decompose_ge(" << n << "," << ndim << ")");
    
#ifdef NEC_MATRIX_CHECK
    // Debug output for matrix verification
    cout << "a = ";
    to_octave(a,n,ndim);
#endif
    
    /* Allocate scratch memory */
    complex_array scm;
    scm.resize(n);
    
    /* Un-transpose the matrix for Gauss elimination */
    for (int i = 1; i < n; i++ ) {
        int64_t i_offset = i * ndim;
        int64_t j_offset = 0;
        for (int j = 0; j < i; j++ ) {
            nec_complex aij = a[i+j_offset];
            a[i+j_offset] = a[j+i_offset];
            a[j+i_offset] = aij;
            
            j_offset += ndim;
        }
    }
    
    bool iflg=false;
    /* step 1 */
    for (int r = 0; r < n; r++ ) {
        int64_t r_offset = r*ndim;
        
        for (int k = 0; k < n; k++ )
            scm[k]= a[k+r_offset];
        
        /* steps 2 and 3 */
        int rm1 = r;
        for (int j = 0; j < rm1; j++ ) {
            int pj= ip[j]-1;
            nec_complex arj = scm[pj];
            a[j+r_offset]= arj;
            scm[pj]= scm[j];
            int jp1 = j+1;

            int64_t j_offset = j*ndim;
            for (int i = jp1; i < n; i++ )
                scm[i] -= a[i+j_offset]* arj;

        }
        
        /* step 4 */
        nec_float dmax = norm(scm[r]);
        
        int rp1 = r+1;
        ip[r]= rp1;
        for (int i = rp1; i < n; i++ ) {
            nec_float elmag = norm(scm[i]);
            if ( elmag >= dmax) {
                dmax = elmag;
                ip[r] = i+1;    // set the permute array element
            }
        }
        
        if ( dmax < 1.e-10)
            iflg=true;
        
        int pr = ip[r]-1;
        a[r+r_offset] = scm[pr];
        scm[pr] = scm[r];
        
        /* step 5 */
        if ( rp1 < n) {
            nec_complex arr = cplx_10() / a[r+r_offset];
            
            for (int i = rp1; i < n; i++ )
                a[i+r_offset]= scm[i]* arr;
        }
        
        if ( true == iflg ) {
            s_output.string("\n    PIVOT(");
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


/*! \brief Solve system of linear equations
        Subroutine to solve the matrix equation lu*x=b where l is a unit
        lower triangular matrix and u is an upper triangular matrix both
        of which are stored in a.    The RHS vector b is input and the
        solution is returned through vector b.     (matrix transposed)
*/
void solve_ge( int64_t n, complex_array& a, int_array& ip,
        complex_array& b, int64_t ndim )
{
    DEBUG_TRACE("solve(" << n << "," << ndim << ")");
    complex_array y(n);
    
    /* forward substitution */
    for (int64_t i = 0; i < n; i++ ) {
        int64_t pivot_index = ip[i]-1;
        y[i]= b[pivot_index];
        b[pivot_index]= b[i];
        int64_t ip1= i+1;
        
        int64_t i_offset = i*ndim;
        for (int64_t j = ip1; j < n; j++ )
            b[j] -= a[j+i_offset] * y[i];
    }
    
    /* backward substitution */
    for (int64_t k = 0; k < n; k++ ) {
        int64_t i= n-k-1;
        nec_complex sum(cplx_00());
        int64_t ip1= i+1;
        
        for (int64_t j = ip1; j < n; j++ )
            sum += a[i+j*ndim]* b[j];
        
        b[i] = (y[i]- sum) / a[i+i*ndim];
    }
}


/*! \brief Use Eigen PartialPivLU to perform LU decomposition.
    The input matrix a_in is in NEC column-major transposed format.
    On return, a_in contains L+U factors and ip contains 1-indexed pivots. */
void lu_decompose(nec_output_file& s_output, int64_t n, complex_array& a_in, int_array& ip, int64_t ndim)
{
    UNUSED(s_output);
    DEBUG_TRACE("lu_decompose(" << n << "," << ndim << ")");
    ASSERT(n <= ndim);

#ifdef NEC_MATRIX_CHECK
    cout << "eigen_a = ";
    to_octave(a_in,n,ndim);
#endif

    /* Un-transpose the matrix for Gauss elimination */
    for (int i = 1; i < n; i++ ) {
        int64_t i_offset = i * ndim;
        int64_t j_offset = 0;
        for (int64_t j = 0; j < i; j++ ) {
            nec_complex aij = a_in[i+j_offset];
            a_in[i+j_offset] = a_in[j+i_offset];
            a_in[j+i_offset] = aij;
            j_offset += ndim;
        }
    }

    /* Factor using Eigen's PartialPivLU through an in-place Map.
       Both nec_complex and Eigen use std::complex<double>, and both
       use column-major storage, so the memory layout is identical. */
    using MatrixXcd = Eigen::Matrix<std::complex<nec_float>, Eigen::Dynamic, Eigen::Dynamic>;
    Eigen::Map<MatrixXcd, 0, Eigen::OuterStride<> > A_map(
        reinterpret_cast<std::complex<nec_float>*>(a_in.data()), n, n,
        Eigen::OuterStride<>(static_cast<int>(ndim)));
    
    Eigen::PartialPivLU<MatrixXcd> lu(A_map);
    
    /* Copy LU factors back into a_in (modifies through the Map). */
    A_map = lu.matrixLU();

    /* Store Eigen's permutation indices in ip (1-based for compatibility).
       solve will use these directly instead of the GE pivot format. */
    auto eigen_indices = lu.permutationP().indices();
    for (int64_t i = 0; i < n; i++) {
        ip[i] = static_cast<int32_t>(eigen_indices[i] + 1);
    }

#ifdef NEC_MATRIX_CHECK
    cout << "eigen_solved = ";
    to_octave(a_in,n,ndim);
    cout << "eigen_ip = ";
    to_octave(ip,n);
#endif
}
    

/* solve reconstructs Eigen PartialPivLU from stored LU and permutation,
   then solves via forward/back substitution. */
void solve( int64_t n, complex_array& a, int_array& ip,
    complex_array& b, int64_t ndim )
{
    DEBUG_TRACE("solve(" << n << "," << ndim << ")");
    
    using MatrixXcd = Eigen::Matrix<std::complex<nec_float>, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXcd = Eigen::Matrix<std::complex<nec_float>, Eigen::Dynamic, 1>;
    
    Eigen::Map<MatrixXcd, 0, Eigen::OuterStride<> > LU(
        reinterpret_cast<std::complex<nec_float>*>(a.data()), n, n,
        Eigen::OuterStride<>(static_cast<int>(ndim)));
    Eigen::Map<VectorXcd> rhs(reinterpret_cast<std::complex<nec_float>*>(b.data()), n);
    
    /* Rebuild the permutation from ip (1-based Eigen indices). */
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(n);
    for (int64_t i = 0; i < n; i++)
        P.indices()[i] = static_cast<int>(ip[i] - 1);
    
    /* Apply permutation: pb = P * rhs */
    VectorXcd pb = P * rhs;
    
    /* Forward substitution: solve L * z = pb */
    VectorXcd z = pb;
    for (int64_t i = 0; i < n; i++) {
        for (int64_t j = i + 1; j < n; j++)
            z[j] -= LU(j, i) * z[i];
    }
    
    /* Backward substitution: solve U * rhs_out = z */
    for (int64_t k = 0; k < n; k++) {
        int64_t i = n - k - 1;
        std::complex<nec_float> sum(0.0, 0.0);
        for (int64_t j = i + 1; j < n; j++)
            sum += LU(i, j) * rhs[j];
        rhs[i] = (z[i] - sum) / LU(i, i);
    }
}
/*-----------------------------------------------------------------------*/

/** factrs

    For symmetric structure, transforms submatricies to form
    matricies of the symmetric modes and calls routine to LU decompose
    matricies.
    
    If no symmetry [nrow = np], the routine is called to LU decompose the
    complete matrix.
*/
void factrs(nec_output_file& s_output,    int64_t np, int64_t nrow, complex_array& a, int_array& ip )
{
    DEBUG_TRACE("factrs(" << np << "," << nrow << ")");
    if (nrow == np) { // no symmetry: skip segment()/loop overhead
        lu_decompose(s_output,    np, a, ip, nrow );
        return;
    }
    
    int num_symmetric_modes = static_cast<int>(nrow / np);
    DEBUG_TRACE("\tnum_symmetric_modes = " << num_symmetric_modes);
    
    for (int mode = 0; mode < num_symmetric_modes; mode++ ) {
        int64_t mode_offset = mode * np;
        
        complex_array a_temp = a.eigen_segment(mode_offset, a.size()-mode_offset);
        int_array ip_temp = ip.eigen_segment(mode_offset, ip.size()-mode_offset);
        
        lu_decompose(s_output,    np, a_temp, ip_temp, nrow );
    }
}

/*-----------------------------------------------------------------------*/



/**
    Subroutine solves, for symmetric structures, handles the
    transformation of the right hand side vector and solution
    of the matrix eq.
    \param neq number of equations?
    \param nrh dimension of right hand    vector?
*/
void solves(complex_array& a, int_array& ip, complex_array& b, int64_t neq,
            int64_t nrh, int64_t np, int64_t n, int64_t mp, int64_t m, int64_t nop, 
            complex_array& symmetry_array)
{
    DEBUG_TRACE("solves(" << neq << "," << nrh << "," << np << "," << n << ")");
    DEBUG_TRACE("            ( nop=" << nop << ")");
    
    /* Allocate some scratch memory */
    complex_array scm;
    scm.resize(n + 2*m);
    
    int64_t npeq= np+ 2*mp;
    nec_float fnop = nec_float(nop);
    nec_float fnorm = 1.0/ fnop;
    int64_t nrow= neq;
    
    if ( nop != 1) {
        for (int ic = 0; ic < nrh; ic++ ) {
            int64_t column_offset = ic*neq;
            if ( (n != 0) && (m != 0) ) {
                for (int64_t i = 0; i < neq; i++ )
                    scm[i]= b[i+column_offset];

                int64_t j= np-1;

                for (int64_t k = 0; k < nop; k++ ) {
                    if ( k != 0 ) {
                        int64_t ia= np-1;
                        for (int64_t i = 0; i < np; i++ ) {
                            ia++;
                            j++;
                            b[j+column_offset]= scm[ia];
                        }

                        if ( k == (nop-1) )
                            continue;
                    } /* if ( k != 0 ) */
    
                    int64_t mp2 = 2*mp;
                    int64_t ib= n-1;
                    for (int64_t i = 0; i < mp2; i++ ) {
                        ib++;
                        j++;
                        b[j+column_offset]= scm[ib];
                    }
                } /* for( k = 0; k < nop; k++ ) */
            } /* if ( (n != 0) && (m != 0) ) */

            /* transform matrix eq. rhs vector according to symmetry modes */
            for (int64_t i = 0; i < npeq; i++ ) {
                nec_complex sum_normal(b[i+column_offset]);
                for (int64_t k = 1; k < nop; k++ ) {
                    int64_t ia= i+ k* npeq;
                    scm[k]= b[ia+column_offset];
                    sum_normal += scm[k];
                }

                b[i+column_offset]= sum_normal * fnorm;

                for (int64_t k = 1; k < nop; k++ ) {
                    int64_t ia= i+ k* npeq;
                    nec_complex sum(scm[0]);
    
                    for (int64_t j = 1; j < nop; j++ )
                        sum += scm[j]* conj( symmetry_array[k+j*nop]);
    
                    b[ia+column_offset]= sum* fnorm;
                }
            } /* for( i = 0; i < npeq; i++ ) */
        } /* for( ic = 0; ic < nrh; ic++ ) */
    } /* if ( nop != 1) */

    /* solve each mode equation */
    for (int64_t kk = 0; kk < nop; kk++ ) {
        int64_t ia= kk* npeq;

        for (int64_t ic = 0; ic < nrh; ic++ ) {
            int64_t column_offset = ic*neq;
            complex_array a_sub = a.eigen_segment(ia, a.size()-ia);
            complex_array b_sub = b.eigen_segment(ia+column_offset, b.size() - (ia+column_offset) );
            int_array ip_sub = ip.eigen_segment(ia, ip.size()-ia);
            solve( npeq, a_sub, ip_sub, b_sub, nrow );
        }
    } /* for( kk = 0; kk < nop; kk++ ) */
    
    if ( nop == 1) {
        return;
    }
    
    /* inverse transform the mode solutions */
    for (int64_t ic = 0; ic < nrh; ic++ ) {
        int64_t column_offset = ic*neq;
        for (int64_t i = 0; i < npeq; i++ ) {
            nec_complex sum_normal(b[i+column_offset]);
            for (int64_t k = 1; k < nop; k++ ) {
                int64_t ia= i+ k* npeq;
                scm[k]= b[ia+column_offset];
                sum_normal += scm[k];
            }

            b[i+column_offset]= sum_normal;
            
            for (int64_t k = 1; k < nop; k++ ) {
                int64_t ia= i+ k* npeq;
                
                nec_complex sum(scm[0]);

                for (int64_t j = 1; j < nop; j++ )
                    sum += scm[j]* symmetry_array[k+j*nop];

                b[ia+column_offset]= sum;
            }
        } /* for( i = 0; i < npeq; i++ ) */

        if ( (n == 0) || (m == 0) )
            continue;

        for (int64_t i = 0; i < neq; i++ )
                        scm[i]= b[i+column_offset];

        int64_t j = np-1;

        for (int64_t k = 0; k < nop; k++ ) {
            if ( k != 0 ) {
                int64_t ia = np-1;
                for (int64_t i = 0; i < np; i++ ) {
                    ia++;
                    j++;
                    b[ia+column_offset]= scm[j];
                }

                if ( k == (nop-1) )
                    continue;
            } /* if ( k != 0 ) */

            int64_t ib = n-1;
            int64_t mp2 = 2* mp;
            for (int64_t i = 0; i < mp2; i++ ) {
                ib++;
                j++;
                b[ib+column_offset]= scm[j];
            }
        } /* for( k = 0; k < nop; k++ ) */
    } /* for( ic = 0; ic < nrh; ic++ ) */
}

