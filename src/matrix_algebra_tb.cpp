#include "catch.hpp"

#include "matrix_algebra.h"
#include "nec_output.h"

/**
 * Set up a matrix for inversion. The check is done with to_octave
 * 
 * a = [3 + 0I, 1 + 0I, -4 + 0I, 2 + 0I; 3 + 0I, 1 + 0I, 0 + 0I, 2 + 0I; 2 + 0I, 13 + 0I, -1 + 0I, 0 + 0I; -2 + 0I, 3 + 0I, -1 + 0I, 4 + 0I];
 * [l,u,p] = lu(a);
 *  
 * check = inv(p)*l*u
 * a == check
 * 
 * the output should equal U+L-eye(size(a))
 * */
void matrix_setup(complex_array& A);
void matrix_setup(complex_array& A) {
  int n = 4;
  A(0,0) = 3.0;
  A(0,1) =1.0;
  A(0,2) = -4.0;
  A(0,3) =2.0;
  
  A(1,0) =3.0;
  A(1,1) =1.0;
  A(1,2) =0.0;
  A(1,3) =2.0;
  
  A(2,0) =2.0;
  A(2,1) =13.0;
  A(2,2) = -1.0;
  A(2,3) =0.0;
  
  A(3,0) = -2.0;
  A(3,1) =3.0;
  A(3,2) = -1.0;
  A(3,3) =4.0;
  for (int i = 1; i < n; i++ )  {
    for (int j = 0; j < i; j++ )
      std::swap(A(i,j),A(j,i));
  }
}

#define REQUIRE_APPROX_EQUAL(a, b) { \
  static nec_float eps = 1e-4; \
  REQUIRE(std::abs(a - b) < eps); }


TEST_CASE( "LU Decomposition Gauss-Doolittle", "[lu_decompose_ge]") {
    int N=4;
    nec_output_file s_output;
    complex_matrix A(N,N);
    int_array piv(N);
    
    matrix_setup(A);
    std::cout << "a = ";
    to_octave(A,N,N);
    
    lu_decompose_ge(s_output, N, A, piv, N);
    std::cout << "lu = ";
    to_octave(A,N,N);
    std::cout << "piv = ";
    to_octave(piv,N);

    REQUIRE_APPROX_EQUAL(A(0,0), 3.0);
    REQUIRE_APPROX_EQUAL(A(0,1), 1.0);
    REQUIRE_APPROX_EQUAL(A(0,2), 0.0);
    REQUIRE_APPROX_EQUAL(A(0,3), 2.0);

    REQUIRE_APPROX_EQUAL(A(1,0), 1.0);
    REQUIRE_APPROX_EQUAL(A(1,1), 12.3333);
    REQUIRE_APPROX_EQUAL(A(1,2), -1.0);
    REQUIRE_APPROX_EQUAL(A(1,3), -1.3333);

    REQUIRE_APPROX_EQUAL(A(2,0), 0.66667);
    REQUIRE_APPROX_EQUAL(A(2,1), 0.0);
    REQUIRE_APPROX_EQUAL(A(2,2), -4.0);
    REQUIRE_APPROX_EQUAL(A(2,3), 0.0);

    REQUIRE_APPROX_EQUAL(A(3,0), -0.66667);
    REQUIRE_APPROX_EQUAL(A(3,1), 0.297297);
    REQUIRE_APPROX_EQUAL(A(3,2), 0.175676);
    REQUIRE_APPROX_EQUAL(A(3,3), 5.72973);
}


TEST_CASE( "LU Decomposition (Eigen)", "[lu_decompose]") {
    // Fill a 4x4 matrix with the same values as matrix_setup(),
    // transpose for NEC storage, factor, solve, and verify.
    // Uses the same pattern as the 12x12 factrs/solves test.
    int N=4;
    nec_output_file s_output;
    complex_matrix A(N,N);
    int_array piv(N);
    
    // Fill conceptual A (same values as matrix_setup)
    A(0,0) = 3.0;  A(0,1) = 1.0;  A(0,2) = -4.0;  A(0,3) = 2.0;
    A(1,0) = 3.0;  A(1,1) = 1.0;  A(1,2) =  0.0;  A(1,3) = 2.0;
    A(2,0) = 2.0;  A(2,1) =13.0;  A(2,2) = -1.0;  A(2,3) = 0.0;
    A(3,0) =-2.0;  A(3,1) = 3.0;  A(3,2) = -1.0;  A(3,3) = 4.0;

    // Known solution x
    complex_array x(N);
    for (int i = 0; i < N; i++)
        x[i] = nec_complex(i + 1, -(i + 1) * 0.5);

    // Compute b = A * x (conceptual, before transpose)
    complex_array b(N);
    for (int i = 0; i < N; i++) {
        b[i] = nec_complex(0.0, 0.0);
        for (int j = 0; j < N; j++)
            b[i] += A(i, j) * x[j];
    }

    // Transpose for NEC column-major storage convention
    for (int i = 1; i < N; i++)
        for (int j = 0; j < i; j++)
            std::swap(A(i, j), A(j, i));

    // Factor (modifies A in place) and solve (modifies b in place)
    lu_decompose(s_output, N, A, piv, N);
    solve(N, A, piv, b, N);

    // Verify solution matches x
    for (int i = 0; i < N; i++) {
        REQUIRE_APPROX_EQUAL(b[i].real(), x[i].real());
        REQUIRE_APPROX_EQUAL(b[i].imag(), x[i].imag());
    }
}


TEST_CASE( "LU Decomposition (Eigen self-check)", "[lu_decompose]") {
  using namespace Eigen;
  MatrixXcd A(4,4);

  A(0,0) = 3.0;  A(0,1) = 1.0;  A(0,2) = -4.0;  A(0,3) = 2.0;
  A(1,0) = 3.0;  A(1,1) = 1.0;  A(1,2) =  0.0;  A(1,3) = 2.0;
  A(2,0) = 2.0;  A(2,1) =13.0;  A(2,2) = -1.0;  A(2,3) = 0.0;
  A(3,0) =-2.0;  A(3,1) = 3.0;  A(3,2) = -1.0;  A(3,3) = 4.0;

  PartialPivLU<MatrixXcd> lu(A);
  MatrixXcd l = MatrixXcd::Identity(4,4);
  l.triangularView<StrictlyLower>() = lu.matrixLU();
  MatrixXcd u = lu.matrixLU().triangularView<Upper>();
  MatrixXcd B = lu.permutationP().inverse() * l * u;

  for (int i=0;i<4;i++) {
    for (int j=0; j<4; j++) {
      REQUIRE_APPROX_EQUAL(B(i,j), A(i,j));
    }
  }
}

/*-----------------------------------------------------------------------*/

/*
 * Comprehensive test: factor a complex 12×12 matrix via factrs(),
 * solve Ax=b via solves(), and verify the solution reconstructs correctly.
 *
 * This exercises the full pipeline: transpose handling, LU decomposition,
 * partial pivoting, forward/back substitution, symmetry dispatch (nop=1,
 * the common no-symmetry path), all with complex arithmetic.
 */
TEST_CASE( "LU Factor+Solve 12x12 complex (factrs/solves)", "[factrs_solves]") {
    const int N = 12;
    nec_output_file s_output;
    complex_matrix A(N, N);
    int_array piv(N);

    // --- Build a well-conditioned 12×12 complex matrix ---
    // Diagonal-dominant to prevent near-singular pivots.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                // Strong diagonal: magnitude grows with row index
                A(i, j) = nec_complex(3.0 * (i + 1), 1.0 * (i + 1));
            } else {
                // Off-diagonal: small magnitude, deterministic pattern
                nec_float re = 0.1 * (i - j);
                nec_float im = 0.05 * ((i + 1) * (j + 1)) / (N * N);
                A(i, j) = nec_complex(re, im);
            }
        }
    }

    // --- Known solution x (deterministic, all entries non-zero) ---
    complex_array x(N);
    for (int i = 0; i < N; i++) {
        x[i] = nec_complex(i + 1, -(i + 1) * 0.5);
    }

    // --- Compute b = A * x using conceptual A (before transpose) ---
    complex_array b(N);
    for (int i = 0; i < N; i++) {
        b[i] = nec_complex(0.0, 0.0);
        for (int j = 0; j < N; j++) {
            b[i] += A(i, j) * x[j];
        }
    }

    // --- Transpose A for NEC column-major storage convention ---
    // (same as matrix_setup() does for the 4×4 LU tests)
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < i; j++) {
            std::swap(A(i, j), A(j, i));
        }
    }

    // --- Factor and solve ---
    // nop=1, mp=m=0 (no symmetry, no patches) — the common path
    factrs(s_output, N, N, A, piv);

    complex_array sym(1);  // dummy symmetry array, unused when nop==1
    solves(A, piv, b, N, 1, N, N, 0, 0, 1, sym);

    // --- Verify solution: b now contains the computed x', should match x ---
    for (int i = 0; i < N; i++) {
        REQUIRE_APPROX_EQUAL(b[i].real(), x[i].real());
        REQUIRE_APPROX_EQUAL(b[i].imag(), x[i].imag());
    }
}
