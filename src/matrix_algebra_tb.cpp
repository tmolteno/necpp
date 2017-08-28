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


TEST_CASE( "LU Decomposition LAPACK", "[lu_decompose]") {
    int N=4;
    nec_output_file s_output;
    complex_matrix A(N,N);
    int_array piv(N);
    
    matrix_setup(A);
    std::cout << "a = ";
    to_octave(A,N,N);
    
    lu_decompose(s_output, N, A, piv, N);
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

#if USING_EIGEN_3VECT

#include <Eigen/Dense>
using namespace Eigen;
#include <iostream>
using namespace std;

TEST_CASE( "LU Decomposition EIGEN", "[lu_decompose]") {
  
  MatrixXcd A(4,4);

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

  cout << "Here is the matrix A:" << endl << A << endl;
  
  Eigen::PartialPivLU<MatrixXcd> lu(A);
  cout << "Here is, up to permutations, its LU decomposition matrix:"
  << endl << lu.matrixLU() << endl;
  cout << "Here is the L part:" << endl;
  MatrixXcd l = MatrixXcd::Identity(4,4);
  l.triangularView<StrictlyLower>() = lu.matrixLU();
  cout << l << endl;
  cout << "Here is the U part:" << endl;
  MatrixXcd u = lu.matrixLU().triangularView<Upper>();
  cout << u << endl;
  cout << "Let us now reconstruct the original matrix m:" << endl;
  MatrixXcd B = lu.permutationP().inverse() * l * u;
  cout << B << endl;
  for (int i=0;i<4;i++) {
    for (int j=0; j<4; j++) {
      REQUIRE_APPROX_EQUAL(B(i,j), A(i,j));
    }
  }
}

#endif /* #if USING_EIGEN_3VECT */
