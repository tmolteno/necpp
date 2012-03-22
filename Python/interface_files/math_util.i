/* these typedefs can be moved to an other interface file. They have been left
here to respect the structure of the original code. */

typedef double nec_float;
typedef std::complex<nec_float> nec_complex;

typedef safe_array<int> int_array;
typedef safe_array<nec_float>  real_array;
typedef safe_array<nec_complex>  complex_array;
