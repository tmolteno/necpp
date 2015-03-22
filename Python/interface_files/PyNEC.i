%module PyNEC

%include <pycomplex.swg>
%include <std_complex.i>

%{
#include "Python.h"
#include "numarray/libnumarray.h"
#include "math_util.h"
#include "nec_context.h"
#include "c_geometry.h"
#include "nec_radiation_pattern.h"
#include "nec_structure_currents.h"
#include "nec_results.h"
#include "nec_ground.h"
#include "safe_array.h"
#include "nec_exception.h"
#include <complex>
%}

/*! Exception handling stuff */

%include exception.i       
%exception
{
	try {
		$action
   	}
	catch (nec_exception* nex)
	{
		SWIG_exception(SWIG_RuntimeError,nex->get_message().c_str());
	}
   	catch (const char* message){
		SWIG_exception(SWIG_RuntimeError,message);
	}
	catch (...){
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
    	}	
}

/*! The following typemaps allow the automatic conversion of vectors and safe_arrays into numarrays */

%typemap (python, out) real_array {
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)($1.get_ptr()), tFloat64, nd, size));
}

%typemap (python, out) int_array {
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)($1.get_ptr()), tLong, nd, size));
}

%typemap (python, out) complex_array {
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)($1.get_ptr()), tComplex64, nd, size));
}

%typemap (python, out) vector<nec_float> {
	vector<double>::pointer ptr = &($1[0]);
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)ptr, tFloat64, nd, size));
}

%typemap (python, out) vector<int> {
	vector<int>::pointer ptr = &($1[0]);
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)ptr, tInt32, nd, size));
}

%typemap (python, out) vector<nec_complex> {
	vector<nec_complex>::pointer ptr = &($1[0]);
	int nd = 1;
	int size = $1.size();
	$result =(PyObject *)(NA_NewArray((void *)ptr, tComplex64, nd, size));
}

/*! The two following interface files have only been created to avoid errors during the wrapping process. */

%import "math_util.i"
%include "safe_array.i"


/*! For each of the following interface files a corresponding python file has been created. The python genrated file has been used as a starting point, then it has been improved to
	provide a more user-friendly module.
*/
%include "nec_context.i"
%include "c_geometry.i"
%include "nec_radiation_pattern.i"
%include "nec_norm_rx_pattern.i"
%include "nec_structure_excitation.i"
%include "nec_antenna_input.i"
%include "nec_near_field_pattern.i"
%include "nec_structure_currents.i"
%include "nec_ground.i"

/*The function below is added to the init function of the wrapped module.
It's mandatory to do so before to use the numpy API*/
%init %{
import_array();
%} 

