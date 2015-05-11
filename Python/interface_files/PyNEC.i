%module PyNEC

%include <pycomplex.swg>
%include <std_complex.i>

%{
#include "Python.h"
#include "numpy/arrayobject.h"
#include "src/math_util.h"
#include "src/matrix_algebra.h"
#include "src/nec_context.h"
#include "src/c_geometry.h"
#include "src/nec_radiation_pattern.h"
#include "src/nec_structure_currents.h"
#include "src/nec_results.h"
#include "src/nec_ground.h"
#include "src/safe_array.h"
#include "src/nec_exception.h"
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

/*! The following typemaps allow the automatic conversion of vectors and safe_arrays into numpy arrays */

%typemap (python, out) real_array {
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_FLOAT64, (void *)($1.data()) ));
}

%typemap (python, out) int_array {
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_INT32, (void *)($1.data()) ));
}

%typemap (python, out) complex_array {
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_COMPLEX64, (void *)($1.data()) ));
}

%typemap (python, out) vector<nec_float> {
        vector<double>::pointer ptr = &($1[0]);
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_FLOAT64, (void *)($1.data()) ));
}

%typemap (python, out) vector<int> {
        vector<int>::pointer ptr = &($1[0]);
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_INT32, (void *)($1.data()) ));
}

%typemap (python, out) vector<nec_complex> {
        vector<nec_complex>::pointer ptr = &($1[0]);
        int nd = 1;
        npy_intp size = $1.size();
        $result =(PyObject *)(PyArray_SimpleNewFromData(nd, &size, NPY_COMPLEX64, (void *)($1.data()) ));
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
It's mandatory to do so before to use the numarray API*/
%init %{
import_array();
%} 

