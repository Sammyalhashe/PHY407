/* file: MontePython_cxx.i */
%module MontePython_cxx
%{
#include "../C++/walker.hpp"
#include "../C++/warray.hpp"
#include "../C++/functions.hpp"
#include "../C++/random.hpp"
#include "../C++/convert.hpp"
#include "../C++/mc.hpp"
#include "../C++/RngStream.h"
#include <Numeric/arrayobject.h>
%}

%feature("autodoc","1"); //requires swig >= 1.3.23

%include "../C++/walker.hpp" 
%include "../C++/warray.hpp" 
%include "../C++/functions.hpp"
%include "../C++/random.hpp"
%include "../C++/convert.hpp"
%include "../C++/mc.hpp"
%include "../C++/RngStream.h"

%extend Warray {

  PyObject* pyInitialize(int no_of_walkers_, int particles_, int dim_,
		    PyObject* a_){

    import_array();

    PyArrayObject* a = (PyArrayObject*) 
      PyArray_ContiguousFromObject(a_, PyArray_DOUBLE, 0, 0);

    self->initialize(no_of_walkers_,particles_,dim_,
		     (double*)a->data,a->dimensions[0]);

    return PyArray_Return(a);
  }

  PyObject* pyWalkers2grid3D(PyObject* a_, double dx, double dy, double dz){

    PyArrayObject* a = (PyArrayObject*)
      PyArray_ContiguousFromObject(a_, PyArray_DOUBLE, 0, 0);
    if (a == NULL) { return NULL; }

    QickArray* ma = new QickArray((double*) a->data, a->nd, a->dimensions);
    self->walkers2grid3D(*ma,dx,dy,dz);

    return PyArray_Return(a);

  }

  PyObject* pyWalkers2grid3DPure(PyObject* a_, 
				 double dx, double dy, double dz,
				 double weight){

    PyArrayObject* a = (PyArrayObject*)
      PyArray_ContiguousFromObject(a_, PyArray_DOUBLE, 0, 0);
    if (a == NULL) { return NULL; }

    QickArray* ma = new QickArray((double*) a->data, a->nd, a->dimensions);
    self->walkers2grid3DPure(*ma,dx,dy,dz,weight);

    return PyArray_Return(a);

  }

}
