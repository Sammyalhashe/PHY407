#include "convert.hpp"
#include <Numeric/arrayobject.h>

Convert_QickArray:: Convert_QickArray()
{
  import_array();
}

PyObject* Convert_QickArray:: my2py(QickArray& a)
{
  //a.dump(std::cout);
  
  PyArrayObject* array =  (PyArrayObject*) \
          PyArray_FromDimsAndData(a.get_dimension_info(), 
				  a.get_dimensions_ptr(), 
				  PyArray_DOUBLE, (char*) a.get_data_ptr()); 
  if (array == NULL) {
    return NULL; /* PyArray_FromDimsAndData raised an exception */ 
  }
  return PyArray_Return(array);
}

PyObject* Convert_QickArray:: my2py_copy(QickArray& a)
{
  //a.dump(std::cout);
  PyArrayObject* array =  (PyArrayObject*) \
          PyArray_FromDims(a.get_dimension_info(), a.get_dimensions_ptr(), 
			   PyArray_DOUBLE); 
  if (array == NULL) {
    return NULL; /* PyArray_FromDims raised an exception */ 
  }
  double* ad = (double*) array->data;
  double* A = a.get_data_ptr();
  for (int i = 0; i < a.no_of_elements(); i++) {
    ad[i] = A[i];
  }
  return PyArray_Return(array);
}

QickArray* Convert_QickArray:: py2my(PyObject* a_)
{
  PyArrayObject* a = (PyArrayObject*) 
    PyArray_ContiguousFromObject(a_, PyArray_DOUBLE, 0, 0);
  if (a == NULL) { return NULL; }
  /*
  PyArrayObject* a = (PyArrayObject*) a_;
  if (a->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "Not an array of double");
    return NULL;
  }
  */
  // borrow the data, but wrap it in MyArray:
  QickArray* ma = new QickArray((double*) a->data, a->nd, a->dimensions);
  return ma;
}

QickArray* Convert_QickArray:: py2my_copy(PyObject* a_)
{
  PyArrayObject* a = (PyArrayObject*) 
    PyArray_ContiguousFromObject(a_, PyArray_DOUBLE, 0, 0);
  if (a == NULL) { return NULL; }
  /*
  PyArrayObject* a = (PyArrayObject*) a_;
  if (a->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "Not an array of double");
    return NULL;
  }
  */
  QickArray* ma = new QickArray();
  if (a->nd == 1) {
    ma->redim(a->dimensions[0]);
  } else if (a->nd == 2) {
    ma->redim(a->dimensions[0], a->dimensions[1]);
  } else if (a->nd == 3) {
    ma->redim(a->dimensions[0], a->dimensions[1], a->dimensions[2]);
  } else if (a->nd == 4) {
    ma->redim(a->dimensions[0], a->dimensions[1], a->dimensions[2],
	      a->dimensions[3]);
  } else if (a->nd == 5) {
    ma->redim(a->dimensions[0], a->dimensions[1], a->dimensions[2],
	      a->dimensions[3], a->dimensions[4]);
  } else if (a->nd == 6) {
    ma->redim(a->dimensions[0], a->dimensions[1], a->dimensions[2],
	      a->dimensions[3], a->dimensions[4], a->dimensions[5]);
  }
  // copy data:
  double* ad = (double*) a->data;
  double* mad = ma->get_data_ptr();
  for (int i = 0; i < ma->no_of_elements(); i++) {
    mad[i] = ad[i];
  }
  return ma;
}

void Convert_QickArray:: dump(QickArray& a)
{
  std::cout << a.python_string() << std::endl;
}
