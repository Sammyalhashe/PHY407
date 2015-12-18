#ifndef CONVERT_HPP
#define CONVERT_HPP
#include <Python.h>
#include "QickArray.hpp"
#include "walker.hpp"

class Convert_QickArray{

public:
  Convert_QickArray();
  // borrow data:
  PyObject*  my2py(QickArray& a);
  QickArray* py2my(PyObject* a);
  // copy data:
  PyObject*  my2py_copy(QickArray& a);
  QickArray* py2my_copy(PyObject* a);
  // print array:
  void dump(QickArray& a);

};

#endif
