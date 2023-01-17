#ifndef __QUSPIN_CORE_NUMPY_INTERFACE_H__
#define __QUSPIN_CORE_NUMPY_INTERFACE_H__

namespace quspin_core_abi {

inline NPY_TYPES npy_typenum(PyArrayObject * array){
  return static_cast<NPY_TYPES>(PyArray_TYPE(array));
}

inline void * npy_data(PyArrayObject * array){
  return PyArray_DATA(array);
}


}

#endif