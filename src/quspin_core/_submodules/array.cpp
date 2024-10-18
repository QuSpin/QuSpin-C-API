#include <algorithm>
#include <complex>
#include <variant>
#include <vector>

#include "_submodules/array.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <quspin/array/array.hpp>

#include <iostream>

namespace py = pybind11;
namespace dt = pybind11::detail;

quspin::DType numpy_to_quspin_dtype(const py::dtype &npy_dtype) {

  switch (npy_dtype.num()) {
  case dt::npy_api::NPY_INT8_:
    return quspin::Int8;
  case dt::npy_api::NPY_UINT8_:
    return quspin::UInt8;
  case dt::npy_api::NPY_INT16_:
    return quspin::Int16;
  case dt::npy_api::NPY_UINT16_:
    return quspin::UInt16;
  case dt::npy_api::NPY_INT32_:
    return quspin::Int32;
  case dt::npy_api::NPY_UINT32_:
    return quspin::UInt32;
  case dt::npy_api::NPY_INT64_:
    return quspin::Int64;
  case dt::npy_api::NPY_UINT64_:
    return quspin::UInt64;
  case dt::npy_api::NPY_FLOAT_:
    return quspin::Float;
  case dt::npy_api::NPY_DOUBLE_:
    return quspin::Double;
  case dt::npy_api::NPY_CFLOAT_:
    return quspin::CFloat;
  case dt::npy_api::NPY_CDOUBLE_:
    return quspin::CDouble;
  default:
    std::string name = py::str(py::object(npy_dtype));
    throw std::invalid_argument("dtype " + name + " not supported.");
  }
}

quspin::Array to_array(py::buffer &buf) {
  py::buffer_info info = buf.request();
  py::dtype dtype = py::dtype(info);
  auto shape = cast_vector<std::size_t>(info.shape);
  auto strides = cast_vector<std::size_t>(info.strides);

  return quspin::Array(shape, strides, numpy_to_quspin_dtype(dtype), info.ptr);
}

const quspin::Array to_array(const py::buffer &buf) {
  py::buffer_info info = buf.request();
  py::dtype dtype = py::dtype(info);
  auto shape = cast_vector<std::size_t>(info.shape);
  auto strides = cast_vector<std::size_t>(info.strides);

  return quspin::Array(shape, strides, numpy_to_quspin_dtype(dtype), info.ptr);
}

template <typename T> std::string quspin_format() {
  char dtype_char = py::dtype::of<T>().char_();
  return std::string(&dtype_char, 1);
}

py::buffer_info to_numpy_style_buffer(quspin::Array &arr) {
  return std::visit(
      [](auto &&arr) {
        using T = typename std::decay_t<decltype(arr)>::value_type;

        auto shape = cast_vector<py::ssize_t>(arr.shape());
        auto strides = cast_vector<py::ssize_t>(arr.strides());

        return py::buffer_info(arr.mut_data(), sizeof(T), quspin_format<T>(),
                               arr.ndim(), shape, strides);
      },
      arr.get_variant_obj());
}

void define_array(py::module_ &array_m) {

  py::class_<quspin::Array>(array_m, "Array", py::buffer_protocol())
      .def(py::init([](py::buffer buf) { return to_array(buf); }))
      .def_buffer([](quspin::Array &arr) -> py::buffer_info {
        return to_numpy_style_buffer(arr);
      })
      .def_property_readonly(
          "dtype",
          [](quspin::Array &arr) {
            return std::visit(
                [](const auto &dtype) {
                  using T = typename std::decay_t<decltype(dtype)>::value_type;
                  return py::dtype::of<T>();
                },
                arr.dtype().get_variant_obj());
          })
      .def_property_readonly("shape",
                             [](quspin::Array &arr) {
                               auto tup = py::tuple(arr.ndim());
                               for (std::size_t i = 0; i < arr.ndim(); i++) {
                                 tup[i] = arr.shape(i);
                               }
                               return tup;
                             })
      .def_property_readonly("strides",
                             [](quspin::Array &arr) {
                               auto tup = py::tuple(arr.ndim());
                               for (std::size_t i = 0; i < arr.ndim(); i++) {
                                 tup[i] = arr.strides(i);
                                 ;
                               }
                               return tup;
                             })
      .def_property_readonly("ndim", &quspin::Array::ndim)
      .def_property_readonly("size", &quspin::Array::size);
}
