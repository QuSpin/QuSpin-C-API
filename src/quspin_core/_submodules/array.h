#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <quspin/array/array.hpp>

namespace py = pybind11;

template <typename To, typename From>
std::vector<To> cast_vector(const std::vector<From> &in) {
  std::vector<To> &&vec = std::vector<To>(in.size());
  std::transform(in.begin(), in.end(), vec.begin(),
                 [](const From &val) { return static_cast<To>(val); });
  return std::move(vec);
}

quspin::DType numpy_to_quspin_dtype(const py::dtype &);
quspin::Array to_array(py::buffer &);
const quspin::Array to_array(const py::buffer &);

void define_array(py::module_ &);
