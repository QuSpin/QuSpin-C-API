#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_array(py::module_ &array_m) {
  array_m.def("hello", []() { return "Hello, World!"; });
}
