#include <pybind11/pybind11.h>

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define VERSION_STR STR(VERSION_INFO)

namespace py = pybind11;

void define_version(py::module_ &version_m) {
  version_m.def("_get_version", []() { return VERSION_STR; });
}
