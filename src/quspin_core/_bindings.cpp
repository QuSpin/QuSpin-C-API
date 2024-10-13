#include <_submodules/array.h>
#include <_submodules/version.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_bindings, m) {

  auto version_m = m.def_submodule("version", "Version submodule");
  define_version(version_m);

  auto array_m = m.def_submodule("array", "Array submodule");
  define_array(array_m);
}
