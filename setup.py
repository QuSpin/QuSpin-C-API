# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import numpy as np
import os
import sys

__version__ = "0.1.0"

extra_compile_args = []
extra_link_args = []

if os.environ.get("COVERAGE", False):
    if sys.platform != "linux":
        raise ValueError("Coverage is only supported on Linux")

    extra_compile_args += [
        "-O0",
        "--coverage",
        "-fno-inline",
        "-fno-inline-small-functions",
        "-fno-default-inline",
    ]
    extra_link_args += ["--coverage"]


if sys.platform == "darwin":
    extra_compile_args += ["-mmacosx-version-min=10.15"]

ext_modules = [
    Pybind11Extension(
        "hamiltonian_core.ext",
        [os.path.join("src", "hamiltonian_core", "ext.cxx")],
        define_macros=[
            ("VERSION_INFO", __version__),
        ],
        include_dirs=[
            os.path.join("src", "hamiltonian_core", "include"),
            np.get_include(),
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)