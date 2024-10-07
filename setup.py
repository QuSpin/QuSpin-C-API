# Available at setup time due to pyproject.toml
import glob
import subprocess
import os
import sys
from typing import List
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.1.0"

LIB_QUISPIN_BUILD = False
LIB_QUISPIN_BUILD_DIR = "libquspin-build"


def extra_compile_args() -> List[str]:
    if sys.platform in ["win32", "cygwin", "win64"]:
        extra_compile_args = ["/openmp", "/std:c++17"]
    if sys.platform in ["darwin"]:
        extra_compile_args = [
            "-DLLVM_ENABLE_PROJECTS",
            "-Xpreprocessor",
            "-fopenmp-version=50" "-fopenmp",
            "--std=c++20",
        ]
    else:
        extra_compile_args = ["-fopenmp", "--std=c++17"]

    if os.environ.get("COVERAGE", False):
        if sys.platform in ["win32", "cygwin", "win64", "darwin"]:
            raise ValueError("Coverage is not supported on Windows or macOS")

        extra_compile_args += [
            "--coverage",
            "-fno-inline",
            "-fno-inline-small-functions",
            "-fno-default-inline",
            "-O0",
        ]

    return extra_compile_args


def extra_link_args() -> List[str]:
    if sys.platform in ["win32", "cygwin", "win64"]:
        extra_link_args = ["/openmp"]
    if sys.platform in ["darwin"]:
        extra_link_args = [
            "-DLLVM_ENABLE_PROJECTS",
            "-Xpreprocessor",
            "-fopenmp-version=50" "-fopenmp",
        ]
    else:
        extra_link_args = ["-fopenmp"]

    if os.environ.get("COVERAGE", False):
        if sys.platform in ["win32", "cygwin", "win64", "darwin"]:
            raise ValueError("Coverage is not supported on Windows or macOS")

        extra_link_args += ["--coverage"]

    return extra_link_args


def build_libquspin():
    if LIB_QUISPIN_BUILD:
        return

    res = subprocess.check_call(
        ["meson", "setup", "libquspin", LIB_QUISPIN_BUILD_DIR, "--reconfigure"]
    )
    if res != 0:
        raise RuntimeError("Failed to build libquspin")

    res = subprocess.check_call(
        ["meson", "compile", "-C", LIB_QUISPIN_BUILD_DIR, "-j", "4"]
    )
    if res != 0:
        raise RuntimeError("Failed to build libquspin")


def generate_extension(name: str) -> Pybind11Extension:
    build_libquspin()

    static_libraries = ["quspin"]
    static_lib_dir = LIB_QUISPIN_BUILD_DIR
    libraries = []
    library_dirs = []
    extra_objects = []

    if sys.platform == "win32":
        libraries.extend(static_libraries)
        library_dirs.append(static_lib_dir)
        extra_objects = []
    else:  # POSIX
        extra_objects = [f"{static_lib_dir}/lib{lib}.a" for lib in static_libraries]

    return Pybind11Extension(
        f"quspin_core.{name}",
        [os.path.join("src", "quspin_core", f"{name}.cpp")],
        define_macros=[
            ("VERSION_INFO", __version__),
        ],
        include_dirs=[
            os.path.join("libquspin", "include"),
        ],
        libraries=libraries,
        library_dirs=library_dirs,
        extra_objects=extra_objects,
        extra_compile_args=extra_compile_args(),
        extra_link_args=extra_link_args(),
    )


def find_extensions() -> List[Pybind11Extension]:
    extensions = []
    for filename in glob.glob("src/quspin_core/*.cpp"):
        name = os.path.splitext(os.path.basename(filename))[0]
        extensions.append(generate_extension(name))

    return extensions


setup(
    ext_modules=find_extensions(),
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
