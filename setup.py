# Available at setup time due to pyproject.toml
import glob
import subprocess
import os
import sys
from typing import Dict, List
from pybind11.setup_helpers import Pybind11Extension, build_ext, ParallelCompile
from setuptools import setup

# Optional multithreaded build
ParallelCompile("NPY_NUM_BUILD_JOBS").install()


__version__ = "0.1.0"

CWD = os.path.dirname(os.path.abspath(__file__))
LIBQUISPIN_BUILD_DIR = os.path.join(CWD, "libquspin-build")
LIBQUSPIN_DIR = os.path.join(CWD, "libquspin")


def run_cmd(cmds: list[str]):
    res = subprocess.run(cmds, stdout=sys.stdout, stderr=sys.stderr, cwd=CWD)

    if res.returncode == 0:
        return

    raise RuntimeError("Failed to build libquspin")


class quspin_build_ext(build_ext):

    def run(self):
        run_cmd(
            [
                "meson",
                "setup",
                "libquspin",
                LIBQUISPIN_BUILD_DIR,
                "--reconfigure",
                "--buildtype=release",
            ]
        )
        run_cmd(["meson", "test", "-C", LIBQUISPIN_BUILD_DIR, "-j4"])
        super().run()


def extra_compile_args() -> List[str]:
    if sys.platform == "win32":
        extra_compile_args = ["/std:c++20"]
    elif sys.platform in ["darwin"]:
        extra_compile_args = [
            "--std=c++20",
            "-mmacosx-version-min=10.14",
        ]
    else:
        extra_compile_args = ["-std=c++2a"]

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
        extra_link_args = []
    elif sys.platform in ["darwin"]:
        extra_link_args = []
    else:
        extra_link_args = []

    if os.environ.get("COVERAGE", False):
        if sys.platform in ["win32", "cygwin", "win64", "darwin"]:
            raise ValueError("Coverage is not supported on Windows or macOS")

        extra_link_args += ["--coverage"]

    return extra_link_args


def setup_quspin_core() -> Dict[str, List[str]]:
    if sys.platform == "win32":
        obj_ext = "obj"
        lib_file = "libquspin.lib"
    elif sys.platform == "darwin":
        obj_ext = "o"
        lib_file = "libquspin.dylib"
    elif sys.platform == "linux":
        obj_ext = "o"
        lib_file = "libquspin.so"
    else:
        raise ValueError(f"Unsupported platform {sys.platform}")

    source_names = []

    for root, _, files in os.walk(os.path.join(LIBQUSPIN_DIR, "src")):
        root = os.path.relpath(root, LIBQUSPIN_DIR)
        files = [f for f in files if f.endswith(".cpp")]
        for file in files:
            source_names.append((*root.split(os.path.sep), file))

    extra_objects = []
    for source_name in source_names:
        object_name = "_".join(source_name) + f".{obj_ext}"
        object_path = os.path.join(LIBQUISPIN_BUILD_DIR, f"{lib_file}.p", object_name)
        print(object_path)
        extra_objects.append(object_path)

    if len(extra_objects) == 0:
        object_folder = os.path.join(LIBQUISPIN_BUILD_DIR, f"{lib_file}.p")
        os.listdir(object_folder)
        raise RuntimeError(f"No object files found in {object_folder}")

    include_dirs = [os.path.join(LIBQUSPIN_DIR, "include")]

    bindings_cpp = os.path.join("src", "quspin_core", "_bindings.cpp")
    submodules_cpp = glob.glob(
        os.path.join("src", "quspin_core", "_submodules", "*.cpp")
    )

    include_dirs.append(os.path.join("src", "quspin_core"))

    return Pybind11Extension(
        "quspin_core._bindings",
        sorted([bindings_cpp, *submodules_cpp]),
        extra_objects=extra_objects,
        extra_compile_args=extra_compile_args(),
        extra_link_args=extra_link_args(),
        include_dirs=include_dirs,
        define_macros=[("VERSION_INFO", __version__)],
    )


setup(
    ext_modules=[setup_quspin_core()],
    cmdclass={"build_ext": quspin_build_ext},
    zip_safe=False,
)
