from setuptools import find_packages, setup, Extension
from distutils.command import build_ext, install
from Cython.Build import cythonize
import numpy as np
import os

def get_extension_kwargs():
    return dict(
        extra_compile_args = "",
        extra_link_args = "",
        include_dirs = [np.get_include(),os.path.join("src","quspin_core","includes")]
    )
    
extension_kwargs = get_extension_kwargs()

with open("README.md", 'r') as f:
    long_description = f.read()

exec(open(os.path.join("src","quspin_core","_version.py")).read())

ext = [
    Extension("quspin_core.basis", [os.path.join("src","quspin_core","basis.pyx")],
        **extension_kwargs
    ),
    Extension("quspin_core.operator", [os.path.join("src","quspin_core","operator.pyx")],
        **extension_kwargs
    ),
]

setup(
    name="quspin-core",
    version=__version__,
    zip_safe=False,
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    author="Phillip Weinberg, Marin Bukov, Markus Schmitt",
    description="Base low-level components for QuSpin.",
    long_description=long_description,
    url="https://github.com/weinbe58/QuSpin-Core",
    ext_modules=cythonize(ext,
        include_path=extension_kwargs["include_dirs"]
    ),
    install_requires=[
        "numpy>=1.19.2",
    ]
)