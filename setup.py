from setuptools import find_packages, setup, Extension
from Cython.Build import cythonize
import numpy as np
import os


def get_includes():
    return [np.get_include(),os.path.join("src","quspin_core","includes")]

if __name__ == "__main__":

    with open("README.md", 'r') as f:
        long_description = f.read()
    
    exec(open(os.path.join("src","quspin_core","_version.py")).read())
    ext = [
        Extension("quspin_core.basis", [os.path.join("src","quspin_core","basis.pyx")],
            include_dirs=get_includes(),
        ),
        Extension("quspin_core.operators", [os.path.join("src","quspin_core","operator.pyx")],
            include_dirs=get_includes(),
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
            include_path=get_includes()
        ),
        install_requires=[
            "numpy>=1.19.2",
        ],
        include_dirs=get_includes(),
    )