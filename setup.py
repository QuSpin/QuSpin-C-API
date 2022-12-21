from setuptools import find_packages, setup, Extension
from Cython.Build import cythonize
import numpy as np
import os,sys

def get_includes():
    [np.get_include(),os.path.join(os.path.abspath(__file__),"src","includes")]
    

if __name__ == "__main__":

    with open("README.md", 'r') as f:
        long_description = f.read()
        
    ext = [
        Extension("quspin_api.basis", ["src/quspin_api/basis.pyx"],
            include_dirs=get_includes(),
            libraries=["quspin_abi"],
            library_dirs=["src/impl"]
        ),
        Extension("quspin_api.operators", ["src/quspin_api/operators.pyx"],
            include_dirs=get_includes(),
            libraries=["quspin_abi"],
            library_dirs=["src/impl"],
            language="c++",
        ),
    ]
    
    setup(
        name="QuSpin-Core",
        version="0.0.1",
        zip_safe=False,
        packages=find_packages(path="src/quspin_api"),
        author="Phillip Weinberg, Marin Bukov, Markus Schmitt",
        description="A classification library using a novel audio-inspired algorithm.",
        long_description=long_description,
        long_description_content_type='text/markdown',
        url="https://github.com/lol-cubes/classification-library",
        ext_modules=cythonize(ext),
                install_requires=[
            'numpy>=1.19.2'
        ]
    )