from setuptools import find_packages, setup, Extension
from distutils.command import build_ext, install
from Cython.Build import cythonize
import numpy as np
import os,subprocess,sys,warnings

def check_for_boost_includes(include_dirs):
    for include_dir in include_dirs:
        if os.path.exists(os.path.join(include_dir,'boost')):
            return True,include_dirs
    
    # check of boost is installed in local workspace
    for root, dirs, files in os.walk(".", topdown=True):
        if 'boost' in dirs:
            include_dirs.append(root)
            return True,include_dirs
            
    return False,include_dirs

def get_include_dirs():
    from sysconfig import get_paths
    data_path = get_paths()["data"]
    
    include_dirs = [np.get_include(),os.path.join('src','quspin_core','includes')]
    
    if sys.platform == 'win32':
        include_dirs.append(os.path.join(data_path,'Library','include'))
    else:
         include_dirs.append(os.path.join(data_path,'include'))
         
    return include_dirs

def get_extension_kwargs(include_dirs,coverage):
    if sys.platform == 'win32':
        if coverage:
             warnings.warn(UserWarning('coverage not supported on OS-X. ignoring coverage request.'))
        return dict(
            extra_compile_args =  ['/std:c++20'],
            extra_link_args = [],
            include_dirs = include_dirs
        )    
    elif sys.platform == 'darwin':
        if coverage:
             warnings.warn(UserWarning('coverage not supported on OS-X. ignoring coverage request.'))
             
        return dict(
            extra_compile_args = ['--std=c++20'],
            extra_link_args = [],
            include_dirs = include_dirs
        )        
    elif sys.platform == 'linux':
        extra_compile_args = ['--std=c++20']
        extra_link_args = []
        if coverage:
            extra_compile_args += [
                '--coverage', 
                '-fno-inline', 
                '-fno-inline-small-functions', 
                '-fno-default-inline', 
                '-O0'
            ]
            extra_link_args += ['--coverage']
            
            
        return dict(
            extra_compile_args = extra_compile_args,
            extra_link_args = extra_link_args,
            include_dirs = include_dirs
        )
 
 
include_dirs = get_include_dirs()
    
if os.environ.get('BOOST_INCLUDE_DIR',False):
    include_dirs.append(os.environ.get('BOOST_INCLUDE_DIR'))

use_boost,include_dirs = check_for_boost_includes(include_dirs)
extension_kwargs = get_extension_kwargs(
    include_dirs,
    os.environ.get('ENABLE_COVERAGE',False)
)

if not os.environ.get("NO_GENERATE_ABI",False):
    subprocess.check_call([sys.executable,
                            os.path.join(os.path.dirname(__file__),
                                'generate_abi.py'),f'{use_boost}'])

with open('README.md', 'r') as f:
    long_description = f.read()

exec(open(os.path.join('src','quspin_core','_version.py')).read())

ext = [
    Extension('quspin_core.basis', [os.path.join('src','quspin_core','basis.pyx')],
        **extension_kwargs
    ),
    Extension('quspin_core.symmetry', [os.path.join('src','quspin_core','symmetry.pyx')],
        **extension_kwargs
    ),
    Extension('quspin_core.operator', [os.path.join('src','quspin_core','operator.pyx')],
        **extension_kwargs
    ),
]

setup(
    name='quspin-core',
    version=__version__,
    zip_safe=False,
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    author='Phillip Weinberg, Marin Bukov, Markus Schmitt',
    description='Base low-level components for QuSpin.',
    long_description=long_description,
    url='https://github.com/weinbe58/QuSpin-Core',
    ext_modules=cythonize(ext,
        include_path=extension_kwargs['include_dirs']
    ),
    install_requires=[
        'numpy>=1.19.2',
    ],
    test_suite="test_python"
)