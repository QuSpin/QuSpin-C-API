# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
###########################
import numpy as np


_allowed_types = {
    'int8':np.int8,
    'int16':np.int16,
    'float32':np.float32,
    'float64':np.float64,
    'complex64':np.complex64,
    'complex128':np.complex128
}

def check_is_perm(perm):
    elements = set(perm)
    if len(elements) != len(perm): 
        raise ValueError('list is not valid permutation, list has repeating elements')

    if set(elements) != set(range(len(perm))):
        raise ValueError('list is not valid permutation, list is missing elements.')


    