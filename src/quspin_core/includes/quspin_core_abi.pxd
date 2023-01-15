from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool
from numpy cimport complex128_t

cdef extern from 'quspin_core_abi/quspin_core_abi.h' namespace 'quspin_core_abi':
    pass

cdef extern from 'quspin_core_abi/symmetry_abi.h' namespace 'quspin_core_abi':
    ctypedef struct bit_perm_args:
        vector[int] perm

    ctypedef struct dit_perm_args:
        vector[int] perm

    ctypedef struct perm_bit_args:
        vector[int] mask

    ctypedef struct perm_dit_args:
        vector[vector[int]] perm
        vector[int] locs

    cdef cppclass symmetry_abi:
        symmetry_abi(
            const int,
            const size_t,
            vector[shared_ptr[void]],
            const complex128_t *,
            vector[shared_ptr[void]],
            const complex128_t *)  except +



