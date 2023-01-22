# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool
from numpy cimport complex128_t, NPY_TYPES

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

cdef extern from 'quspin_core_abi/operator_abi.h' namespace 'quspin_core_abi':

    cdef enum OPERATOR_TYPES:
        OP_STRING
        OP_TWO_BODY

    ctypedef struct operator_string_args:
        const int nlocs
        int * locs
        int * perms
        void * datas

    ctypedef struct N_body_op_args:
        int * locs
        void * data

    cdef cppclass operator_abi:
        operator_abi(
            NPY_TYPES,
            OPERATOR_TYPES,
            const int,
            vector[shared_ptr[void]]) except +




cdef extern from 'quspin_core_abi/numpy_interface.h' namespace 'quspin_core_abi':
    NPY_TYPES npy_typenum(ndarray)
    void * npy_data(ndarray)


