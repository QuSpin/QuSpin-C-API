# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
from quspin_core_abi cimport *

__all__ = ["bit_basis"]

cdef class bit_basis:
    cdef object _name
    def __init__(self):
        self._name = "bit_basis"

    @property
    def name(self):
        return self._name