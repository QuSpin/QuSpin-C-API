# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
from quspin_core_abi cimport *

__all__ = ["operator_string"]

cdef class operator_string:
    cdef object _name
    def __init__(self):
        self._name = "operator_string"

    @property
    def name(self):
        return self._name