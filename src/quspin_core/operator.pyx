# cython: language_level=3
# distutils: language=c++
from operator cimport *

__all__ = ["operator_string"]

cdef class operator_string:
    cdef object _name
    def __init__(self):
        self._name = "operator_string"

    @property
    def name(self):
        return self._name