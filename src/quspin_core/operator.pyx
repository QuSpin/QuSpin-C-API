# cython: language_level=3
# distutils: language=c++
from operator cimport *

cdef class operator:
    cdef object _name
    def __init__(self):
        self._name = "operator"

    @property
    def name(self):
        return self._name