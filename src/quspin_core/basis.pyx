# cython: language_level=3
# distutils: language=c++
from basis cimport *

cdef class bit_basis:
    cdef object _name
    def __init__(self):
        self._name = "basis"