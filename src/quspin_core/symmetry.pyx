# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
###########################
from quspin_core_abi cimport symmetry_abi,bit_perm_args,perm_bit_args,dit_perm_args,perm_dit_args
from numpy cimport complex128_t, ndarray, PyArray_DATA
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr,make_shared,reinterpret_pointer_cast
###########################
import numpy as np
from ._utils import check_is_perm

cdef class BitPerm:
    cdef vector[int] perm

    def __init__(self,perm):
        perm = [int(ele) for ele in perm]
        check_is_perm(perm)
        self.perm = perm

    @property
    def lhss(self):
        return 2

    def __hash__(self):
        return hash((BitPerm,2,tuple(self.perm)))

    cdef shared_ptr[void] get_arg(self):
        return reinterpret_pointer_cast[void,bit_perm_args](
            make_shared[bit_perm_args](self.perm)
        )

cdef class PermBit:
    cdef vector[int] mask

    def __init__(self,mask):
        mask = [int(ele) for ele in mask]
        if any((ele > 1 or ele < 0) for ele in mask):
            raise ValueError('mask must be collection of bools or {0,1}.')

        self.mask = mask

    @property
    def lhss(self):
        return 2

    def __hash__(self):
        return hash((PermBit,2,tuple(self.mask)))

    cdef shared_ptr[void] get_arg(self):
        return reinterpret_pointer_cast[void,perm_bit_args](
            make_shared[perm_bit_args](self.mask)
        )

cdef class DitPerm:
    cdef vector[int] perm
    cdef int lhss

    def __init__(self,lhss,perm):
        perm = [int(ele) for ele in perm]
        check_is_perm(perm)
        self.perm = perm
        self.lhss = lhss

    @property
    def lhss(self):
        return self.lhss

    def __hash__(self):
        return hash((DitPerm,self.lhss,tuple(self.perm)))

    cdef shared_ptr[void] get_arg(self):
        return reinterpret_pointer_cast[void,dit_perm_args](
            make_shared[dit_perm_args](self.perm)
        )

cdef class PermDit:
    cdef vector[vector[int]] perms
    cdef vector[int] locs
    cdef int lhss

    def __init__(self,lhss,locs,perms):
        locs = [int(ele) for ele in locs]
        perms = [[int(ele) for ele in perm] for perm in perms]

        if len(locs) != len(perms):
            raise ValueError('number of elements in locs must equal the number of permutations in perms.')

        if any(len(perm) != lhss for perm in perms):
            raise ValueError('length of permutations must be equal to "lhss"')

        cdef vector[int] perm_vec

        self.lhss = lhss
        self.locs = locs
        for perm in perms:
            check_is_perm(perm)
            perm_vec = perm
            self.perms.push_back(perm_vec)

    @property
    def lhss(self):
        return self.lhss

    def __hash__(self):
        perms = tuple(tuple( ele for ele in perm) for perm in self.perms)
        locs = tuple(self.locs)
        return hash((PermDit,self.lhss,locs,perms))

    cdef shared_ptr[void] get_arg(self):
        return reinterpret_pointer_cast[void,perm_dit_args](
            make_shared[perm_dit_args](self.perms,self.locs)
        )
 
cdef class SymmetryAPI:
    """An Extension class that constructs a low-level QuSpin symmetry_abi object."""

    cdef symmetry_abi * symm

    def __cinit__(
        self,
        int lhss, 
        int bits,
        dict lat_args = {},
        dict loc_args = {}):

        if len(lat_args) > 0:
            lat_symm,lat_chars = list(zip(*lat_args.items()))
        else:
            lat_symm,lat_chars = [],[]
        
        if len(loc_args) > 0:
            loc_symm,loc_chars = list(zip(*loc_args.items()))
        else:
            loc_symm,loc_chars = [],[]

        for symm in lat_symm:

            if (type(symm) not in [BitPerm,DitPerm]) or (type(symm) != type(lat_symm[0])):
                raise TypeError(
                    'all symmetry objects in "lat_args" must be one '
                    'of the same type: "BitPerm" or "DitPerm"'
                )

            if symm.lhss != lhss:
                raise ValueError(
                    'all symmetry in "lat_args" must have '
                    'the same local hilbert space size.'
                    f'found {symm.lhss}, expecting {lhss}'
                )

        for symm in loc_symm:
            if (type(symm) not in [PermBit,PermDit]) or (type(symm) != type(loc_symm[0])):
                raise TypeError(
                    'all symmetry objects in "loc_args" must be one '
                    'of the same type: "PermBit" or "PermDit"'
                )

            if symm.lhss != lhss:
                raise ValueError(
                    'all symmetry in "loc_args" must have '
                    'the same local hilbert space size.'
                    f'found {symm.lhss}, expecting {lhss}'
                )

        cdef vector[shared_ptr[void]] lat_args_vec
        cdef vector[shared_ptr[void]] loc_args_vec

        cdef PermBit perm_bit
        cdef BitPerm bit_perm
        cdef PermDit perm_dit
        cdef DitPerm dit_perm

        if len(lat_args) > 0 and len(loc_args) > 0:
            if (type(lat_symm[0]) == BitPerm and type(loc_symm[0]) == PermBit):
                for bit_perm in lat_symm:
                    lat_args_vec.push_back(bit_perm.get_arg())

                for perm_bit in loc_symm:
                    loc_args_vec.push_back(perm_bit.get_arg())

            elif (type(lat_symm[0]) == DitPerm and type(loc_symm[0]) == PermDit):
                for dit_perm in lat_symm:
                    lat_args_vec.push_back(dit_perm.get_arg())

                for perm_dit in loc_symm:
                    loc_args_vec.push_back(perm_dit.get_arg())
            else:
                raise TypeError('cannot mix dit and bit symmetries.')

        cdef ndarray _lat_chars = np.ascontiguousarray(lat_chars,dtype=np.complex128)
        cdef ndarray _loc_chars = np.ascontiguousarray(loc_chars,dtype=np.complex128)
        
        self.symm = new symmetry_abi(
            lhss,bits,
            lat_args_vec,
            <complex128_t*>PyArray_DATA(_lat_chars),
            loc_args_vec,
            <complex128_t*>PyArray_DATA(_loc_chars)
        )

    def __dealloc__(self):
        del self.symm
