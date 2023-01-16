# cython: language_level=3
# distutils: language=c++
###########################
from quspin_core_abi cimport symmetry_abi,bit_perm_args,perm_bit_args,dit_perm_args,perm_dit_args
from numpy cimport complex128_t
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr,make_shared,reinterpret_pointer_cast
###########################
from ._symmetry_utils import check_args



cdef class symmetry_api:
    """An Extension class that constructs a low-level QuSpin symmetry_abi object."""

    cdef symmetry_abi * symm

    def __cinit__(
        self,
        int lhss, 
        int bits,
        lat_args = [],
        lat_chars = [],
        loc_args = [],
        loc_chars = []):

        (
            lat_args, lat_chars,
            loc_args, loc_chars,
        ) = check_args(
            lhss, lat_args, lat_chars,
            loc_args, loc_chars
        )

        cdef vector[shared_ptr[void]] _lat_args
        cdef vector[shared_ptr[void]] _loc_args

        cdef vector[int] int_vec
        cdef vector[vector[int]] int_vec_vec

        if lhss == 2:
            for perm in lat_args:
                int_vec = perm
                _lat_args.push_back(
                    reinterpret_pointer_cast[void,bit_perm_args](
                        make_shared[bit_perm_args](int_vec)
                    )
                )
            
            for mask in loc_args:
                int_vec = mask
                _loc_args.push_back(
                    reinterpret_pointer_cast[void,perm_bit_args](
                        make_shared[perm_bit_args](int_vec)
                    )
                )

        elif lhss > 2:
            for perm in lat_args:
                int_vec = perm
                _lat_args.push_back(
                    reinterpret_pointer_cast[void,dit_perm_args](
                        make_shared[dit_perm_args](int_vec)
                    )
                )

            for perms,locs in loc_args:
                int_vec = locs
                int_vec_vec = perms
                _loc_args.push_back(
                    reinterpret_pointer_cast[void,perm_dit_args](
                        make_shared[perm_dit_args](int_vec_vec,int_vec)
                    )
                )

        cdef complex128_t[::1] _lat_chars = lat_chars
        cdef complex128_t[::1] _loc_chars = loc_chars

        self.symm = new symmetry_abi(lhss,bits,_lat_args,&_lat_chars[0],_loc_args,&_loc_chars[0])

    def __dealloc__(self):
        del self.symm
