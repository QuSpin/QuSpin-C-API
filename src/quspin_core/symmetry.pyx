# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
###########################
from quspin_core_abi cimport symmetry_abi,bit_perm_args,perm_bit_args,dit_perm_args,perm_dit_args
from numpy cimport complex128_t
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr,make_shared,reinterpret_pointer_cast
###########################
import numpy as _np

def _to_int_array(a):
    return _np.asarray(list(a)).astype(int,casting='equiv')

def _to_bool_array(a):
    return _np.asarray(list(a)).astype(bool,casting='equiv')

def _to_complex_array(a):
    return _np.asarray(list(a)).astype(_np.complex128,casting='safe')

def _check_lat_args(args,chars):
    args = _to_int_array(args) # cast to integers
    args = [[int(ele) for ele in perm] for perm in args] # extract as list of lists
    
    if len(args) != len(chars): 
        raise ValueError("number of lattice symmetries must equal the number of characters.")
    
    return args, _to_complex_array(chars)

def _check_bit_lat_args(args,chars):
    return _check_lat_args(args,chars)

def _check_dit_lat_args(lhss,args,chars):
    return _check_lat_args(args,chars)

def _check_bit_loc_args(args,chars):
    args = _to_bool_array(args) # cast to integers
    args = [[bool(ele) for ele in perm] for perm in args] # extract as list of lists
    
    if len(args) != len(chars): 
        raise ValueError("number of local symmetries must equal the number of characters.")
    
    return args, _to_complex_array(chars)

def _check_dit_loc_args(lhss,args,chars):
    new_args = []
    for (perms,locs) in args:
        perms = _to_int_array(perms)
        locs = _to_int_array(locs)
        
        if perms.ndim != 2: raise ValueError('local dit permutation must be given for each location as a 2d array-like object')
        if locs.ndim != 1: raise ValueError('locations for dit permutations must be given as a 1d array-like object')
        if perms.shape[0] != locs.shape[0]: raise ValueError('the number of locations must match the number of permutations.')
        if perms.shape[1] != lhss: raise ValueError
        perms = [[int(v) for v in perm] for perm in perms[:]]
        locs = [int(v) for v in locs]
        
        new_args.append((perms,locs))
        
    if len(new_args) != len(chars): 
        raise ValueError("number of local symmetries must equal the number of characters.")            
    
    return new_args,_to_complex_array(chars)

def _check_bit_args(lat_list,lat_chars,loc_list,loc_chars):
    return (_check_bit_lat_args(lat_list,lat_chars) +
            _check_bit_loc_args(loc_list,loc_chars))

def _check_dit_args(lhss,lat_list,lat_chars,loc_list,loc_chars):
    return (_check_dit_lat_args(lhss,lat_list,lat_chars) +
            _check_dit_loc_args(lhss,loc_list,loc_chars))

def check_args(lhss,lat_args,loc_args):
    lat_list,lat_chars = list(zip(*lat_args.items()))
    loc_list,loc_chars = list(zip(*loc_args.items()))
    
    if lhss < 2:
        raise ValueError('expecting 1 < lhss < 256.')
    elif lhss == 2:
        return _check_bit_args(lat_list,lat_chars,loc_list,loc_chars)
    else:
        return _check_dit_args(lhss,lat_list,lat_chars,loc_list,loc_chars)



cdef class symmetry_api:
    """An Extension class that constructs a low-level QuSpin symmetry_abi object."""

    cdef symmetry_abi * symm

    def __cinit__(
        self,
        int lhss, 
        int bits,
        dict lat_args = {},
        dict loc_args = {}):

        (
            lat_list, lat_chars,
            loc_list, loc_chars,
        ) = check_args(lhss, lat_args, loc_args)

        cdef vector[shared_ptr[void]] _lat_list
        cdef vector[shared_ptr[void]] _loc_list

        cdef vector[int] int_vec
        cdef vector[vector[int]] int_vec_vec

        if lhss == 2:
            for perm in lat_list:
                int_vec = perm
                _lat_list.push_back(
                    reinterpret_pointer_cast[void,bit_perm_args](
                        make_shared[bit_perm_args](int_vec)
                    )
                )
            
            for mask in loc_list:
                int_vec = mask
                _loc_list.push_back(
                    reinterpret_pointer_cast[void,perm_bit_args](
                        make_shared[perm_bit_args](int_vec)
                    )
                )

        elif lhss > 2:
            for perm in lat_list:
                int_vec = perm
                _lat_list.push_back(
                    reinterpret_pointer_cast[void,dit_perm_args](
                        make_shared[dit_perm_args](int_vec)
                    )
                )

            for perms,locs in loc_list:
                int_vec = locs
                int_vec_vec = perms
                _loc_list.push_back(
                    reinterpret_pointer_cast[void,perm_dit_args](
                        make_shared[perm_dit_args](int_vec_vec,int_vec)
                    )
                )

        cdef complex128_t[::1] _lat_chars = lat_chars
        cdef complex128_t[::1] _loc_chars = loc_chars

        self.symm = new symmetry_abi(lhss,bits,_lat_list,&_lat_chars[0],_loc_list,&_loc_chars[0])

    def __dealloc__(self):
        del self.symm
