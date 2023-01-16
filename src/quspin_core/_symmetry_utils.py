import numpy as _np

def to_int_array(a):
    return _np.asarray(a).astype(int,casting='equiv')

def to_bool_array(a):
    return _np.asarray(a).astype(bool,casting='equiv')

def to_complex_array(a):
    return _np.asarray(a).astype(_np.complex128,casting='equiv')

def check_lat_args(args,chars):
    args = to_int_array(args) # cast to integers
    args = [[int(ele) for ele in perm] for perm in args] # extract as list of lists
    chars = _np.asarray(chars, dtype=_np.complex128)
    
    if len(args) != len(chars): 
        raise ValueError("number of lattice symmetries must equal the number of characters.")
    
    return args, chars

def check_bit_lat_args(args,chars):
    return check_lat_args(args,chars)

def check_dit_lat_args(lhss,args,chars):
    return check_lat_args(args,chars)

def check_bit_loc_args(args,chars):
    args = to_bool_array(args) # cast to integers
    args = [[bool(ele) for ele in perm] for perm in args] # extract as list of lists
    chars = _np.asarray(chars, dtype=_np.complex128)
    
    if len(args) != len(chars): 
        raise ValueError("number of local symmetries must equal the number of characters.")
    
    return args, chars

def check_dit_loc_args(lhss,args,chars):
    new_args = []
    for (perms,locs) in args:
        perms = to_int_array(perms)
        locs = to_int_array(locs)
        
        if perms.ndim != 2: raise ValueError('local dit permutation must be given for each location as a 2d array-like object')
        if locs.ndim != 1: raise ValueError('locations for dit permutations must be given as a 1d array-like object')
        if perms.shape[0] != locs.shape[0]: raise ValueError('the number of locations must match the number of permutations.')
        if perms.shape[1] != lhss: raise ValueError
        perms = [[int(v) for v in perm] for perm in perms[:]]
        locs = [int(v) for v in locs]
        
        new_args.append((perms,locs))
        
    if len(new_args) != len(chars): 
        raise ValueError("number of local symmetries must equal the number of characters.")            
    
    return new_args,_np.asarray(chars, dtype=_np.complex128)

def check_bit_args(lat_args,lat_chars,loc_args,loc_chars):
    return (check_bit_lat_args(lat_args,lat_chars) +
            check_bit_loc_args(loc_args,loc_chars))

def check_dit_args(lhss,lat_args,lat_chars,loc_args,loc_chars):
    return (check_dit_lat_args(lhss,lat_args,lat_chars) +
            check_dit_loc_args(lhss,loc_args,loc_chars))

def check_args(lhss,lat_args,lat_chars,loc_args,loc_chars):

    if lhss < 2:
        raise ValueError('expecting 1 < lhss < 256.')
    elif lhss == 2:
        return check_bit_args(lat_args,lat_chars,loc_args,loc_chars)
    else:
        return check_dit_args(lhss,lat_args,lat_chars,loc_args,loc_chars)
        
    
    
    
    return (
        lat_args,
        _np.asarray(lat_chars, dtype=_np.complex128),
        loc_args,
        _np.asarray(loc_chars, dtype=_np.complex128)
    )
