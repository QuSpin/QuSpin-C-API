# cython: language_level=3
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++
from quspin_core_abi cimport *
from numpy cimport ndarray,PyArray_DATA,PyArray_TYPE
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr,make_shared,reinterpret_pointer_cast
############################################################################
from ._utils import check_is_perm, _allowed_types
import numpy as np



cdef class OperatorString:
    cdef ndarray datas
    cdef ndarray perms
    cdef ndarray locs
    cdef int nlocs
    cdef int lhss

    def __init__(self,locs, perms, datas):
        self.locs = np.ascontiguousarray(locs,dtype=np.intc)
        self.datas = np.atleast_2d(np.ascontiguousarray(datas))
        self.perms = np.atleast_2d(np.ascontiguousarray(perms,dtype=np.intc))

        if np.dtype(self.datas.dtype).name not in _allowed_types:
            types = ', '.join(_allowed_types.keys())
            raise TypeError(f'dtype must be one of: {types}')

        if self.perms.ndim > 2:
            raise ValueError('"perms" must be a 2d array of integers.')

        if self.datas.ndim > 2:
            raise ValueError('"datas" must be a 2d array of integers.')

        if (self.perms.shape[0] != self.datas.shape[0] or self.perms.shape[1] != self.datas.shape[1]):
            raise ValueError('"perms" and "datas" must be the same shape.')

        if self.locs.ndim > 1:
            raise ValueError('"locs" must be the an integer or a 1d array of integers.')

        if self.locs.size != self.perms.shape[0]:
            raise ValueError('number of rows in "perms" must equal the number of elements in "locs"')

        for perm in self.perms[:]:
            check_is_perm(perm)

        self.nlocs = self.locs.size
        self.lhss = self.perms.shape[1]

    @property
    def datas(self):
        arr = self.datas
        arr.flags.writeable=False
        return arr

    @property
    def perms(self):
        arr = self.perms
        arr.flags.writeable=False
        return arr

    @property
    def locs(self):
        arr = self.locs
        arr.flags.writeable=False
        return arr

    @property
    def lhss(self):
        return self.lhss

    def astype(self,dtype,**kwargs):    
        return OperatorString(self.locs, self.perms, self.datas.astype(dtype,**kwargs))

    cdef OPERATOR_TYPES get_term_type(self):
        return OP_STRING

    cdef NPY_TYPES get_typenum(self):
        return <NPY_TYPES>PyArray_TYPE(self.datas)

    cdef shared_ptr[void] get_arg(self):
        return reinterpret_pointer_cast[void,operator_string_args](make_shared[operator_string_args](
                self.nlocs, 
                <int*>PyArray_DATA(self.locs), 
                <int*>PyArray_DATA(self.perms), 
                <void*>PyArray_DATA(self.datas)
            )
        )

cdef class NBodyOperator:
    cdef ndarray data
    cdef ndarray locs
    cdef int lhss

    def __init__(self,locs,data):
        self.locs = np.ascontiguousarray(locs,dtype=np.intc)
        self.data = np.atleast_2d(np.ascontiguousarray(data))

        if np.dtype(self.data.dtype).name not in _allowed_types:
            types = ', '.join(_allowed_types.keys())
            raise TypeError(f'dtype must be one of: {types}')

        if self.data.ndim > 2 or self.data.shape[0] != self.data.shape[1]:
            raise ValueError('"data" must be square array. ')

        if self.locs.ndim > 1:
            raise ValueError('"locs" must be a 1d array of integers')

        if len(self.locs) not in  [2]:
            raise ValueError('NBodyOperator only supports two-body terms in the QuSpin-ABI')

        n_body = len(self.locs)

        self.lhss = int((self.data.shape[0])**(1.0/n_body))

        if self.lhss**n_body !=self.data.shape[0]:
            raise ValueError(f'number of rows in "data" must must be the {n_body}th power of an integer, expecting {self.lhss**n_body}, got {self.data.shape[0]}.')

    def astype(self,dtype,**kwargs):
        return NBodyOperator(self.locs,self.data.astype(dtype,**kwargs))


    @property
    def data(self):
        arr = self.data
        arr.flags.writeable=False
        return arr

    @property
    def locs(self):
        arr = self.locs
        arr.flags.writeable=False
        return arr

    @property
    def lhss(self):
        return self.lhss

    cdef OPERATOR_TYPES get_term_type(self):
        return OP_TWO_BODY
            
    
    cdef NPY_TYPES get_typenum(self):
        return <NPY_TYPES>PyArray_TYPE(self.data)

    cdef shared_ptr[void] get_arg(self):

        return reinterpret_pointer_cast[void,N_body_op_args](make_shared[N_body_op_args](
                <int*>PyArray_DATA(self.locs),
                <void*>PyArray_DATA(self.data)
            )
        )


cdef class OperatorAPI:

    cdef operator_abi * terms

    def __cinit__(self,terms,dtype):


        terms = list(terms)

        if len(terms) == 0:
            raise ValueError("OperatorAPI requires at least one Term")

        for term in terms:

            if (type(term) not in [OperatorString,NBodyOperator]) or (type(term) != type(terms[0])):
                raise ValueError(
                    'all operators in "terms" must be one '
                    'of the same type: "OperatorString" or "NBodyOperator"'
                )

            if term.lhss != terms[0].lhss:
                raise ValueError(
                    'all operators in "terms" must have '
                    'the same local hilbert space size.'
                    f'found {term.lhss}, expecting {terms[0].lhss}'
                )

        terms = [term.astype(dtype,copy=True) for term in terms]

        cdef vector[shared_ptr[void]] args
        cdef NBodyOperator n_body_op
        cdef OperatorString op_string
        cdef int lhss
        cdef NPY_TYPES T_typenum
        cdef OPERATOR_TYPES term_type

        if type(terms[0]) == OperatorString:
            for op_string in terms:
                args.push_back(op_string.get_arg())
            
            op_string = terms[0]
            T_typenum = op_string.get_typenum()
            term_type = op_string.get_term_type()
            lhss = op_string.lhss
        else:
            for n_body_op in terms:
                args.push_back(n_body_op.get_arg())

            n_body_op = terms[0]
            T_typenum = n_body_op.get_typenum()
            term_type = n_body_op.get_term_type()
            lhss = n_body_op.lhss


        self.terms = new operator_abi(T_typenum,term_type,lhss,args)

    


        
    
        
        

