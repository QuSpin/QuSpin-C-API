import numpy as np
from numpy import int8,int16,int32,int64,float32,float64,complex64,complex128
import os
from typing import Optional
from itertools import product



numpy_ctypes={float32:"float",float64:"double",complex64:"npy_cfloat_wrapper",complex128:"npy_cdouble_wrapper",
                int32:"npy_int32",int64:"npy_int64",int8:"npy_int8",int16:"npy_int16"}

numpy_numtypes={float32:["NPY_FLOAT32"],float64:["NPY_FLOAT64"],complex64:["NPY_COMPLEX64"],complex128:["NPY_COMPLEX128"],
                int32:["NPY_INT32"],int64:["NPY_INT64"],int16:["NPY_INT16","NPY_SHORT"],int8:["NPY_INT8","NPY_BYTE"]}

real_dtypes={float32:float32, float64:float64, complex64:float32, complex128:float64}

index_types = [int32,int64]
matrix_types = [int8,int16,float32,float64,complex64,complex128]
symmetric_matrix_types = [float32,float64,complex64,complex128]

on_the_fly_types = [float32,float64,complex64,complex128]

class cpp:
    
    @staticmethod
    def emit_declare(name:str,type:str) -> str:
        return f'{type} {name}'

    @staticmethod
    def emit_var(name:str,type:str) -> str:
        return cpp.emit_declare(name,type)+';'

    @staticmethod
    def emit_typedef(name:str,type:str) -> str:
        return f'typedef {type} {name};'
    
    @staticmethod
    def emit_using(name:str,type:str,template:Optional[str] = None) -> str:
        if template == None:
            return f'using {name} = {type};'
        else:
            return f'{template}\nusing {name} = {type};'

    @staticmethod
    def emit_case(var_name,cases:dict,default):
        case_list = []
        for (case_value,case_code) in cases.items():
            case_code.replace('\n','\n    ') # indent
            case_list.append(f'    case {case_value}:\n{case_code}'.replace('\n','\n        '))
        case_list.append(f'    default:\n        {default}\n')
        case_body  = '\n'.join(case_list)
        return f'switch({var_name})\n{{\n{case_body}\n}}'            
            
    @staticmethod
    def emit_method(name:str,ret:str,args:list[str],body:str,const_method:bool = False) -> str:
        body = body.replace('\n','\n    ') # indent
        args = ','.join(args)
        if const_method:
            return f'{ret} {name}({args}) const \n{{\n    {body}\n}}'
        else:
            return f'{ret} {name}({args})\n{{\n    {body}\n}}'

    @staticmethod
    def emit_destructor(name:str,body:str="") -> str:
        body = body.replace('\n','\n    ')
        return f'~{name}()\n{{\n    {body}\n}}'

    @staticmethod
    def emit_constructor(name:str,args:list[str] = [], body:str = "",preconstruct : Optional[str] = None, ) -> str:
        
        body = "    "+body.replace('\n','\n    ')
        args = ','.join(args)
        if type(preconstruct) == str:
            return f'{name}({args}) : {preconstruct}\n{{\n{body}\n}}'
        else:
           return  f'{name}({args})\n{{\n{body}\n}}'

    
    @staticmethod
    def emit_class(name:str,constructor:str,destructor:str,public_list:list[str]=[],private_list:list[str]=[]) -> str:
        if private_list:
            private_attr = '    ' + ('\n\n    ').join([b.replace('\n','\n    ') for b in private_list])
        else:
            private_attr = ''
        
        if public_list:
            public_attr = '    ' + ('\n\n    ').join([b.replace('\n','\n    ') for b in public_list])
        else:
            public_attr = ''
        constructor  = '    ' + constructor.replace('\n','\n    ')
        destructor   = '    ' + destructor.replace('\n','\n    ')

        return f'class {name}\n{{\nprivate:\n{private_attr}\npublic:\n\n{constructor}\n\n{destructor}\n\n{public_attr}\n}};'

class basis:
    @staticmethod
    def emit_term_switch_code_generator(int_types:list,boost_types:list,index_types:list,matrix_types:list):
        cases = {}
        switch_code = 0
        switch_code_ctypes = {}

        for basis_ctype,basis_jtype,basis_ktype,bits in int_types:
            bits_case_body = 'if(0){}\n'

            for index_type in index_types:
                index_type_num = numpy_numtypes[index_type][0]
                index_ctype = numpy_ctypes[index_type]
                
                matrix_type_cases = {}
                for matrix_type in matrix_types:
                    matrix_ctype = numpy_ctypes[matrix_type]
                    for matrix_type_num in numpy_numtypes[matrix_type]:
                        matrix_type_cases[matrix_type_num] = f'return {switch_code};'
                    
                    switch_code_ctypes[switch_code] = (bits,basis_ctype,index_ctype,matrix_ctype)
                    switch_code += 1
                
                matrix_type_case = cpp.emit_case("T_typenum",matrix_type_cases,
                    default='return -1;'
                )
                matrix_type_case = matrix_type_case.replace('\n','\n    ')
                
                bits_case_body += f'else if(PyArray_EquivTypenums(J_typenum,{index_type_num}))\n{{\n    {matrix_type_case}\n}}\n'

            cases[bits] = bits_case_body + "else {return -1;}"
        
        method_body = cpp.emit_case("bits",cases,'return -1;')
        args = [
            cpp.emit_declare('bits','const size_t'),
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES')
        ]
        
        return switch_code_ctypes,cpp.emit_method("term_switch_code_generator","static size_t",args,method_body,const_method=False)

    @staticmethod
    def emit_on_the_fly_switch_code_generator(int_types:list,boost_types:list,matrix_types:list):
        cases = {}
        switch_code = 0
        switch_code_ctypes = {}

        iterate_typenums = [(T,T_typenum) for T in matrix_types for T_typenum in numpy_numtypes[T]]

        for basis_ctype,basis_jtype,basis_ktype,bits in int_types:
            T_cases =  {}
            for (T,T_typenum) in iterate_typenums:
                X_cases = {}
                for (X,X_typenum) in iterate_typenums:
                    R = np.result_type(X,T)
                    for Y in matrix_types:
                        if R == Y:
                            Y_conditional = '||'.join([f'Y_typenum == {Y_typenum}' for Y_typenum in numpy_numtypes[Y]])
                            X_cases[X_typenum] = f'return ( {Y_conditional} ? {switch_code} : -1);'
                            switch_code_ctypes[switch_code] = (bits, basis_ctype, numpy_ctypes[T], numpy_ctypes[X], numpy_ctypes[Y])
                            switch_code += 1
                            
                T_cases[T_typenum] = cpp.emit_case("X_typenum",X_cases,"return -1;")
            cases[bits] = cpp.emit_case("T_typenum",T_cases,"return -1;")
        
        method_body = cpp.emit_case("bits",cases,'return -1;')
        args = [
            cpp.emit_declare('bits','const size_t'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('X_typenum','NPY_TYPES'),
            cpp.emit_declare('Y_typenum','NPY_TYPES'),
        ]
        
        return switch_code_ctypes,cpp.emit_method("on_the_fly_switch_code_generator","static size_t",args,method_body,const_method=False)

    @staticmethod
    def emit_on_the_fly(codes:dict, bitbasis:str,operator:str,template:str):
        cases = {}
        
        for code,(bits,I,T,X,Y) in codes.items():
            bit_basis = f'{bitbasis}_{bits}'
            tmp = template.format(T=T)
            cases[code] = (
                f'std::reinterpret_pointer_cast<{bit_basis}>(basis_ptr)->on_the_fly'\
                f'((const {tmp}*)terms, '\
                f'nterms, '\
                f'*(const {Y}*)a, '\
                f'(const {X}*)input, '\
                f'*(const {Y}*)b, '\
                f'({Y}*)output);'
            )
            
        args = [
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('X_typenum','NPY_TYPES'),
            cpp.emit_declare('Y_typenum','NPY_TYPES'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('a','void*'),
            cpp.emit_declare('input','void*'),
            cpp.emit_declare('b','void*'),
            cpp.emit_declare('output','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = on_the_fly_switch_code_generator(bits,T_typenum,X_typenum,Y_typenum);\n{switch}'
        return cpp.emit_method(f'on_the_fly_{operator}','size_t',args,body,const_method=True)

    @staticmethod
    def emit_calc_rowptr(codes:dict, bitbasis:str,operator:str,template:str):
        cases = {}
        
        for code,(bits,I,J,T) in codes.items():
            tmp = template.format(T=T)
            cases[code] = (
                f'std::reinterpret_pointer_cast<{bitbasis}_{bits}>(basis_ptr)->calc_rowptr'\
                f'((const {tmp}*)terms, '\
                f'nterms, '\
                f'({J}*)rowptr);'
            )
        args = [
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('rowptr','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);\n{switch}'
        return cpp.emit_method(f'calc_rowptr_{operator}','size_t',args,body,const_method=True)

    @staticmethod
    def emit_calc_matrix(codes:dict, bitbasis:str,operator:str,template:str):
        cases = {}
        
        for code,(bits,I,J,T) in codes.items():
            tmp = template.format(T=T)
            cases[code] = (
                f'std::reinterpret_pointer_cast<{bitbasis}_{bits}>(basis_ptr)->calc_matrix'\
                f'((const {tmp}*)terms, '\
                f'nterms, '\
                f'({T}*)values, '\
                f'({J}*)indices, '\
                f'({J}*)rowptr);'
            )
        args = [
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('values',"void*"),
            cpp.emit_declare('indices','void*'),
            cpp.emit_declare('rowptr','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);\n{switch}'
        return cpp.emit_method(f'calc_matrix_{operator}','size_t',args,body,const_method=True)

    @staticmethod
    def emit_build_subspace(codes:dict, bitbasis:str,operator:str,template:str):
        cases = {}
        
        for code,(bits,I,J,T) in codes.items():
            tmp = template.format(T=T)
            cases[code] = (
                f'std::reinterpret_pointer_cast<{bitbasis}_{bits}>(basis_ptr)->build_subspace'\
                f'((const {tmp}*)terms, '\
                f'nterms, '\
                f'seed_state, lhss);'
            )
        args = [
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('seed_state','const std::vector<int>&')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);\n{switch}'
        return cpp.emit_method(f'build_subspace_{operator}','size_t',args,body,const_method=False)

    @staticmethod
    def emit_get_state(int_types:list,boost_types:list, bitbasis:str):
        args = [
            cpp.emit_declare('state_index','const size_t'),
            cpp.emit_declare('output', 'std::vector<int>&'),
            cpp.emit_declare('length = 0','const int'),
        ]
        cases = {}
        for  ctype,jtype,ktype,bits in int_types:
            cases[bits] = (
            '{\n'
            f'    auto state = reinterpret_pointer_cast<{bitbasis}_{bits}>(basis_ptr)->space->get_state(state_index);\n'\
            f'    auto state_vec = state.to_vector(length);\n'
            f'    output.insert(output.end(),state_vec.begin(),state_vec.end());\n'
            f'    break;\n'
            '}'
        )
        body = cpp.emit_case('bits',cases,'throw std::runtime_error("unreachanble reached.");')
        return cpp.emit_method("get_state","void",args,body,const_method=True)

        
def emit_symmetric_bitbasis_attr(int_types:list,boost_types:list) -> tuple:
    on_the_fly_switch_codes,on_the_fly_switch_code_generator = basis.emit_on_the_fly_switch_code_generator(int_types,boost_types,symmetric_matrix_types)
    term_switch_codes,term_switch_code_generator = basis.emit_term_switch_code_generator(int_types,boost_types,index_types,symmetric_matrix_types)
        
    
    private_list = [
        cpp.emit_var('basis_ptr','std::shared_ptr<void>'),
        cpp.emit_var('bits','const size_t'),
        cpp.emit_var('lhss = 2','static const int'),
        term_switch_code_generator,
        on_the_fly_switch_code_generator,
    ]
    
    public_list = [
        basis.emit_calc_rowptr(term_switch_codes,'symmetric_bitbasis','operator_string','quspin::operator_string<{T}>'),
        basis.emit_calc_matrix(term_switch_codes,'symmetric_bitbasis','operator_string','quspin::operator_string<{T}>'),        
        basis.emit_on_the_fly(on_the_fly_switch_codes,'symmetric_bitbasis','operator_string','quspin::operator_string<{T}>'),
        basis.emit_build_subspace(term_switch_codes,'symmetric_bitbasis','operator_string','quspin::operator_string<{T}>'),
        # basis.emit_calc_rowptr(term_switch_codes,'symmetric_bitbasis','two_body','quspin::N_body_bit_op<{T},2>'),
        # basis.emit_calc_matrix(term_switch_codes,'symmetric_bitbasis','two_body','quspin::N_body_bit_op<{T},2>'),
        # basis.emit_on_the_fly(on_the_fly_switch_codes,'symmetric_bitbasis','two_body','quspin::N_body_bit_op<{T},2>'),
        # basis.emit_build_subspace(term_switch_codes,'symmetric_bitbasis','two_body','quspin::N_body_bit_op<{T},2>'),
        basis.emit_get_state(int_types,boost_types,'symmetric_bitbasis')
    ]
    
    return private_list,public_list


def emit_symmetric_bitbasis_constructor(int_types:list,boost_types:list) -> str:
    args = [
        cpp.emit_declare('_bits','const size_t'),
        cpp.emit_declare('symmetry','std::shared_ptr<void>'),
        cpp.emit_declare('Ns_est = 0','const size_t'),
    ]
    # body = "if(0){}\n"    
    # for  ctype,jtype,ktype,bits in int_types:
    #     body += (
    #     f'else if(_bits == {bits}){{\n'\
    #     f'    std::shared_ptr<bit_subspace_{bits}> _space = std::make_shared<bit_subspace_{bits}>(Ns_est);\n'\
    #     f'    std::shared_ptr<bit_symmetry<{ctype}>> _symmetry = std::reinterpret_pointer_cast<bit_symmetry<{ctype}>>(symmetry);\n'\
    #     f'    std::shared_ptr<symmetric_bitbasis_{bits}> _basis_ptr = std::make_shared<symmetric_bitbasis_{bits}>(*_symmetry,_space);\n'\
    #     f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'\
    #     f'}}\n'
    # )

    # body += 'else{throw std::runtime_error("number of bits not supported");}\n'
    
    cases = {}
    for  ctype,jtype,ktype,bits in int_types:
        cases[bits] = (
        '{\n'
        f'    std::shared_ptr<bit_subspace_{bits}> _space = std::make_shared<bit_subspace_{bits}>(Ns_est);\n'\
        f'    std::shared_ptr<bit_symmetry<{ctype}>> _symmetry = std::reinterpret_pointer_cast<bit_symmetry<{ctype}>>(symmetry);\n'\
        f'    std::shared_ptr<symmetric_bitbasis_{bits}> _basis_ptr = std::make_shared<symmetric_bitbasis_{bits}>(*_symmetry,_space);\n'\
        f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'\
        f'    break;\n'\
        '}'
        )

    body = cpp.emit_case("_bits",cases,'throw std::runtime_error("number of bits not supported");')

    return cpp.emit_constructor('symmetric_bitbasis_abi',args=args,body=body,
                                preconstruct='bits(_bits)')


def emit_symmetric_bitbasis_class(int_types:list,boost_types:list) -> str:
    name = 'symmetric_bitbasis_abi'
    constructor = emit_symmetric_bitbasis_constructor(int_types,boost_types)
    destructor = cpp.emit_destructor(name)
    private_list,public_list = emit_symmetric_bitbasis_attr(int_types,boost_types)

    return cpp.emit_class(name,constructor,destructor,public_list,private_list)


def emit_symmetric_ditbasis_attr(int_types:list,boost_types:list) -> tuple:
    on_the_fly_switch_codes,on_the_fly_switch_code_generator = basis.emit_on_the_fly_switch_code_generator(int_types,boost_types,symmetric_matrix_types)
    term_switch_codes,term_switch_code_generator = basis.emit_term_switch_code_generator(int_types,boost_types,index_types,symmetric_matrix_types)
        
    
    private_list = [
        cpp.emit_var('basis_ptr','std::shared_ptr<void>'),
        cpp.emit_var('bits','const size_t'),
        cpp.emit_var('lhss ','const int'),
        term_switch_code_generator,
        on_the_fly_switch_code_generator,
    ]
    
    public_list = [
        basis.emit_calc_rowptr(term_switch_codes,'symmetric_ditbasis','operator_string','quspin::operator_string<{T}>'),
        basis.emit_calc_matrix(term_switch_codes,'symmetric_ditbasis','operator_string','quspin::operator_string<{T}>'),
        basis.emit_on_the_fly(on_the_fly_switch_codes,'symmetric_ditbasis','operator_string','quspin::operator_string<{T}>'),
        basis.emit_build_subspace(term_switch_codes,'symmetric_ditbasis','operator_string','quspin::operator_string<{T}>'),
        # basis.emit_calc_rowptr(term_switch_codes,'symmetric_ditbasis','two_body','quspin::N_body_dit_op<{T},2>'),
        # basis.emit_calc_matrix(term_switch_codes,'symmetric_ditbasis','two_body','quspin::N_body_dit_op<{T},2>'),
        # basis.emit_on_the_fly(on_the_fly_switch_codes,'symmetric_ditbasis','two_body','quspin::N_body_dit_op<{T},2>'),
        # basis.emit_build_subspace(term_switch_codes,'symmetric_ditbasis','two_body','quspin::N_body_dit_op<{T},2>'),
        basis.emit_get_state(int_types,boost_types,'symmetric_ditbasis')

    ]
    
    return private_list,public_list


def emit_symmetric_ditbasis_constructor(int_types:list,boost_types:list) -> str:
    args = [
        cpp.emit_declare('_bits','const size_t'),
        cpp.emit_declare('_lhss','const int'),
        cpp.emit_declare('symmetry','std::shared_ptr<void>'),
        cpp.emit_declare('Ns_est = 0','const size_t'),
    ]

    body = "if(0){}\n"    
    for  ctype,jtype,ktype,bits in int_types:
        body += (
        f'else if(_bits == {bits}){{\n'\
        f'    auto _space = std::make_shared<dit_subspace_{bits}>(Ns_est,_lhss);\n'\
        f'    auto _symmetry = std::reinterpret_pointer_cast<dit_symmetry<{ctype}>>(symmetry);\n'\
        f'    auto _basis_ptr = std::make_shared<symmetric_ditbasis_{bits}>(*_symmetry,_space);\n'\
        f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'\
        f'}}\n'
    )

    body +='else{throw std::runtime_error("number of bits not supported");}\n'
    

    return cpp.emit_constructor('symmetric_ditbasis_abi',args=args,body=body,
                                preconstruct='bits(_bits), lhss(_lhss)')


def emit_symmetric_ditbasis_class(int_types:list,boost_types:list) -> str:
    name = 'symmetric_ditbasis_abi'
    constructor = emit_symmetric_ditbasis_constructor(int_types,boost_types)
    destructor = cpp.emit_destructor(name)
    private_list,public_list = emit_symmetric_ditbasis_attr(int_types,boost_types)

    return cpp.emit_class(name,constructor,destructor,public_list,private_list)


def emit_bitbasis_class(int_types:list,boost_types:list) -> str:
    return ''


def emit_ditbasis_class(int_types:list,boost_types:list) -> str:
    return ''


def emit_typedefs(int_types:list,boost_types:list):
    typedefs = [
        cpp.emit_using('bit_perm','quspin::basis::bit_perm<I>','template<typename I>'),
        cpp.emit_using('perm_bit','quspin::basis::perm_bit<I>','template<typename I>'),
        cpp.emit_using('bit_set','quspin::basis::bit_set<I>','template<typename I>'),

        cpp.emit_using('bit_symmetry',
                       'quspin::basis::symmetry<bit_perm<I>,perm_bit<I>,bit_set<I>,npy_cdouble_wrapper>',
                       'template<typename I>'),
        cpp.emit_using('dit_perm','quspin::basis::dit_perm<I>','template<typename I>'),
        cpp.emit_using('perm_dit','quspin::basis::perm_dit<I>','template<typename I>'),
        cpp.emit_using('dit_set','quspin::basis::dit_set<I>','template<typename I>'),
        cpp.emit_using('dit_symmetry',
                       'quspin::basis::symmetry<dit_perm<I>,perm_dit<I>,dit_set<I>,npy_cdouble_wrapper>',
                       'template<typename I>'),
        cpp.emit_using('operator_string','quspin::operator_string<T>','template<typename T>'),
        cpp.emit_using('two_body','quspin::N_body_bit_op<T,2>','template<typename T>'),
        cpp.emit_using('two_body_dit_op','quspin::N_body_dit_op<T,2>','template<typename T>'),
        "// concrete definitions",
    ]
    for ctype,jtype,ktype,bits in (int_types):
       typedefs += [
            cpp.emit_using(f'bit_subspace_{bits}',f'quspin::basis::bit_subspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'bit_fullspace_{bits}',f'quspin::basis::bit_fullspace<{ctype},{jtype}>'),
            cpp.emit_using(f'dit_subspace_{bits}',f'quspin::basis::dit_subspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'dit_fullspace_{bits}',f'quspin::basis::dit_fullspace<{ctype},{jtype}>'),
       ]
    for ctype,jtype,ktype,bits in (boost_types):
       typedefs += [
            cpp.emit_using(f'bit_subspace_{bits}',f'quspin::basis::bit_subspace<{ctype},{jtype},{ktype}>'),
            # cpp.emit_using(f'bit_fullspace_{bits}',f'quspin::basis::bit_fullspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'dit_subspace_{bits}',f'quspin::basis::dit_subspace<{ctype},{jtype},{ktype}>'),
            # cpp.emit_using(f'dit_fullspace_{bits}',f'quspin::basis::dit_fullspace<{ctype},{jtype},{ktype}>'),
       ]
       

    for ctype,jtype,ktype,bits in (int_types+boost_types):
       typedefs += [
            cpp.emit_using(f'symmetric_bitbasis_{bits}',f'quspin::basis::symmetric_basis<bit_subspace_{bits},bit_symmetry<{ctype}>>'),
            cpp.emit_using(f'symmetric_ditbasis_{bits}',f'quspin::basis::symmetric_basis<dit_subspace_{bits},dit_symmetry<{ctype}>>'),
       ]  




    
    return '\n\n'.join(typedefs)

def emit_basis_abi_body(int_types:list,boost_types:list) -> str:
    typedefs=emit_typedefs(int_types,boost_types)
    symmetric_bitbasis_class = emit_symmetric_bitbasis_class(int_types,boost_types)
    symmetric_ditbasis_class = emit_symmetric_ditbasis_class(int_types,boost_types)

    bitbasys_class = emit_bitbasis_class(int_types,boost_types)
    ditbasis_class = emit_ditbasis_class(int_types,boost_types)
    
    return f"""
{typedefs}

// abi class definitions
{symmetric_bitbasis_class}

{symmetric_ditbasis_class}

{bitbasys_class}

{ditbasis_class}

"""

def emit_basis_abi_source(use_boost:bool) -> str:
    if use_boost:
        boost_types = [
            ('quspin::basis::uint1024_t' ,"npy_intp", "int"                   ,1024 ),
            ('quspin::basis::uint4096_t' ,"npy_intp", "int"                   ,4096 ),
            ('quspin::basis::uint16384_t',"npy_intp", "int"                   ,16384),
        ]
    else:
        boost_types = []

    int_types = [
        ('quspin::basis::uint32_t',"npy_intp", "quspin::basis::uint8_t",32),
        ('quspin::basis::uint64_t',"npy_intp", "quspin::basis::uint8_t",64),
    ]    
    basis_abi_body = emit_basis_abi_body(int_types,boost_types)
    if use_boost:
        boost_header = "#define USE_BOOST"
    else:
        boost_header = ""
        
    return f"""#ifndef __QUSPIN_BASIS_ABI__
#define __QUSPIN_BASIS_ABI__
#define __QUSPIN_VERSION__ "{__version__}"
{boost_header}

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <quspin_abi/complex_ops.h>
#include <memory>

namespace quspin {{
    using namespace quspin_abi;
}}

#include <quspin/quspin.h>

namespace quspin_abi {{
{basis_abi_body}
}}
#endif"""

if __name__ == '__main__':
    use_boost=False
    pwd = os.path.split(os.path.abspath(__file__))[0]
    exec(open(os.path.join(pwd,"src","quspin_core","_version.py")).read())
    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_abi','basis_abi.h'),'w') as IO:
        IO.write(emit_basis_abi_source(use_boost))
        