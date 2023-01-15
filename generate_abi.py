import numpy as np
from numpy import int8,int16,int32,int64,float32,float64,complex64,complex128
import os,sys
from typing import NoReturn, Optional
from itertools import product




numpy_ctypes={float32:"float",float64:"double",complex64:"npy_cfloat_wrapper",complex128:"npy_cdouble_wrapper",
                int32:"npy_int32",int64:"npy_int64",int8:"npy_int8",int16:"npy_int16"}


numpy_numtypes={float32:"NPY_FLOAT32",float64:"NPY_FLOAT64",complex64:"NPY_COMPLEX64",complex128:"NPY_COMPLEX128",
                int32:"NPY_INT32",int64:"NPY_INT64",int16:"NPY_INT16",int8:"NPY_INT8"}

real_dtypes={float32:float32, float64:float64, complex64:float32, complex128:float64}

index_types = [int32,int64]
matrix_types = [int8,int16,float32,float64,complex64,complex128]
symmetric_matrix_types = [float32,float64,complex64,complex128]

on_the_fly_types = [float32,float64,complex64,complex128]

class cpp:
    """A collection of methods useful for generating C++ code."""
    
    @staticmethod
    def emit_domain_error(name:str, domain:str)-> str:
        return f'throw std::domain_error("expecting value of \'{name}\' to be in range: {domain}");'

    @staticmethod
    def emit_invalid_argument(name:str, domain:str)-> str:
        return f'throw std::invalid_argument("expecting value of \'{name}\' to be in: [{domain}]");'

    @staticmethod
    def emit_bad_typeid(name:str, excpected_types:str) -> str:
        return f'throw std::bad_typeid("expecting type of \'{name}\' to be one of: {excpected_types}");'
    
    @staticmethod
    def emit_out_of_range(name:str, collection:str) -> str:
        return f'throw std::out_of_range("value of \'{name}\' not found in \'{collection}\'");'

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
        
        body = (f'\n{{\n    {body}\n}}' if len(body)>0 else '{}')
        args = ('\n    '+(',\n    '.join(args)) if args else '')
        if const_method:
            return f'{ret} {name}({args}) const {body}'
        else:
            return f'{ret} {name}({args}){body}'

    @staticmethod
    def emit_destructor(name:str,body:str="") -> str:
        body = body.replace('\n','\n    ')
        body = (f'\n{{\n    {body}\n}}' if len(body)>0 else '{}')
        return f'~{name}(){body}'

    @staticmethod
    def emit_constructor(name:str,args:list[str] = [], body:str = "",preconstruct : Optional[str] = None ) -> str:
        
        body = "    "+body.replace('\n','\n    ')
        args = '\n    '+(',\n    '.join(args))
        if type(preconstruct) == str:
            return f'{name}({args}) : {preconstruct}\n{{\n{body}\n}}'
        else:
           return  f'{name}({args})\n{{\n{body}\n}}'

    
    @staticmethod
    def emit_class(name:str,constructor:str,destructor:str,public_list:list[str]=[],private_list:list[str]=[]) -> str:
        if private_list:
            private_attr = '    ' + ('\n    ').join([b.replace('\n','\n    ') for b in private_list])+ '\n'
        else:
            private_attr = ''
        
        if public_list:
            public_attr = '    ' + ('\n    ').join([b.replace('\n','\n    ') for b in public_list]) + '\n'
        else:
            public_attr = ''
        constructor  = '    ' + constructor.replace('\n','\n    ')
        destructor   = '    ' + destructor.replace('\n','\n    ')

        return f'class {name}\n{{\nprivate:\n{private_attr}\npublic:\n\n{constructor}\n\n{destructor}\n\n{public_attr}\n}};'

class bosonic_basis:
    """Class to generate basis ABI for bosonic basis types
    
    this class is used to generate the quspin-core bosonic_basis ABI (application backend interface).
    The high idea is to put all possible basis+(integer types) + (numpy types) into a single
    C++ class. functions defined as static methods of the cless give a switch_code that is used
    to dispatch the possible template and basis type methods to remove the usage of explicitly 
    basis types as well as the costly-to-compile cython composite types. 
    
    In cython there will only be a thin class wrapper to handle the python/numpy inputs.
    
    bosonic_basis ABI implementation design:
    
    ```C++
    
    ```
    
    """
    def __init__(self,int_types:list,boost_types:list) -> NoReturn:
        self.name = 'bosonic_basis_abi'
        self.int_types = int_types
        self.boost_types = boost_types
        self.bit_term = [
            ('OP_STRING','quspin::operator_string<{}>'),
            ('OP_TWO_BODY','quspin::N_body_bit_op<{},2>')
        ]
        self.dit_term = [
            ('OP_STRING','quspin::operator_string<{}>'),
            ('OP_TWO_BODY','quspin::N_body_dit_op<{},2>')
        ]        
        
    def emit(self)->str:
        return cpp.emit_class(
            self.name,
            self.emit_constructor(),
            self.emit_destructor(),
            **self.emit_attr()
        )  

    def emit_constructor(self)->str:
        args = [
            cpp.emit_declare('_bits','const size_t'),
            cpp.emit_declare('_lhss','const int'),
            cpp.emit_declare('fullspace','const bool'),
            cpp.emit_declare('_symmetry','symmetry_abi&'),
            cpp.emit_declare('Ns = 0','const size_t'),
        ]
        
        generator_body,case_types = self.get_basis_switch_code_body()
        cases = {}
        for  switch_code,(space_type,symmetry_type,basis_type) in case_types.items():
            if len(symmetry_type)>0: 
                cases[switch_code] = (
                '{\n'
                f'    std::shared_ptr<{space_type}> space = std::make_shared<{space_type}>(Ns,_lhss);\n'
                f'    std::shared_ptr<{symmetry_type}> symmetry = std::reinterpret_pointer_cast<{symmetry_type}>(_symmetry.data());\n'
                f'    std::shared_ptr<{basis_type}> _basis_ptr = std::make_shared<{basis_type}>(*symmetry,space);\n'
                f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'
                f'    break;\n'
                '}'
            )
            else:
                cases[switch_code] = (
                '{\n'
                f'    std::shared_ptr<{space_type}> space = std::make_shared<{space_type}>(Ns,_lhss);\n'
                f'    std::shared_ptr<{basis_type}> _basis_ptr = std::make_shared<{basis_type}>(space);\n'
                f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'
                f'    break;\n'
                '}'
            )

        body = cpp.emit_case('basis_switch_code',cases,'throw std::runtime_error("cannot interpret args");')

        return cpp.emit_constructor(self.name,args=args,body=body,
                                    preconstruct='lhss(_lhss),basis_switch_code(generate_basis_switch_code(_bits,_lhss,_symmetry.get(),fullspace))')

    def emit_destructor(self)->str:
        return cpp.emit_destructor(self.name)

    def emit_attr(self)->dict[list,list]:
        private_list = [
            cpp.emit_var('lhss','const int'),
            cpp.emit_var('basis_switch_code','const size_t'),
            cpp.emit_var('basis_ptr','std::shared_ptr<void>'),
            self.emit_generate_basis_switch_code(),
            self.emit_generate_term_switch_code(),
            self.emit_generate_otf_switch_code(),
            self.emit_generate_build_subspace_switch_code(),
        ]
        
        public_list = [
            self.emit_calc_rowptr(),
            self.emit_calc_matrix(),
            self.emit_on_the_fly(),
            self.emit_build_subspace(),
            self.emit_get_state(),
            self.emit_get_index(),
        ]
        return dict(private_list=private_list,public_list=public_list)

    def get_basis_switch_code_body(self) -> tuple[str,dict]:
        
        switch_code = 0
        cases = {}
        case_types = {}
        for ctype,jtype,ktype,bits in self.int_types:
            
            # (space_type,symmetry_type,basis_type) 
            bit_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'bit_subspace_{bits}',f'bit_symmetry<{ctype}>',f'symmetric_bitbasis_{bits}')
            switch_code += 1
            
            bit_no_symmetry_body = f'return (full_space ? {switch_code} : {switch_code+1});'
            case_types[switch_code] = (f'bit_fullspace_{bits}',f'',f'fullspace_bitbasis_{bits}')
            case_types[switch_code+1] = (f'bit_subspace_{bits}',f'',f'subspace_bitbasis_{bits}')
            switch_code += 2

            dit_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'dit_subspace_{bits}',f'dit_symmetry<{ctype}>',f'symmetric_ditbasis_{bits}')
            switch_code += 1
            
            dit_no_symmetry_body = f'return (full_space ? {switch_code} : {switch_code+1});'
            case_types[switch_code] = (f'dit_fullspace_{bits}',f'',f'fullspace_ditbasis_{bits}')
            case_types[switch_code+1] = (f'dit_subspace_{bits}',f'',f'subspace_ditbasis_{bits}')
            switch_code += 2
            
            bit_body = (
                f'if(symmetry)\n'
                f'{{\n'
                f'    {bit_symmetry_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {bit_no_symmetry_body}\n'
                f'}}\n'                    
            )
            
            dit_body = (
                f'if(symmetry)\n'
                f'{{\n'
                f'    {dit_symmetry_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {dit_no_symmetry_body}\n'
                f'}}\n'  
            )
            
            
            
            bit_body = bit_body.replace('\n','\n    ')
            dit_body = dit_body.replace('\n','\n    ')

            cases[bits] = (
                f'if(lhss<2){{return -1;}}\n'
                f'else if(lhss==2)\n'
                f'{{\n'
                f'    {bit_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {dit_body}\n'
                f'}}'
            )
        
        for ctype,jtype,ktype,bits in self.boost_types: 
            # (space_type,symmetry_type,basis_type) 
            bit_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'bit_subspace_{bits}',f'bit_symmetry<{ctype}>',f'symmetric_bitbasis_{bits}')
            switch_code += 1
            
            bit_no_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'bit_subspace_{bits}',f'',f'subspace_bitbasis_{bits}')
            switch_code += 1

            dit_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'dit_subspace_{bits}',f'dit_symmetry<{ctype}>',f'symmetric_ditbasis_{bits}')
            switch_code += 1
            
            dit_no_symmetry_body = f'return (full_space ? -1 : {switch_code});'
            case_types[switch_code] = (f'dit_subspace_{bits}',f'',f'subspace_ditbasis_{bits}')
            switch_code += 1
            
            bit_body = (
                f'if(symmetry)\n'
                f'{{\n'
                f'    {bit_symmetry_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {bit_no_symmetry_body}\n'
                f'}}\n'                    
            )
            
            dit_body = (
                f'if(symmetry)\n'
                f'{{\n'
                f'    {dit_symmetry_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {dit_no_symmetry_body}\n'
                f'}}\n'  
            )
            
            bit_body = bit_body.replace('\n','\n    ')
            dit_body = dit_body.replace('\n','\n    ')

            cases[bits] = (
                f'if(lhss<2){{return -1;}}\n'
                f'else if(lhss==2)\n'
                f'{{\n'
                f'    {bit_body}\n'
                f'}}\n'
                f'else\n'
                f'{{\n'
                f'    {dit_body}\n'
                f'}}'
            )
        
        return cpp.emit_case('bits',cases,'return -1;'),case_types
    
    def get_term_switch_code_body(self) -> tuple[str,dict]:
        cases = {}
        switch_code = 0
        switch_code_types = {}
        
        _,basis_case_types = self.get_basis_switch_code_body()

        for basis_switch_code,(_,_,basis_type) in basis_case_types.items():
            bits_case_body = 'if(0){}\n'

            for index_type in index_types:
                index_type_num = numpy_numtypes[index_type]
                index_ctype = numpy_ctypes[index_type]
                
                matrix_type_cases = {}
                for matrix_type in matrix_types:
                    matrix_ctype = numpy_ctypes[matrix_type]
                    matric_typenum = numpy_numtypes[matrix_type]
                    
                    
                    if 'symmetric' in basis_type and np.issubdtype(matrix_type,np.integer):
                        continue
                        # matrix_type_cases[matric_typenum] = 'return -1;'
                    else:
                        term_cases = {}
                            
                        terms = (self.bit_term if 'bit' in basis_type else self.dit_term)
                        for op_type,term_type in terms:
                            term_cases[op_type] = f'return {switch_code};'
                            switch_code_types[switch_code] = (basis_type,index_ctype,matrix_ctype,term_type.format(matrix_ctype))
                            switch_code += 1
                            
                        matrix_type_cases[matric_typenum] = cpp.emit_case("op_type",term_cases,'return -1;')
                        
                matrix_type_case = cpp.emit_case("T_typenum",matrix_type_cases,
                    default='return -1;'
                )
                matrix_type_case = matrix_type_case.replace('\n','\n    ')
                
                bits_case_body += f'else if(PyArray_EquivTypenums(J_typenum,{index_type_num}))\n{{\n    {matrix_type_case}\n}}\n'

            cases[basis_switch_code] = bits_case_body + "else {return -1;}"
        
        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')

        
        return method_body,switch_code_types

    def get_otf_switch_code_body(self) -> tuple[str,dict]:
        cases = {}
        switch_code = 0
        switch_code_types = {}

        _,basis_case_types = self.get_basis_switch_code_body()
        iterate_typenums = [(T,numpy_numtypes[T]) for T in symmetric_matrix_types]

        for basis_switch_code,(_,_,basis_type) in basis_case_types.items():
            T_cases =  {}
            for (T,T_typenum) in iterate_typenums:
                X_cases = {}
                for (X,X_typenum) in iterate_typenums:
                    R = np.result_type(X,T)
                    for Y,Y_typenum in iterate_typenums:
                        if R == Y:
                            Y_conditional = f'Y_typenum == {Y_typenum}'
                            term_cases = {}
                            
                            terms = (self.bit_term if 'bit' in basis_type else self.dit_term)
                            for op_type,term_type in terms:
                                term_cases[op_type] = f'return {switch_code};'
                                switch_code_types[switch_code] = (basis_type,numpy_ctypes[T], numpy_ctypes[X], numpy_ctypes[Y],term_type.format(numpy_ctypes[T]))
                                switch_code += 1
                            
                            term_body = cpp.emit_case('op_type',term_cases,'return -1;')
                            term_body = term_body.replace('\n','\n    ')
                            X_cases[X_typenum] = (
                                f'if({Y_conditional})\n'
                                f'{{\n'
                                f'    {term_body}\n'
                                f'}}\n'
                                f'else{{return -1;}}'
                            )
                T_cases[T_typenum] = cpp.emit_case("X_typenum",X_cases,"return -1;")
            cases[basis_switch_code] = cpp.emit_case("T_typenum",T_cases,"return -1;")
        
        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')
        return method_body,switch_code_types

    def get_build_subspace_switch_code_body(self) -> tuple[str,dict]:
        cases = {}
        switch_code = 0
        switch_code_types = {}
        
        _,basis_case_types = self.get_basis_switch_code_body()
        for basis_switch_code,(_,_,basis_type) in basis_case_types.items():
            
            if "fullspace" in basis_type: # skip fullspace as they can't be generated. 
                cases[basis_switch_code] = 'return -1;'
                continue
            
            matrix_type_cases = {}
            for matrix_type in matrix_types:
                matrix_ctype = numpy_ctypes[matrix_type]
                matric_typenum = numpy_numtypes[matrix_type]
                
                if 'symmetric' in basis_type and np.issubdtype(matrix_type,np.integer):
                    continue
                    #  matrix_type_cases[matric_typenum] = 'return -1;'
                else:
                    term_cases = {}
                    
                    terms = (self.bit_term if 'bit' in basis_type else self.dit_term)
                    for op_type,term_type in terms:
                        term_cases[op_type] = f'return {switch_code};'
                        switch_code_types[switch_code] = (basis_type,matrix_ctype,term_type.format(matrix_ctype))
                        switch_code += 1
                        
                    matrix_type_cases[matric_typenum] = cpp.emit_case("op_type",term_cases,'return -1;')
                    
            cases[basis_switch_code] = cpp.emit_case("T_typenum",matrix_type_cases,
                default='return -1;'
            )

        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')
        return method_body,switch_code_types

    def emit_generate_basis_switch_code(self)->str:
        body,case_types = self.get_basis_switch_code_body()
        args = [
            cpp.emit_declare('bits','const size_t'),
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('symmetry','void *'),
            cpp.emit_declare('full_space','const bool'),
        ]
        
        return cpp.emit_method('generate_basis_switch_code','static size_t',args,body)

    def emit_generate_term_switch_code(self)->str:
        body,case_types = self.get_term_switch_code_body()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('op_type','OPERATOR_TYPES'),
        ]
        
        return cpp.emit_method('generate_term_switch_code','static size_t',args,body)

    def emit_generate_otf_switch_code(self)->str:
        body,case_types = self.get_otf_switch_code_body()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('X_typenum','NPY_TYPES'),
            cpp.emit_declare('Y_typenum','NPY_TYPES'),
            cpp.emit_declare('op_type','OPERATOR_TYPES'),
        ]
        
        return cpp.emit_method('generate_otf_switch_code','static size_t',args,body)
    
    def emit_generate_build_subspace_switch_code(self)->str:
        body,case_types = self.get_build_subspace_switch_code_body()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('op_type','OPERATOR_TYPES'),
        ]
        
        return cpp.emit_method('generate_build_subspace_switch_code','static size_t',args,body)
    
    def emit_calc_rowptr(self)->str:
        cases = {}
        _,codes = self.get_term_switch_code_body()
        
        for code,(basis_type,J,T,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->calc_rowptr'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'({J}*)rowptr);'
            )
        args = [
            cpp.emit_declare('op','operator_abi'),
            cpp.emit_declare('npy_rowptr','PyArrayObject *'),

        ]
        
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = (
            f'NPY_TYPES T_typenum = op.get_T_typenum();\n'
            f'NPY_TYPES J_typenum = npy_typenum(npy_rowptr);\n'
            f'OPERATOR_TYPES op_type = op.get_op_type();\n'
            f'const int nterms = op.get_nterms();\n'
            f'void * terms = op.data();\n'
            f'void * rowptr = npy_data(npy_rowptr);\n'
            f'const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,op_type);'\
            f'\n{switch}'
        )
        return cpp.emit_method(f'calc_rowptr','size_t',args,body,const_method=True)
    
    def emit_calc_matrix(self)->str:
        cases = {}
        _,codes = self.get_term_switch_code_body()

        
        for code,(basis_type,J,T,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->calc_matrix'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'({T}*)values, '\
                f'({J}*)indices, '\
                f'({J}*)rowptr);'
            )
        args = [
            cpp.emit_declare('op','operator_abi'),
            cpp.emit_declare('npy_values','PyArrayObject *'),
            cpp.emit_declare('npy_indices','PyArrayObject *'),
            cpp.emit_declare('npy_rowptr','PyArrayObject *'),
        ]
        
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = (
            f'NPY_TYPES J_typenum = npy_typenum(npy_indices);\n'
            f'NPY_TYPES T_typenum = npy_typenum(npy_values);\n'
            f'OPERATOR_TYPES op_type = op.get_op_type();\n'
            'const int nterms = op.get_nterms();\n'
            f'void * terms = op.data();\n'
            f'void * values = npy_data(npy_values);\n'
            f'void * indices = npy_data(npy_indices);\n'
            f'void * rowptr = npy_data(npy_rowptr);\n'
            f'const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,op_type);\n{switch}'
        )
        return cpp.emit_method(f'calc_matrix','size_t',args,body,const_method=True)

    def emit_on_the_fly(self)->str:
        cases = {}
        _,codes = self.get_otf_switch_code_body()

        for code,(basis_type,T,X,Y,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->on_the_fly'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'*(const {Y}*)a, '\
                f'(const {X}*)input, '\
                f'*(const {Y}*)b, '\
                f'({Y}*)output);'
            )
            
        args = [
            cpp.emit_declare('op','operator_abi'),
            cpp.emit_declare('npy_a','PyArrayObject *'),
            cpp.emit_declare('npy_input','PyArrayObject *'),
            cpp.emit_declare('npy_b','PyArrayObject *'),
            cpp.emit_declare('npy_output','PyArrayObject *'),
        ]
        
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = (
            f'NPY_TYPES T_typenum = op.get_T_typenum();\n'
            f'NPY_TYPES X_typenum = npy_typenum(npy_input);\n'
            f'NPY_TYPES Y_typenum = npy_typenum(npy_output);\n'
            f'OPERATOR_TYPES op_type = op.get_op_type();\n'
            f'const int nterms = op.get_nterms();\n'
            f'void * terms = op.data();\n'
            f'void * input = npy_data(npy_input);\n'
            f'void * output = npy_data(npy_output);\n'
            f'void * a = npy_data(npy_a);\n'
            f'void * b = npy_data(npy_b);\n'
            f'const size_t switch_code = generate_otf_switch_code(basis_switch_code,T_typenum,X_typenum,Y_typenum,op_type);\n'
            f'{switch}'
        )
        return cpp.emit_method(f'on_the_fly','size_t',args,body,const_method=True)
      
    def emit_build_subspace(self)->str:
        cases = {}
        _,codes = self.get_build_subspace_switch_code_body()

        for code,(basis_type,T,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->build_subspace'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'seed_state, lhss);'
            )
        args = [
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('op_type','OPERATOR_TYPES'),        
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('seed_state','const std::vector<int>&')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = generate_build_subspace_switch_code(basis_switch_code,T_typenum,op_type);\n{switch}'
        return cpp.emit_method(f'build_subspace','size_t',args,body,const_method=False)

    def emit_get_state(self)->str:
        args = [
            cpp.emit_declare('state_index','const npy_intp'),
        ]
        _,case_types = self.get_basis_switch_code_body()
        cases = {}
        for  switch_code,(space_type,symmetry_type,basis_type) in case_types.items():
            cases[switch_code] = (
                f'return std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->space->get_state(state_index).to_vector();'
            )
        body = cpp.emit_case('basis_switch_code',cases,'throw std::runtime_error("invalid basis type");')
        return cpp.emit_method('get_state','std::vector<quspin::basis::dit_integer_t>',args,body,const_method=True)

    def emit_get_index(self)->str:
        args = [
            cpp.emit_declare('state_vector','const std::vector<quspin::basis::dit_integer_t>&'),
        ]

        _,case_types = self.get_basis_switch_code_body()
        cases = {}
        for  switch_code,(space_type,symmetry_type,basis_type) in case_types.items():
            cases[switch_code] = (
                f'return std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->space->get_index(state_vector);'
            )
        body = cpp.emit_case('basis_switch_code',cases,'throw std::runtime_error("invalid basis type");')
        return cpp.emit_method('get_index','npy_intp',args,body,const_method=True)

class anyonic_basis:
    def __init__(self,*args):
        pass

    def emit(self) -> str:
        return ''

def emit_basis_abi_typedefs(int_types:list,boost_types:list):
    typedefs = []
    for ctype,jtype,ktype,bits in (int_types):
       typedefs += [
            cpp.emit_using(f'bit_subspace_{bits}',f'quspin::basis::bit_subspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'dit_subspace_{bits}',f'quspin::basis::dit_subspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'bit_fullspace_{bits}',f'quspin::basis::bit_fullspace<{ctype},{jtype}>'),
            cpp.emit_using(f'dit_fullspace_{bits}',f'quspin::basis::dit_fullspace<{ctype},{jtype}>'),
       ]
    for ctype,jtype,ktype,bits in (boost_types):
       typedefs += [
            cpp.emit_using(f'bit_subspace_{bits}',f'quspin::basis::bit_subspace<{ctype},{jtype},{ktype}>'),
            cpp.emit_using(f'dit_subspace_{bits}',f'quspin::basis::dit_subspace<{ctype},{jtype},{ktype}>'),
       ]
       

    for ctype,jtype,ktype,bits in (int_types+boost_types):
       typedefs += [
            cpp.emit_using(f'symmetric_bitbasis_{bits}',f'quspin::basis::symmetric_basis<bit_subspace_{bits},bit_symmetry<{ctype}>>'),
            cpp.emit_using(f'symmetric_ditbasis_{bits}',f'quspin::basis::symmetric_basis<dit_subspace_{bits},dit_symmetry<{ctype}>>'),
            cpp.emit_using(f'subspace_bitbasis_{bits}',f'quspin::basis::basis<bit_subspace_{bits}>'),
            cpp.emit_using(f'subspace_ditbasis_{bits}',f'quspin::basis::basis<dit_subspace_{bits}>'),
       ]  

    for ctype,jtype,ktype,bits in (int_types):
       typedefs += [
            cpp.emit_using(f'subspace_bitbasis_{bits}',f'quspin::basis::basis<bit_subspace_{bits}>'),
            cpp.emit_using(f'subspace_ditbasis_{bits}',f'quspin::basis::basis<dit_subspace_{bits}>'),
            cpp.emit_using(f'fullspace_bitbasis_{bits}',f'quspin::basis::basis<bit_fullspace_{bits}>'),
            cpp.emit_using(f'fullspace_ditbasis_{bits}',f'quspin::basis::basis<dit_fullspace_{bits}>'),
       ]  

    
    return '\n\n'.join(typedefs)

def emit_basis_abi_body(int_types:list,boost_types:list) -> str:
    typedefs=emit_basis_abi_typedefs(int_types,boost_types)
    bosonic_basis_class = bosonic_basis(int_types,boost_types).emit()
    anyonic_basis_class = anyonic_basis(int_types,boost_types).emit()
    
    return f"""
{typedefs}

// abi class definitions
{bosonic_basis_class}

{anyonic_basis_class}

"""

def emit_basis_abi_source(use_boost:bool) -> str:
    boost_types = get_boost_types(use_boost)

    int_types = [
        ('quspin::basis::uint32_t',"npy_intp", "quspin::basis::uint8_t",32),
        ('quspin::basis::uint64_t',"npy_intp", "quspin::basis::uint8_t",64),
    ]    
    basis_abi_body = emit_basis_abi_body(int_types,boost_types)
   
    return f"""#ifndef __QUSPIN_CORE_BASIS_ABI__
#define __QUSPIN_CORE_BASIS_ABI__

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <quspin_core_abi/complex_ops.h>
#include <quspin_core_abi/symmetry_abi.h>
#include <quspin_core_abi/operator_abi.h>
#include <memory>

#include <quspin/quspin.h>


namespace quspin_core_abi {{
{basis_abi_body}
}}
#endif"""


class operator:
    r"""class to generate QuSpin-Core operator ABI
    
    this class is used to generate the quspin-core operator ABI (application backend interface).
    The high idea is to put all possible (term types) + (numpy types) into a single
    C++ class. functions defined as static methods of the cless give a switch_code that is used
    to dispatch the possible template and term type methods to remove the usage of explicitly 
    term types as well as the costly-to-compile cython composite types. 
    
    In cython there will only be a thin class wrapper to handle the python/numpy inputs. 
    
    operator ABI class design
     
    ```C++
    class operator_abi
    {
    private:
        NPY_TYPES T_typenum;
        OPERATOR_TYPES op_type;
        const int nterms;
        
        const size_t type_switch_code;
        std::vector<quspin::operator_string<npy_int8> operator_string_npy_int8;
        std::vector<quspin::N_body_bit_op<npy_int8,2> two_body_bit_op_npy_int8;
        std::vector<quspin::N_body_dit_op<npy_int8,2> two_body_dit_op_npy_int8;
        std::vector<quspin::operator_string<npy_int16> operator_string_npy_int16;
        std::vector<quspin::N_body_bit_op<npy_int16,2> two_body_bit_op_npy_int16;
        std::vector<quspin::N_body_dit_op<npy_int16,2> two_body_dit_op_npy_int16;
        //... all operator types

        static size_t generate_type_switch_code(
            const int op_type,
            const int lhss, // doesn't matter for operator_string
            NPY_TYPE T_typenum)
        {
            // ... generate switch codes
        }

    public:
        // define two constructors, one for operator_strings and one for two_body_... terms
        operator_abi(
            NPY_TYPES _T_typenum,
            OPERATOR_TYPES _op_type,
            const int lhss,
            std::vector<void*> _op_args //
            ) : T_typenum(_T_typenum), op_type(_op_type), nterms(_op_args.size()), type_switch_code(_op_type,lhss,_T_typenum)
        {
            switch(type_switch_code)
            {
                // initialize for each case
            }
        }
        

        
        ~operator_abi()
        {
            
        }
        
        static NPY_TYPES get_T_typenum() {return T_typenum;}
        static int get_op_type() {return op_type;}
        static int get_nterms(){ return nterms;}         
        void* terms_ptr() {
            switch(type_switch_code)
            {
                // return the appropriate (void*)vector.data();
            }
        }

    };
    ```

    """
    
    def __init__(self) -> NoReturn:
        self.name = 'operator_abi'
        
        self.operator_template = {
            'operator_string':'quspin::operator_string<{}>', 
            'two_body_bit_op':'quspin::N_body_bit_op<{},2>',
            'two_body_dit_op':'quspin::N_body_dit_op<{},2>'
        }
        
        self.operator_typenum = {
            'operator_string':'OP_STRING',
            'two_body_bit_op':'OP_TWO_BODY',
            'two_body_dit_op':'OP_TWO_BODY'
        }
        self.operator_args = {
            'operator_string':'operator_string_args<{}>',
            'two_body_bit_op':'N_body_op_args<{}>',
            'two_body_dit_op':'N_body_op_args<{}>'
        }

    def emit(self)->str:
        return cpp.emit_class(
            self.name,
            self.emit_constructor(),
            self.emit_destructor(),
            **self.emit_attr()
        )  

    def emit_constructor(self):
        args = [
            cpp.emit_declare('_T_typenum','NPY_TYPES'),
            cpp.emit_declare('_op_type','OPERATOR_TYPES'),
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('_op_args','std::vector<void*>')
        ]

        _,case_types = self.get_generate_type_switch_code_body()
        
        cases = {}
        for switch_code,(term_name,arg_type,vec_type,vec_name) in case_types.items():
            
            if 'operator_string' in term_name:
                cases[switch_code] = (
                    f'for(void *  _op_arg : _op_args){{\n'
                    f'    {arg_type} * op_arg = ({arg_type}*)_op_arg;\n'
                    f'    {vec_name}.emplace_back(op_arg->locs,op_arg->perms,op_arg->datas);\n'
                    f'}}'
                )
            elif 'two_body_bit_op' in term_name:
                cases[switch_code] = (
                    f'for(void *  _op_arg : _op_args){{\n'
                    f'    {arg_type} * op_arg = ({arg_type}*)_op_arg;'\
                    f'    {vec_name}.emplace_back(op_arg->locs,op_arg->data);'
                    f'}}'
                )
            elif 'two_body_dit_op' in term_name:
                cases[switch_code] = (
                    f'for(void *  _op_arg : _op_args){{\n'
                    f'    {arg_type} * op_arg = ({arg_type}*)_op_arg;'\
                    f'    {vec_name}.emplace_back(lhss,op_arg->locs,op_arg->data);'
                    f'}}'                )
            else:
                raise NotImplementedError()                

        body = cpp.emit_case('type_switch_code',cases,'throw std::runtime_error("cannot interpret args");')
        
        return cpp.emit_constructor(self.name,args,body,
            preconstruct=(
                '\nT_typenum(_T_typenum),'\
                '\nop_type(_op_type),'\
                '\nnterms(_op_args.size()),'\
                '\ntype_switch_code(generate_type_switch_code(_op_type,lhss,_T_typenum))'
            )
        )
    
    def emit_destructor(self)->str:
        return cpp.emit_destructor(self.name)

    def emit_attr(self)->dict:
        private_list = [
            self.emit_generate_type_switch_code(),       
            cpp.emit_var('T_typenum','NPY_TYPES'),
            cpp.emit_var('op_type','OPERATOR_TYPES'),
            cpp.emit_var('nterms','const int'),
            cpp.emit_var('type_switch_code','const size_t'),  
        ]
        
        _,case_switch = self.get_generate_type_switch_code_body()
        
        for switch_code,(_,_,vec_type,vec_name) in case_switch.items():
            private_list.append(cpp.emit_var(vec_name,vec_type))

        
        public_list = [
            cpp.emit_method('get_T_typenum','inline NPY_TYPES',[],'return T_typenum;',const_method=True),
            cpp.emit_method('get_op_type','inline OPERATOR_TYPES',[],'return op_type;',const_method=True),
            cpp.emit_method('get_nterms','inline int',[],'return nterms;',const_method=True),
            self.emit_data()
        ]
        return dict(private_list=private_list,public_list=public_list)

    def get_generate_type_switch_code_body(self):            
            
        
        op_cases = {}
        case_types = {}
        switch_code = 0
        
        term_name = 'operator_string'
        term_typenum = self.operator_typenum[term_name]
        type_cases = {}
        for matrix_type in matrix_types:
            matrix_typenum = numpy_numtypes[matrix_type]
            ctype = numpy_ctypes[matrix_type]
            
            term_template = self.operator_template[term_name].format(ctype)
            arg_type = self.operator_args[term_name].format(ctype)
            vec_type = f'std::vector<{term_template}>'
            vec_name = f'{term_name}_{ctype}'
            case_types[switch_code] = (term_name,arg_type,vec_type,vec_name)
            type_cases[matrix_typenum] = f'return {switch_code};'
            switch_code += 1
            
        op_cases[term_typenum] = cpp.emit_case('T_typenum',type_cases,'return -1;')
        
        
        
        term_name = 'two_body_bit_op'
        term_typenum = self.operator_typenum[term_name]
        bit_cases = {}
        for matrix_type in matrix_types:
            matrix_typenum = numpy_numtypes[matrix_type]
            ctype = numpy_ctypes[matrix_type]
            
            term_template = self.operator_template[term_name].format(ctype)
            arg_type = self.operator_args[term_name].format(ctype)
            vec_type = f'std::vector<{term_template}>'
            vec_name = f'{term_name}_{ctype}'
            case_types[switch_code] = (term_name,arg_type,vec_type,vec_name)
            bit_cases[matrix_typenum] = f'return {switch_code};'
            switch_code += 1
            
        term_name = 'two_body_dit_op'
        term_typenum = self.operator_typenum[term_name]
        dit_cases = {}
        for matrix_type in matrix_types:
            matrix_typenum = numpy_numtypes[matrix_type]
            ctype = numpy_ctypes[matrix_type]
            
            term_template = self.operator_template[term_name].format(ctype)
            arg_type = self.operator_args[term_name].format(ctype)
            vec_type = f'std::vector<{term_template}>'
            vec_name = f'{term_name}_{ctype}'
            case_types[switch_code] = (term_name,arg_type,vec_type,vec_name)
            dit_cases[matrix_typenum] = f'return {switch_code};'
            switch_code += 1
            
        bit_body = cpp.emit_case('T_typenum',bit_cases,'return -1;')
        dit_body = cpp.emit_case('T_typenum',dit_cases,'return -1;')

        bit_body = bit_body.replace('\n','\n    ')
        dit_body = dit_body.replace('\n','\n    ')
        
        op_cases[term_typenum] = (
            f'if(lhss==2)\n{{\n'
            f'    {bit_body}\n'
            f'}}\n'
            f'else{{\n'
            f'    {dit_body}\n'
            f'}}'
        )
        
        
        
        method_body = f'if(lhss < 2){{return -1;}}\n'+cpp.emit_case('op_type',op_cases,'return -1;')
            
        

        return method_body,case_types
    
    def emit_generate_type_switch_code(self):
        args = [
            cpp.emit_declare('op_type','OPERATOR_TYPES'),
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
        ]
        body,_ = self.get_generate_type_switch_code_body()
        
        return cpp.emit_method('generate_type_switch_code','static size_t',args,body)

    def emit_data(self):        
        _,case_types = self.get_generate_type_switch_code_body()
        cases = {}
        for switch_code,(term_name,arg_type,vec_type,vec_name) in case_types.items():
            cases[switch_code] = (f'return (void *){vec_name}.data();')

        return cpp.emit_method('data','void *',[],cpp.emit_case('type_switch_code',cases,'return nullptr;'))
  
def emit_operator_abi_typedefs():
    return """

template<typename T>
struct operator_string_args { // OP_STRING
    std::vector<int> locs;
    std::vector<std::vector<T>> datas;
    std::vector<std::vector<int>> perms;
    
};

template<typename T>
struct N_body_op_args { // TWO_BODY
    // store each the data as a void*. 
    // data will be recast to the appropriate
    // pointer type and copied later.
    std::vector<int> locs;
    std::vector<T> data;
};

enum OPERATOR_TYPES {OP_STRING, OP_TWO_BODY};

"""

def emit_operator_abi_body() -> str:
    typedefs = emit_operator_abi_typedefs()
    operator_class = operator().emit()
    return f"""
{typedefs}
// abi class definitions
{operator_class}
"""
    
def emit_operator_abi_source(use_boost:bool) -> str:
    return f"""#ifndef __QUSPIN_CORE_OPERATOR_ABI__
#define __QUSPIN_CORE_OPERATOR_ABI__

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <quspin_core_abi/complex_ops.h>

#include <memory>
#include <vector>

#include <quspin/quspin.h>


namespace quspin_core_abi {{
{emit_operator_abi_body()}
}}
#endif"""    


class symmetry:
    def __init__(self,int_types):
        self.name = 'symmetry_abi'
        self.int_types = int_types
        
    def emit(self):
        return cpp.emit_class(
            self.name,
            self.emit_constructor(),
            self.emit_destructor(),
            **self.emit_attr()
        )
        
    def emit_constructor(self):
        _,case_types = self.emit_type_switch_code_body()
        
        cases = {}
        for switch_code,(symmetry_type,lat_symm,loc_symm,lat_args,loc_args) in case_types.items():
            if 'dit' in symmetry_type:
                cases[switch_code] = (
                    f'{{\n'
                    f'    std::vector<{lat_symm}> lat_symm;\n'
                    f'    std::vector<npy_complex_wrapper> lat_char(_lat_char,_lat_char+lat_symm.size());'
                    f'    std::vector<{loc_symm}> loc_symm;\n'
                    f'    std::vector<npy_complex_wrapper> lat_char(_lat_char,_lat_char+lat_symm.size());'
                    f'    for(void * _lat_arg : _lat_args){{\n'
                    f'        {lat_args} * lat_arg =  ({lat_args} *)_lat_arg;\n'
                    f'        lat_symm.emplace_back(lhss,lat_arg->perm);\n'
                    f'    }}\n'
                    f'    for(void * _loc_arg : _loc_args){{\n'
                    f'        {loc_args} *loc_arg =  ({loc_args} *)_loc_arg;\n'
                    f'        loc_symm.emplace_back(loc_arg->perm,loc_arg->locs);\n'
                    f'    }}\n'
                    f'    std::shared_ptr<{symmetry_type}> symmetry = std::make_shared<{symmetry_type}>(lat_symm,lat_char,loc_symm,loc_char);\n'
                    f'    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);\n'
                    f'}}'
                )
            else:
                cases[switch_code] = (
                    f'{{\n'
                    f'    std::vector<{lat_symm}> lat_symm;\n'
                    f'    std::vector<npy_complex_wrapper> lat_char(_lat_char,_lat_char+lat_symm.size());'
                    f'    std::vector<{loc_symm}> loc_symm;\n'
                    f'    std::vector<npy_complex_wrapper> lat_char(_lat_char,_lat_char+lat_symm.size());'
                    f'    for(void * _lat_arg : _lat_args){{\n'
                    f'        {lat_args} * lat_arg =  ({lat_args} *)_lat_arg;\n'
                    f'        lat_symm.emplace_back(lat_arg->perm);\n'
                    f'    }}\n'
                    f'    for(void * _loc_arg : _loc_args){{\n'
                    f'        {loc_args} * loc_arg =  ({loc_args} *)_loc_arg;\n'
                    f'        loc_symm.emplace_back(loc_arg->mask);\n'
                    f'    }}\n'
                    f'    std::shared_ptr<{symmetry_type}> symmetry = std::make_shared<{symmetry_type}>(lat_symm,lat_char,loc_symm,loc_char);\n'
                    f'    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);\n'
                    f'}}'
                )                

        args = [
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('bits','const size_t'),
            cpp.emit_declare('_lat_args','std::vector<void*>'),
            cpp.emit_declare('lat_char','npy_cdouble *'),
            cpp.emit_declare('_loc_args','std::vector<void*>'),
            cpp.emit_declare('loc_char','npy_cdouble *'),

        ]

        switch = cpp.emit_case('type_switch_code',cases,'throw std::runtime_error("cannot parse arguments.");')
        body = (
            f'const size_t type_switch_code = generate_type_switch_code(lhss,bits);\n'
            f'if(_lat_args.size() == 0 || _loc_args.size() == 0){{return;}}\n'
            f'{switch}'
        )
        return cpp.emit_constructor(self.name,args,body=body)
    
    def emit_destructor(self):
        return cpp.emit_destructor(self.name)
    
    def emit_attr(self):
        private_list = [
            cpp.emit_var('symmetry_ptr','std::shared_ptr<void>'),
            self.emit_generate_type_switch_code()
        ]
        public_list = [
            cpp.emit_method('data','std::shared_ptr<void>',[],'return symmetry_ptr;'),
            cpp.emit_method('get','void*',[],'return symmetry_ptr.get();')

        ]
        return dict(
            private_list = private_list,
            public_list = public_list
        )
        
    def emit_type_switch_code_body(self):
        switch_code = 0
        case_types = {}
        bit_cases = {}
        bit_values = ', '.join([str(bits) for _,_,_,bits in self.int_types])

        for ctype,jtype,ktype,bits in self.int_types:
            bit_cases[bits] = f'return {switch_code};'
            case_types[switch_code] = (
                f'bit_symmetry<{ctype}>',
                f'bit_perm<{ctype}>', # lat
                f'perm_bit<{ctype}>',
                f'bit_perm_args',
                f'perm_bit_args',
            )
            switch_code += 1
                        
        dit_cases = {}
        for ctype,jtype,ktype,bits in self.int_types:
            dit_cases[bits] = f'return {switch_code};'
            case_types[switch_code] = (
                f'dit_symmetry<{ctype}>',
                f'dit_perm<{ctype}>',
                f'perm_dit<{ctype}>',
                f'dit_perm_args',
                f'perm_dit_args',
            )          
            switch_code += 1
  
            
        bit_body = cpp.emit_case('bits',bit_cases,cpp.emit_invalid_argument('bits',bit_values))
        dit_body = cpp.emit_case('bits',dit_cases,cpp.emit_invalid_argument('bits',bit_values))

        bit_body = bit_body.replace('\n','\n    ')
        dit_body = dit_body.replace('\n','\n    ')
        
        lhss_exception = cpp.emit_domain_error('lhss','lhss >= 2')
        method_body = (
            f'if(lhss<2){{{lhss_exception}}}\n'
            f'else if(lhss==2)\n'
            f'{{\n'
            f'    {bit_body}\n'
            f'}}\n'
            f'else\n'
            f'{{\n'
            f'    {dit_body}\n'
            f'}}'
        )
        
        return method_body,case_types

    def emit_generate_type_switch_code(self):
        args = [
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('bits','const size_t')
        ]
        
        method_body,_ = self.emit_type_switch_code_body()
        return cpp.emit_method('generate_type_switch_code','static size_t',args,method_body)        

def emit_symmetry_abi_typedefs():
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
        f'struct bit_perm_args{{std::vector<int> perm;}};',
        f'struct perm_bit_args{{std::vector<quspin::basis::dit_integer_t> mask;}};',
        f'struct dit_perm_args{{std::vector<int> perm;}};',
        f'struct perm_dit_args{{std::vector<std::vector<int>> perm;std::vector<int> locs;}};',

    ]
    return '\n\n'.join(typedefs)

def emit_symmetry_abi_body(use_boost:bool):
    boost_types = get_boost_types(use_boost)

    int_types = [
        ('quspin::basis::uint32_t',"npy_intp", "quspin::basis::uint8_t",32),
        ('quspin::basis::uint64_t',"npy_intp", "quspin::basis::uint8_t",64),
    ]    
    typedefs = emit_symmetry_abi_typedefs()
    symmetry_class = symmetry(int_types+boost_types).emit()
    return f"""
{typedefs}

// abi class definitions
{symmetry_class}

"""

def emit_symmetry_abi_source(use_boost):
    symmetry_abi_body = emit_symmetry_abi_body(use_boost)
    return f"""#ifndef __QUSPIN_CORE_SYMMETRY_ABI__
#define __QUSPIN_CORE_SYMMETRY_ABI__

#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <memory>
#include <vector>

#include <quspin/quspin.h>


namespace quspin_core_abi {{
{symmetry_abi_body}
}}
#endif""" 

def emit_abi_source(use_boost):
    boost_flag = ('#define USE_BOOST\n' if use_boost else '')
        
    
    return f"""#ifndef __QUSPIN_CORE_ABI__
#define __QUSPIN_CORE_ABI__
#define __QUSPIN_CORE_VERSION__ "{__version__}"
{boost_flag}
#include <quspin_core_abi/complex_ops.h>
#include <quspin_core_abi/symmetry_abi.h>
#include <quspin_core_abi/operator_abi.h>
#include <quspin_core_abi/basis_abi.h>

#endif
"""

def get_boost_types(use_boost:bool):
    if use_boost:
        boost_types = [
            ('quspin::basis::uint128_t'  ,"npy_intp", "quspin::basis::uint8_t",128 ),
            ('quspin::basis::uint1024_t' ,"npy_intp", "int"                   ,1024 ),
            ('quspin::basis::uint4096_t' ,"npy_intp", "int"                   ,4096 ),
            ('quspin::basis::uint16384_t',"npy_intp", "int"                   ,16384),
        ]
    else:
        boost_types = []
        
    return boost_types



if __name__ == '__main__':
    try:
        use_boost=eval(sys.argv[1])
    except IndexError:
        use_boost = False
        
    pwd = os.path.split(os.path.abspath(__file__))[0]
    exec(open(os.path.join(pwd,"src","quspin_core","_version.py")).read())
    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_core_abi','basis_abi.h'),'w') as IO:
        IO.write(emit_basis_abi_source(use_boost))
        
    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_core_abi','operator_abi.h'),'w') as IO:
        IO.write(emit_operator_abi_source(use_boost))

    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_core_abi','symmetry_abi.h'),'w') as IO:
        IO.write(emit_symmetry_abi_source(use_boost))
    
    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_core_abi','quspin_core_abi.h'),'w') as IO:
        IO.write(emit_abi_source(use_boost))    