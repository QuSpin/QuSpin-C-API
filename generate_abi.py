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
        args = '\n    '+(',\n    '.join(args))
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
    def __init__(self,int_types:list,boost_types:list) -> NoReturn:
        self.name = 'bosonic_basis_abi'
        self.int_types = int_types
        self.boost_types = boost_types
        self.bit_term = [
            (0,'quspin::operator_string<{T}>'),
            (2,'quspin::N_body_bit_op<{T},2>')
        ]
        self.dit_term = [
            (0,'quspin::operator_string<{T}>'),
            (2,'quspin::N_body_dit_op<{T},2>')
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
            cpp.emit_declare('_symmetry','std::shared_ptr<void>'),
            cpp.emit_declare('Ns = 0','const size_t'),
        ]
        
        generator_body,case_types = self.get_basis_switch_code_switch()
        cases = {}
        for  switch_code,(space_type,symmetry_type,basis_type) in case_types.items():
            if len(symmetry_type)>0: 
                cases[switch_code] = (
                '{\n'
                f'    std::shared_ptr<{space_type}> space = std::make_shared<{space_type}>(Ns,_lhss);\n'\
                f'    std::shared_ptr<{symmetry_type}> symmetry = std::reinterpret_pointer_cast<{symmetry_type}>(_symmetry);\n'\
                f'    std::shared_ptr<{basis_type}> _basis_ptr = std::make_shared<{basis_type}>(*symmetry,space);\n'\
                f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'\
                f'    break;\n'\
                '}'
            )
            else:
                cases[switch_code] = (
                '{\n'
                f'    std::shared_ptr<{space_type}> space = std::make_shared<{space_type}>(Ns,_lhss);\n'\
                f'    std::shared_ptr<{basis_type}> _basis_ptr = std::make_shared<{basis_type}>(space);\n'\
                f'    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);\n'\
                f'    break;\n'\
                '}'
            )

        body = cpp.emit_case('basis_switch_code',cases,'throw std::runtime_error("cannot interpret args");')

        return cpp.emit_constructor(self.name,args=args,body=body,
                                    preconstruct='lhss(_lhss),basis_switch_code(generate_basis_switch_code(_bits,_lhss,_symmetry.get(),fullspace))')

    def emit_destructor(self)->str:
        return cpp.emit_destructor(self.name)

    def emit_attr(self)->tuple[list,list]:
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
        ]
        return dict(private_list=private_list,public_list=public_list)

    def get_basis_switch_code_switch(self):
        
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
                f'if(symmetry)\n'\
                f'{{\n'\
                f'    {bit_symmetry_body}\n'
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {bit_no_symmetry_body}\n'\
                f'}}\n'                    
            )
            
            dit_body = (
                f'if(symmetry)\n'\
                f'{{\n'\
                f'    {dit_symmetry_body}\n'
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {dit_no_symmetry_body}\n'\
                f'}}\n'  
            )
            
            
            
            bit_body = bit_body.replace('\n','\n    ')
            dit_body = dit_body.replace('\n','\n    ')

            cases[bits] = (
                f'if(lhss<2){{return -1;}}\n'\
                f'else if(lhss==2)\n'
                f'{{\n'\
                f'    {bit_body}\n'\
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {dit_body}\n'\
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
                f'if(symmetry)\n'\
                f'{{\n'\
                f'    {bit_symmetry_body}\n'
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {bit_no_symmetry_body}\n'\
                f'}}\n'                    
            )
            
            dit_body = (
                f'if(symmetry)\n'\
                f'{{\n'\
                f'    {dit_symmetry_body}\n'
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {dit_no_symmetry_body}\n'\
                f'}}\n'  
            )
            
            bit_body = bit_body.replace('\n','\n    ')
            dit_body = dit_body.replace('\n','\n    ')

            cases[bits] = (
                f'if(lhss<2){{return -1;}}\n'\
                f'else if(lhss==2)\n'
                f'{{\n'\
                f'    {bit_body}\n'\
                f'}}\n'\
                f'else\n'\
                f'{{\n'\
                f'    {dit_body}\n'\
                f'}}'
            )
        
        return cpp.emit_case('bits',cases,'return -1;'),case_types
    
    def get_term_switch_code_switch(self):
        cases = {}
        switch_code = 0
        switch_code_types = {}
        
        _,basis_case_types = self.get_basis_switch_code_switch()

        for basis_switch_code,(_,_,basis_type) in basis_case_types.items():
            bits_case_body = 'if(0){}\n'

            for index_type in index_types:
                index_type_num = numpy_numtypes[index_type]
                index_ctype = numpy_ctypes[index_type]
                
                matrix_type_cases = {}
                for matrix_type in matrix_types:
                    matrix_ctype = numpy_ctypes[matrix_type]
                    matric_typenum = numpy_numtypes[matrix_type]
                    
                    
                    if 'symmetric' in basis_type:
                        matrix_type_cases[matric_typenum] = 'return -1;'
                    else:
                        term_cases = {}
                            
                        terms = (self.bit_term if 'bit' in basis_type else self.dit_term)
                        for length,term_type in terms:
                            term_cases[length] = f'return {switch_code};'
                            switch_code_types[switch_code] = (basis_type,index_ctype,matrix_ctype,term_type.format(T=matrix_ctype))
                            switch_code += 1
                            
                        matrix_type_cases[matric_typenum] = cpp.emit_case("length",term_cases,'return -1;')
                        
                matrix_type_case = cpp.emit_case("T_typenum",matrix_type_cases,
                    default='return -1;'
                )
                matrix_type_case = matrix_type_case.replace('\n','\n    ')
                
                bits_case_body += f'else if(PyArray_EquivTypenums(J_typenum,{index_type_num}))\n{{\n    {matrix_type_case}\n}}\n'

            cases[basis_switch_code] = bits_case_body + "else {return -1;}"
        
        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')

        
        return method_body,switch_code_types

    def get_otf_switch_code_switch(self):
        cases = {}
        switch_code = 0
        switch_code_types = {}

        _,basis_case_types = self.get_basis_switch_code_switch()
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
                            for length,term_type in terms:
                                term_cases[length] = f'return {switch_code};'
                                switch_code_types[switch_code] = (basis_type,numpy_ctypes[T], numpy_ctypes[X], numpy_ctypes[Y],term_type.format(T=numpy_ctypes[T]))
                                switch_code += 1
                            
                            term_body = cpp.emit_case('length',term_cases,'return -1;')
                            term_body = term_body.replace('\n','\n    ')
                            X_cases[X_typenum] = (
                                f'if({Y_conditional})\n'\
                                f'{{\n'\
                                f'    {term_body}\n'
                                f'}}\n'
                                f'else{{return -1;}}'
                            )
                T_cases[T_typenum] = cpp.emit_case("X_typenum",X_cases,"return -1;")
            cases[basis_switch_code] = cpp.emit_case("T_typenum",T_cases,"return -1;")
        
        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')
        return method_body,switch_code_types

    def get_build_subspace_switch_code_switch(self):
        cases = {}
        switch_code = 0
        switch_code_types = {}
        
        _,basis_case_types = self.get_basis_switch_code_switch()
        for basis_switch_code,(_,_,basis_type) in basis_case_types.items():
            
            if "fullspace" in basis_type: # skip fullspace as they can't be generated. 
                cases[basis_switch_code] = 'return -1;'
                continue
            
            matrix_type_cases = {}
            for matrix_type in matrix_types:
                matrix_ctype = numpy_ctypes[matrix_type]
                matric_typenum = numpy_numtypes[matrix_type]
                
                if 'symmetric' in basis_type:
                     matrix_type_cases[matric_typenum] = 'return -1;'
                else:
                    term_cases = {}
                    
                    terms = (self.bit_term if 'bit' in basis_type else self.dit_term)
                    for length,term_type in terms:
                        term_cases[length] = f'return {switch_code};'
                        switch_code_types[switch_code] = (basis_type,matrix_ctype,term_type.format(T=matrix_ctype))
                        switch_code += 1
                        
                    matrix_type_cases[matric_typenum] = cpp.emit_case("length",term_cases,'return -1;')
                    
            cases[basis_switch_code] = cpp.emit_case("T_typenum",matrix_type_cases,
                default='return -1;'
            )

        method_body = cpp.emit_case("basis_switch_code",cases,'return -1;')
        return method_body,switch_code_types

    def emit_generate_basis_switch_code(self)->str:
        body,case_types = self.get_basis_switch_code_switch()
        args = [
            cpp.emit_declare('bits','const size_t'),
            cpp.emit_declare('lhss','const int'),
            cpp.emit_declare('symmetry','void *'),
            cpp.emit_declare('full_space','const bool'),
        ]
        
        return cpp.emit_method('generate_basis_switch_code','static size_t',args,body)

    def emit_generate_term_switch_code(self)->str:
        body,case_types = self.get_term_switch_code_switch()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int',)
        ]
        
        return cpp.emit_method('generate_term_switch_code','static size_t',args,body)

    def emit_generate_otf_switch_code(self)->str:
        body,case_types = self.get_otf_switch_code_switch()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('X_typenum','NPY_TYPES'),
            cpp.emit_declare('Y_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int',)
        ]
        
        return cpp.emit_method('generate_otf_switch_code','static size_t',args,body)
    
    def emit_generate_build_subspace_switch_code(self)->str:
        body,case_types = self.get_build_subspace_switch_code_switch()
        args = [
            cpp.emit_declare('basis_switch_code','const size_t'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int',)
        ]
        
        return cpp.emit_method('generate_build_subspace_switch_code','static size_t',args,body)
    
    def emit_calc_rowptr(self)->str:
        cases = {}
        _,codes = self.get_term_switch_code_switch()
        
        for code,(basis_type,J,T,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->calc_rowptr'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'({J}*)rowptr);'
            )
        args = [
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('rowptr','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,length);\n{switch}'
        return cpp.emit_method(f'calc_rowptr','size_t',args,body,const_method=True)
    
    def emit_calc_matrix(self):
        cases = {}
        _,codes = self.get_term_switch_code_switch()

        
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
            cpp.emit_declare('J_typenum','NPY_TYPES'),
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('values',"void*"),
            cpp.emit_declare('indices','void*'),
            cpp.emit_declare('rowptr','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,length);\n{switch}'
        return cpp.emit_method(f'calc_matrix','size_t',args,body,const_method=True)

    def emit_on_the_fly(self):
        cases = {}
        _,codes = self.get_otf_switch_code_switch()

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
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('X_typenum','NPY_TYPES'),
            cpp.emit_declare('Y_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int'),
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('a','void*'),
            cpp.emit_declare('input','void*'),
            cpp.emit_declare('b','void*'),
            cpp.emit_declare('output','void*')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = generate_otf_switch_code(basis_switch_code,T_typenum,X_typenum,Y_typenum,length);\n{switch}'
        return cpp.emit_method(f'on_the_fly','size_t',args,body,const_method=True)
      
    def emit_build_subspace(self):
        cases = {}
        _,codes = self.get_build_subspace_switch_code_switch()

        for code,(basis_type,T,term_type) in codes.items():
            cases[code] = (
                f'std::reinterpret_pointer_cast<{basis_type}>(basis_ptr)->build_subspace'\
                f'((const {term_type}*)terms, '\
                f'nterms, '\
                f'seed_state, lhss);'
            )
        args = [
            cpp.emit_declare('T_typenum','NPY_TYPES'),
            cpp.emit_declare('length','const int'),        
            cpp.emit_declare('terms','void*'),
            cpp.emit_declare('nterms','const int'),
            cpp.emit_declare('seed_state','const std::vector<int>&')
        ]
        switch = cpp.emit_case('switch_code',cases,'return -1;')
        body = f'const size_t switch_code = generate_build_subspace_switch_code(basis_switch_code,T_typenum,length);\n{switch}'
        return cpp.emit_method(f'build_subspace','size_t',args,body,const_method=False)

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
        "// concrete definitions",
    ]
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
    typedefs=emit_typedefs(int_types,boost_types)
    bosonic_basis_class = bosonic_basis(int_types,boost_types).emit()
    
    return f"""
{typedefs}

// abi class definitions
{bosonic_basis_class}


"""

def emit_basis_abi_source(use_boost:bool) -> str:
    if use_boost:
        boost_types = [
            ('quspin::basis::uint128_t'  ,"npy_intp", "quspin::basis::uint8_t",128 ),
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
    try:
        use_boost=eval(sys.argv[1])
    except IndexError:
        use_boost = False
        
    pwd = os.path.split(os.path.abspath(__file__))[0]
    exec(open(os.path.join(pwd,"src","quspin_core","_version.py")).read())
    with open(os.path.join(pwd,'src','quspin_core','includes','quspin_abi','basis_abi.h'),'w') as IO:
        IO.write(emit_basis_abi_source(use_boost))
        