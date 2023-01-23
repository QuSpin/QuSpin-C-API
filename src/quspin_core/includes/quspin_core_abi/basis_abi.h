#ifndef __QUSPIN_CORE_BASIS_ABI__
    #define __QUSPIN_CORE_BASIS_ABI__

    #include <numpy/ndarrayobject.h>
    #include <numpy/ndarraytypes.h>
    #include <quspin_core_abi/numpy_interface.h>
    #include <quspin_core_abi/complex_ops.h>
    #include <quspin_core_abi/symmetry_abi.h>
    #include <quspin_core_abi/operator_abi.h>
    #include <memory>

    #include <quspin/quspin.h>


    namespace quspin_core_abi {
    
    using bit_subspace_32 = quspin::basis::bit_subspace<quspin::basis::uint32_t,npy_intp,quspin::basis::uint8_t>;

using dit_subspace_32 = quspin::basis::dit_subspace<quspin::basis::uint32_t,npy_intp,quspin::basis::uint8_t>;

using bit_fullspace_32 = quspin::basis::bit_fullspace<quspin::basis::uint32_t,npy_intp>;

using dit_fullspace_32 = quspin::basis::dit_fullspace<quspin::basis::uint32_t,npy_intp>;

using bit_subspace_64 = quspin::basis::bit_subspace<quspin::basis::uint64_t,npy_intp,quspin::basis::uint8_t>;

using dit_subspace_64 = quspin::basis::dit_subspace<quspin::basis::uint64_t,npy_intp,quspin::basis::uint8_t>;

using bit_fullspace_64 = quspin::basis::bit_fullspace<quspin::basis::uint64_t,npy_intp>;

using dit_fullspace_64 = quspin::basis::dit_fullspace<quspin::basis::uint64_t,npy_intp>;

using symmetric_bitbasis_32 = quspin::basis::symmetric_basis<bit_subspace_32,bit_symmetry<quspin::basis::uint32_t>>;

using symmetric_ditbasis_32 = quspin::basis::symmetric_basis<dit_subspace_32,dit_symmetry<quspin::basis::uint32_t>>;

using subspace_bitbasis_32 = quspin::basis::basis<bit_subspace_32>;

using subspace_ditbasis_32 = quspin::basis::basis<dit_subspace_32>;

using symmetric_bitbasis_64 = quspin::basis::symmetric_basis<bit_subspace_64,bit_symmetry<quspin::basis::uint64_t>>;

using symmetric_ditbasis_64 = quspin::basis::symmetric_basis<dit_subspace_64,dit_symmetry<quspin::basis::uint64_t>>;

using subspace_bitbasis_64 = quspin::basis::basis<bit_subspace_64>;

using subspace_ditbasis_64 = quspin::basis::basis<dit_subspace_64>;

using subspace_bitbasis_32 = quspin::basis::basis<bit_subspace_32>;

using subspace_ditbasis_32 = quspin::basis::basis<dit_subspace_32>;

using fullspace_bitbasis_32 = quspin::basis::basis<bit_fullspace_32>;

using fullspace_ditbasis_32 = quspin::basis::basis<dit_fullspace_32>;

using subspace_bitbasis_64 = quspin::basis::basis<bit_subspace_64>;

using subspace_ditbasis_64 = quspin::basis::basis<dit_subspace_64>;

using fullspace_bitbasis_64 = quspin::basis::basis<bit_fullspace_64>;

using fullspace_ditbasis_64 = quspin::basis::basis<dit_fullspace_64>;

    // abi class definitions
    class bosonic_basis_abi
{
private:
    const int lhss;
    const size_t basis_switch_code;
    std::shared_ptr<void> basis_ptr;
    static size_t generate_basis_switch_code(
        const size_t bits,
        const int lhss,
        void * symmetry,
        const bool full_space)
    {
        switch(bits)
        {
            case 32:
                if(lhss<2){throw std::domain_error("expecting value of lhss to be in range: 1 < lhss < 255");}
                else if(lhss==2)
                {
                    if(symmetry)
                    {
                        return (full_space ? -1 /* error handled in constructor */ : 0);
                    }
                    else
                    {
                        return (full_space ? 1 : 2);
                    }
                    
                }
                else
                {
                    if(symmetry)
                    {
                        return (full_space ? -1 /* error handled in constructor */ : 3);
                    }
                    else
                    {
                        return (full_space ? 4 : 5);
                    }
                    
                }
            case 64:
                if(lhss<2){throw std::domain_error("expecting value of lhss to be in range: 1 < lhss < 255");}
                else if(lhss==2)
                {
                    if(symmetry)
                    {
                        return (full_space ? -1 /* error handled in constructor */ : 6);
                    }
                    else
                    {
                        return (full_space ? 7 : 8);
                    }
                    
                }
                else
                {
                    if(symmetry)
                    {
                        return (full_space ? -1 /* error handled in constructor */ : 9);
                    }
                    else
                    {
                        return (full_space ? 10 : 11);
                    }
                    
                }
            default:
                throw std::invalid_argument("expecting value of bits to be in: [32, 64, 128, 1024, 4096, 16384]");
        
        }
    }
    static size_t generate_term_switch_code(
        const size_t basis_switch_code,
        NPY_TYPES J_typenum,
        NPY_TYPES T_typenum,
        OPERATOR_TYPES op_type)
    {
        switch(basis_switch_code)
        {
            case 0:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 0;
                                case OP_TWO_BODY:
                                    return 1;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 2;
                                case OP_TWO_BODY:
                                    return 3;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 4;
                                case OP_TWO_BODY:
                                    return 5;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 6;
                                case OP_TWO_BODY:
                                    return 7;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 8;
                                case OP_TWO_BODY:
                                    return 9;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 10;
                                case OP_TWO_BODY:
                                    return 11;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 12;
                                case OP_TWO_BODY:
                                    return 13;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 14;
                                case OP_TWO_BODY:
                                    return 15;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 1:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 16;
                                case OP_TWO_BODY:
                                    return 17;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 18;
                                case OP_TWO_BODY:
                                    return 19;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 20;
                                case OP_TWO_BODY:
                                    return 21;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 22;
                                case OP_TWO_BODY:
                                    return 23;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 24;
                                case OP_TWO_BODY:
                                    return 25;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 26;
                                case OP_TWO_BODY:
                                    return 27;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 28;
                                case OP_TWO_BODY:
                                    return 29;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 30;
                                case OP_TWO_BODY:
                                    return 31;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 32;
                                case OP_TWO_BODY:
                                    return 33;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 34;
                                case OP_TWO_BODY:
                                    return 35;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 36;
                                case OP_TWO_BODY:
                                    return 37;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 38;
                                case OP_TWO_BODY:
                                    return 39;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 2:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 40;
                                case OP_TWO_BODY:
                                    return 41;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 42;
                                case OP_TWO_BODY:
                                    return 43;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 44;
                                case OP_TWO_BODY:
                                    return 45;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 46;
                                case OP_TWO_BODY:
                                    return 47;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 48;
                                case OP_TWO_BODY:
                                    return 49;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 50;
                                case OP_TWO_BODY:
                                    return 51;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 52;
                                case OP_TWO_BODY:
                                    return 53;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 54;
                                case OP_TWO_BODY:
                                    return 55;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 56;
                                case OP_TWO_BODY:
                                    return 57;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 58;
                                case OP_TWO_BODY:
                                    return 59;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 60;
                                case OP_TWO_BODY:
                                    return 61;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 62;
                                case OP_TWO_BODY:
                                    return 63;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 3:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 64;
                                case OP_TWO_BODY:
                                    return 65;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 66;
                                case OP_TWO_BODY:
                                    return 67;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 68;
                                case OP_TWO_BODY:
                                    return 69;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 70;
                                case OP_TWO_BODY:
                                    return 71;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 72;
                                case OP_TWO_BODY:
                                    return 73;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 74;
                                case OP_TWO_BODY:
                                    return 75;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 76;
                                case OP_TWO_BODY:
                                    return 77;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 78;
                                case OP_TWO_BODY:
                                    return 79;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 4:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 80;
                                case OP_TWO_BODY:
                                    return 81;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 82;
                                case OP_TWO_BODY:
                                    return 83;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 84;
                                case OP_TWO_BODY:
                                    return 85;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 86;
                                case OP_TWO_BODY:
                                    return 87;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 88;
                                case OP_TWO_BODY:
                                    return 89;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 90;
                                case OP_TWO_BODY:
                                    return 91;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 92;
                                case OP_TWO_BODY:
                                    return 93;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 94;
                                case OP_TWO_BODY:
                                    return 95;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 96;
                                case OP_TWO_BODY:
                                    return 97;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 98;
                                case OP_TWO_BODY:
                                    return 99;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 100;
                                case OP_TWO_BODY:
                                    return 101;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 102;
                                case OP_TWO_BODY:
                                    return 103;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 5:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 104;
                                case OP_TWO_BODY:
                                    return 105;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 106;
                                case OP_TWO_BODY:
                                    return 107;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 108;
                                case OP_TWO_BODY:
                                    return 109;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 110;
                                case OP_TWO_BODY:
                                    return 111;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 112;
                                case OP_TWO_BODY:
                                    return 113;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 114;
                                case OP_TWO_BODY:
                                    return 115;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 116;
                                case OP_TWO_BODY:
                                    return 117;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 118;
                                case OP_TWO_BODY:
                                    return 119;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 120;
                                case OP_TWO_BODY:
                                    return 121;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 122;
                                case OP_TWO_BODY:
                                    return 123;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 124;
                                case OP_TWO_BODY:
                                    return 125;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 126;
                                case OP_TWO_BODY:
                                    return 127;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 6:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 128;
                                case OP_TWO_BODY:
                                    return 129;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 130;
                                case OP_TWO_BODY:
                                    return 131;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 132;
                                case OP_TWO_BODY:
                                    return 133;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 134;
                                case OP_TWO_BODY:
                                    return 135;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 136;
                                case OP_TWO_BODY:
                                    return 137;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 138;
                                case OP_TWO_BODY:
                                    return 139;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 140;
                                case OP_TWO_BODY:
                                    return 141;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 142;
                                case OP_TWO_BODY:
                                    return 143;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 7:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 144;
                                case OP_TWO_BODY:
                                    return 145;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 146;
                                case OP_TWO_BODY:
                                    return 147;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 148;
                                case OP_TWO_BODY:
                                    return 149;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 150;
                                case OP_TWO_BODY:
                                    return 151;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 152;
                                case OP_TWO_BODY:
                                    return 153;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 154;
                                case OP_TWO_BODY:
                                    return 155;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 156;
                                case OP_TWO_BODY:
                                    return 157;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 158;
                                case OP_TWO_BODY:
                                    return 159;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 160;
                                case OP_TWO_BODY:
                                    return 161;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 162;
                                case OP_TWO_BODY:
                                    return 163;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 164;
                                case OP_TWO_BODY:
                                    return 165;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 166;
                                case OP_TWO_BODY:
                                    return 167;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 8:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 168;
                                case OP_TWO_BODY:
                                    return 169;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 170;
                                case OP_TWO_BODY:
                                    return 171;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 172;
                                case OP_TWO_BODY:
                                    return 173;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 174;
                                case OP_TWO_BODY:
                                    return 175;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 176;
                                case OP_TWO_BODY:
                                    return 177;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 178;
                                case OP_TWO_BODY:
                                    return 179;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 180;
                                case OP_TWO_BODY:
                                    return 181;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 182;
                                case OP_TWO_BODY:
                                    return 183;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 184;
                                case OP_TWO_BODY:
                                    return 185;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 186;
                                case OP_TWO_BODY:
                                    return 187;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 188;
                                case OP_TWO_BODY:
                                    return 189;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 190;
                                case OP_TWO_BODY:
                                    return 191;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 9:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 192;
                                case OP_TWO_BODY:
                                    return 193;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 194;
                                case OP_TWO_BODY:
                                    return 195;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 196;
                                case OP_TWO_BODY:
                                    return 197;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 198;
                                case OP_TWO_BODY:
                                    return 199;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 200;
                                case OP_TWO_BODY:
                                    return 201;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 202;
                                case OP_TWO_BODY:
                                    return 203;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 204;
                                case OP_TWO_BODY:
                                    return 205;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 206;
                                case OP_TWO_BODY:
                                    return 207;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 10:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 208;
                                case OP_TWO_BODY:
                                    return 209;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 210;
                                case OP_TWO_BODY:
                                    return 211;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 212;
                                case OP_TWO_BODY:
                                    return 213;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 214;
                                case OP_TWO_BODY:
                                    return 215;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 216;
                                case OP_TWO_BODY:
                                    return 217;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 218;
                                case OP_TWO_BODY:
                                    return 219;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 220;
                                case OP_TWO_BODY:
                                    return 221;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 222;
                                case OP_TWO_BODY:
                                    return 223;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 224;
                                case OP_TWO_BODY:
                                    return 225;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 226;
                                case OP_TWO_BODY:
                                    return 227;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 228;
                                case OP_TWO_BODY:
                                    return 229;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 230;
                                case OP_TWO_BODY:
                                    return 231;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            case 11:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 232;
                                case OP_TWO_BODY:
                                    return 233;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 234;
                                case OP_TWO_BODY:
                                    return 235;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 236;
                                case OP_TWO_BODY:
                                    return 237;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 238;
                                case OP_TWO_BODY:
                                    return 239;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 240;
                                case OP_TWO_BODY:
                                    return 241;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 242;
                                case OP_TWO_BODY:
                                    return 243;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 244;
                                case OP_TWO_BODY:
                                    return 245;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_INT16:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 246;
                                case OP_TWO_BODY:
                                    return 247;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT32:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 248;
                                case OP_TWO_BODY:
                                    return 249;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_FLOAT64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 250;
                                case OP_TWO_BODY:
                                    return 251;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX64:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 252;
                                case OP_TWO_BODY:
                                    return 253;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        case NPY_COMPLEX128:
                            switch(op_type)
                            {
                                case OP_STRING:
                                    return 254;
                                case OP_TWO_BODY:
                                    return 255;
                                default:
                                    throw std::runtime_error("this message should not show up.");
                            
                            }
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else {throw std::invalid_argument("expecting value of Index type to be in: [int32, int64]");}
            default:
                throw std::runtime_error("this message should not show up.");
        
        }
    }
    static size_t generate_otf_switch_code(
        const size_t basis_switch_code,
        NPY_TYPES T_typenum,
        NPY_TYPES X_typenum,
        NPY_TYPES Y_typenum,
        OPERATOR_TYPES op_type)
    {
        switch(basis_switch_code)
        {
            case 0:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 0;
                                        case OP_TWO_BODY:
                                            return 1;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 2;
                                        case OP_TWO_BODY:
                                            return 3;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 4;
                                        case OP_TWO_BODY:
                                            return 5;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 6;
                                        case OP_TWO_BODY:
                                            return 7;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 8;
                                        case OP_TWO_BODY:
                                            return 9;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 10;
                                        case OP_TWO_BODY:
                                            return 11;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 12;
                                        case OP_TWO_BODY:
                                            return 13;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 14;
                                        case OP_TWO_BODY:
                                            return 15;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 16;
                                        case OP_TWO_BODY:
                                            return 17;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 18;
                                        case OP_TWO_BODY:
                                            return 19;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 20;
                                        case OP_TWO_BODY:
                                            return 21;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 22;
                                        case OP_TWO_BODY:
                                            return 23;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 24;
                                        case OP_TWO_BODY:
                                            return 25;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 26;
                                        case OP_TWO_BODY:
                                            return 27;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 28;
                                        case OP_TWO_BODY:
                                            return 29;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 30;
                                        case OP_TWO_BODY:
                                            return 31;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 1:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 32;
                                        case OP_TWO_BODY:
                                            return 33;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 34;
                                        case OP_TWO_BODY:
                                            return 35;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 36;
                                        case OP_TWO_BODY:
                                            return 37;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 38;
                                        case OP_TWO_BODY:
                                            return 39;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 40;
                                        case OP_TWO_BODY:
                                            return 41;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 42;
                                        case OP_TWO_BODY:
                                            return 43;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 44;
                                        case OP_TWO_BODY:
                                            return 45;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 46;
                                        case OP_TWO_BODY:
                                            return 47;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 48;
                                        case OP_TWO_BODY:
                                            return 49;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 50;
                                        case OP_TWO_BODY:
                                            return 51;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 52;
                                        case OP_TWO_BODY:
                                            return 53;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 54;
                                        case OP_TWO_BODY:
                                            return 55;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 56;
                                        case OP_TWO_BODY:
                                            return 57;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 58;
                                        case OP_TWO_BODY:
                                            return 59;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 60;
                                        case OP_TWO_BODY:
                                            return 61;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 62;
                                        case OP_TWO_BODY:
                                            return 63;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 2:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 64;
                                        case OP_TWO_BODY:
                                            return 65;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 66;
                                        case OP_TWO_BODY:
                                            return 67;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 68;
                                        case OP_TWO_BODY:
                                            return 69;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 70;
                                        case OP_TWO_BODY:
                                            return 71;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 72;
                                        case OP_TWO_BODY:
                                            return 73;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 74;
                                        case OP_TWO_BODY:
                                            return 75;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 76;
                                        case OP_TWO_BODY:
                                            return 77;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 78;
                                        case OP_TWO_BODY:
                                            return 79;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 80;
                                        case OP_TWO_BODY:
                                            return 81;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 82;
                                        case OP_TWO_BODY:
                                            return 83;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 84;
                                        case OP_TWO_BODY:
                                            return 85;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 86;
                                        case OP_TWO_BODY:
                                            return 87;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 88;
                                        case OP_TWO_BODY:
                                            return 89;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 90;
                                        case OP_TWO_BODY:
                                            return 91;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 92;
                                        case OP_TWO_BODY:
                                            return 93;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 94;
                                        case OP_TWO_BODY:
                                            return 95;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 3:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 96;
                                        case OP_TWO_BODY:
                                            return 97;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 98;
                                        case OP_TWO_BODY:
                                            return 99;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 100;
                                        case OP_TWO_BODY:
                                            return 101;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 102;
                                        case OP_TWO_BODY:
                                            return 103;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 104;
                                        case OP_TWO_BODY:
                                            return 105;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 106;
                                        case OP_TWO_BODY:
                                            return 107;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 108;
                                        case OP_TWO_BODY:
                                            return 109;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 110;
                                        case OP_TWO_BODY:
                                            return 111;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 112;
                                        case OP_TWO_BODY:
                                            return 113;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 114;
                                        case OP_TWO_BODY:
                                            return 115;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 116;
                                        case OP_TWO_BODY:
                                            return 117;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 118;
                                        case OP_TWO_BODY:
                                            return 119;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 120;
                                        case OP_TWO_BODY:
                                            return 121;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 122;
                                        case OP_TWO_BODY:
                                            return 123;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 124;
                                        case OP_TWO_BODY:
                                            return 125;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 126;
                                        case OP_TWO_BODY:
                                            return 127;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 4:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 128;
                                        case OP_TWO_BODY:
                                            return 129;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 130;
                                        case OP_TWO_BODY:
                                            return 131;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 132;
                                        case OP_TWO_BODY:
                                            return 133;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 134;
                                        case OP_TWO_BODY:
                                            return 135;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 136;
                                        case OP_TWO_BODY:
                                            return 137;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 138;
                                        case OP_TWO_BODY:
                                            return 139;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 140;
                                        case OP_TWO_BODY:
                                            return 141;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 142;
                                        case OP_TWO_BODY:
                                            return 143;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 144;
                                        case OP_TWO_BODY:
                                            return 145;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 146;
                                        case OP_TWO_BODY:
                                            return 147;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 148;
                                        case OP_TWO_BODY:
                                            return 149;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 150;
                                        case OP_TWO_BODY:
                                            return 151;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 152;
                                        case OP_TWO_BODY:
                                            return 153;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 154;
                                        case OP_TWO_BODY:
                                            return 155;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 156;
                                        case OP_TWO_BODY:
                                            return 157;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 158;
                                        case OP_TWO_BODY:
                                            return 159;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 5:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 160;
                                        case OP_TWO_BODY:
                                            return 161;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 162;
                                        case OP_TWO_BODY:
                                            return 163;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 164;
                                        case OP_TWO_BODY:
                                            return 165;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 166;
                                        case OP_TWO_BODY:
                                            return 167;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 168;
                                        case OP_TWO_BODY:
                                            return 169;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 170;
                                        case OP_TWO_BODY:
                                            return 171;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 172;
                                        case OP_TWO_BODY:
                                            return 173;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 174;
                                        case OP_TWO_BODY:
                                            return 175;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 176;
                                        case OP_TWO_BODY:
                                            return 177;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 178;
                                        case OP_TWO_BODY:
                                            return 179;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 180;
                                        case OP_TWO_BODY:
                                            return 181;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 182;
                                        case OP_TWO_BODY:
                                            return 183;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 184;
                                        case OP_TWO_BODY:
                                            return 185;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 186;
                                        case OP_TWO_BODY:
                                            return 187;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 188;
                                        case OP_TWO_BODY:
                                            return 189;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 190;
                                        case OP_TWO_BODY:
                                            return 191;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 6:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 192;
                                        case OP_TWO_BODY:
                                            return 193;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 194;
                                        case OP_TWO_BODY:
                                            return 195;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 196;
                                        case OP_TWO_BODY:
                                            return 197;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 198;
                                        case OP_TWO_BODY:
                                            return 199;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 200;
                                        case OP_TWO_BODY:
                                            return 201;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 202;
                                        case OP_TWO_BODY:
                                            return 203;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 204;
                                        case OP_TWO_BODY:
                                            return 205;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 206;
                                        case OP_TWO_BODY:
                                            return 207;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 208;
                                        case OP_TWO_BODY:
                                            return 209;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 210;
                                        case OP_TWO_BODY:
                                            return 211;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 212;
                                        case OP_TWO_BODY:
                                            return 213;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 214;
                                        case OP_TWO_BODY:
                                            return 215;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 216;
                                        case OP_TWO_BODY:
                                            return 217;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 218;
                                        case OP_TWO_BODY:
                                            return 219;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 220;
                                        case OP_TWO_BODY:
                                            return 221;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 222;
                                        case OP_TWO_BODY:
                                            return 223;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 7:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 224;
                                        case OP_TWO_BODY:
                                            return 225;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 226;
                                        case OP_TWO_BODY:
                                            return 227;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 228;
                                        case OP_TWO_BODY:
                                            return 229;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 230;
                                        case OP_TWO_BODY:
                                            return 231;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 232;
                                        case OP_TWO_BODY:
                                            return 233;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 234;
                                        case OP_TWO_BODY:
                                            return 235;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 236;
                                        case OP_TWO_BODY:
                                            return 237;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 238;
                                        case OP_TWO_BODY:
                                            return 239;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 240;
                                        case OP_TWO_BODY:
                                            return 241;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 242;
                                        case OP_TWO_BODY:
                                            return 243;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 244;
                                        case OP_TWO_BODY:
                                            return 245;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 246;
                                        case OP_TWO_BODY:
                                            return 247;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 248;
                                        case OP_TWO_BODY:
                                            return 249;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 250;
                                        case OP_TWO_BODY:
                                            return 251;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 252;
                                        case OP_TWO_BODY:
                                            return 253;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 254;
                                        case OP_TWO_BODY:
                                            return 255;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 8:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 256;
                                        case OP_TWO_BODY:
                                            return 257;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 258;
                                        case OP_TWO_BODY:
                                            return 259;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 260;
                                        case OP_TWO_BODY:
                                            return 261;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 262;
                                        case OP_TWO_BODY:
                                            return 263;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 264;
                                        case OP_TWO_BODY:
                                            return 265;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 266;
                                        case OP_TWO_BODY:
                                            return 267;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 268;
                                        case OP_TWO_BODY:
                                            return 269;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 270;
                                        case OP_TWO_BODY:
                                            return 271;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 272;
                                        case OP_TWO_BODY:
                                            return 273;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 274;
                                        case OP_TWO_BODY:
                                            return 275;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 276;
                                        case OP_TWO_BODY:
                                            return 277;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 278;
                                        case OP_TWO_BODY:
                                            return 279;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 280;
                                        case OP_TWO_BODY:
                                            return 281;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 282;
                                        case OP_TWO_BODY:
                                            return 283;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 284;
                                        case OP_TWO_BODY:
                                            return 285;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 286;
                                        case OP_TWO_BODY:
                                            return 287;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 9:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 288;
                                        case OP_TWO_BODY:
                                            return 289;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 290;
                                        case OP_TWO_BODY:
                                            return 291;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 292;
                                        case OP_TWO_BODY:
                                            return 293;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 294;
                                        case OP_TWO_BODY:
                                            return 295;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 296;
                                        case OP_TWO_BODY:
                                            return 297;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 298;
                                        case OP_TWO_BODY:
                                            return 299;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 300;
                                        case OP_TWO_BODY:
                                            return 301;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 302;
                                        case OP_TWO_BODY:
                                            return 303;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 304;
                                        case OP_TWO_BODY:
                                            return 305;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 306;
                                        case OP_TWO_BODY:
                                            return 307;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 308;
                                        case OP_TWO_BODY:
                                            return 309;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 310;
                                        case OP_TWO_BODY:
                                            return 311;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 312;
                                        case OP_TWO_BODY:
                                            return 313;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 314;
                                        case OP_TWO_BODY:
                                            return 315;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 316;
                                        case OP_TWO_BODY:
                                            return 317;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 318;
                                        case OP_TWO_BODY:
                                            return 319;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 10:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 320;
                                        case OP_TWO_BODY:
                                            return 321;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 322;
                                        case OP_TWO_BODY:
                                            return 323;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 324;
                                        case OP_TWO_BODY:
                                            return 325;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 326;
                                        case OP_TWO_BODY:
                                            return 327;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 328;
                                        case OP_TWO_BODY:
                                            return 329;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 330;
                                        case OP_TWO_BODY:
                                            return 331;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 332;
                                        case OP_TWO_BODY:
                                            return 333;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 334;
                                        case OP_TWO_BODY:
                                            return 335;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 336;
                                        case OP_TWO_BODY:
                                            return 337;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 338;
                                        case OP_TWO_BODY:
                                            return 339;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 340;
                                        case OP_TWO_BODY:
                                            return 341;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 342;
                                        case OP_TWO_BODY:
                                            return 343;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 344;
                                        case OP_TWO_BODY:
                                            return 345;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 346;
                                        case OP_TWO_BODY:
                                            return 347;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 348;
                                        case OP_TWO_BODY:
                                            return 349;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 350;
                                        case OP_TWO_BODY:
                                            return 351;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 11:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT32)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 352;
                                        case OP_TWO_BODY:
                                            return 353;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float32]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 354;
                                        case OP_TWO_BODY:
                                            return 355;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 356;
                                        case OP_TWO_BODY:
                                            return 357;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 358;
                                        case OP_TWO_BODY:
                                            return 359;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 360;
                                        case OP_TWO_BODY:
                                            return 361;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_FLOAT64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 362;
                                        case OP_TWO_BODY:
                                            return 363;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [float64]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 364;
                                        case OP_TWO_BODY:
                                            return 365;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 366;
                                        case OP_TWO_BODY:
                                            return 367;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 368;
                                        case OP_TWO_BODY:
                                            return 369;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 370;
                                        case OP_TWO_BODY:
                                            return 371;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX64)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 372;
                                        case OP_TWO_BODY:
                                            return 373;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex64]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 374;
                                        case OP_TWO_BODY:
                                            return 375;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 376;
                                        case OP_TWO_BODY:
                                            return 377;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_FLOAT64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 378;
                                        case OP_TWO_BODY:
                                            return 379;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX64:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 380;
                                        case OP_TWO_BODY:
                                            return 381;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            case NPY_COMPLEX128:
                                if(Y_typenum == NPY_COMPLEX128)
                                {
                                    switch(op_type)
                                    {
                                        case OP_STRING:
                                            return 382;
                                        case OP_TWO_BODY:
                                            return 383;
                                        default:
                                            throw std::runtime_error("this message should not show up.");
                                    
                                    }
                                }
                                else{throw std::invalid_argument("expecting value of output dtype to be in: [complex128]");}
                            default:
                                throw std::invalid_argument("expecting value of input dtype to be in: [complex128]");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    static size_t generate_build_subspace_switch_code(
        const size_t basis_switch_code,
        NPY_TYPES T_typenum,
        OPERATOR_TYPES op_type)
    {
        switch(basis_switch_code)
        {
            case 0:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 0;
                            case OP_TWO_BODY:
                                return 1;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 2;
                            case OP_TWO_BODY:
                                return 3;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 4;
                            case OP_TWO_BODY:
                                return 5;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 6;
                            case OP_TWO_BODY:
                                return 7;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 1:
                throw std::runtime_error("there is no subspace to build.");
            case 2:
                switch(T_typenum)
                {
                    case NPY_INT8:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 8;
                            case OP_TWO_BODY:
                                return 9;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_INT16:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 10;
                            case OP_TWO_BODY:
                                return 11;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 12;
                            case OP_TWO_BODY:
                                return 13;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 14;
                            case OP_TWO_BODY:
                                return 15;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 16;
                            case OP_TWO_BODY:
                                return 17;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 18;
                            case OP_TWO_BODY:
                                return 19;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 3:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 20;
                            case OP_TWO_BODY:
                                return 21;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 22;
                            case OP_TWO_BODY:
                                return 23;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 24;
                            case OP_TWO_BODY:
                                return 25;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 26;
                            case OP_TWO_BODY:
                                return 27;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 4:
                throw std::runtime_error("there is no subspace to build.");
            case 5:
                switch(T_typenum)
                {
                    case NPY_INT8:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 28;
                            case OP_TWO_BODY:
                                return 29;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_INT16:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 30;
                            case OP_TWO_BODY:
                                return 31;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 32;
                            case OP_TWO_BODY:
                                return 33;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 34;
                            case OP_TWO_BODY:
                                return 35;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 36;
                            case OP_TWO_BODY:
                                return 37;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 38;
                            case OP_TWO_BODY:
                                return 39;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 6:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 40;
                            case OP_TWO_BODY:
                                return 41;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 42;
                            case OP_TWO_BODY:
                                return 43;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 44;
                            case OP_TWO_BODY:
                                return 45;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 46;
                            case OP_TWO_BODY:
                                return 47;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 7:
                throw std::runtime_error("there is no subspace to build.");
            case 8:
                switch(T_typenum)
                {
                    case NPY_INT8:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 48;
                            case OP_TWO_BODY:
                                return 49;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_INT16:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 50;
                            case OP_TWO_BODY:
                                return 51;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 52;
                            case OP_TWO_BODY:
                                return 53;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 54;
                            case OP_TWO_BODY:
                                return 55;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 56;
                            case OP_TWO_BODY:
                                return 57;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 58;
                            case OP_TWO_BODY:
                                return 59;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case 9:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 60;
                            case OP_TWO_BODY:
                                return 61;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 62;
                            case OP_TWO_BODY:
                                return 63;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 64;
                            case OP_TWO_BODY:
                                return 65;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 66;
                            case OP_TWO_BODY:
                                return 67;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [float32, float64, complex64, or complex128]");
                
                }
            case 10:
                throw std::runtime_error("there is no subspace to build.");
            case 11:
                switch(T_typenum)
                {
                    case NPY_INT8:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 68;
                            case OP_TWO_BODY:
                                return 69;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_INT16:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 70;
                            case OP_TWO_BODY:
                                return 71;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT32:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 72;
                            case OP_TWO_BODY:
                                return 73;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_FLOAT64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 74;
                            case OP_TWO_BODY:
                                return 75;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX64:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 76;
                            case OP_TWO_BODY:
                                return 77;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    case NPY_COMPLEX128:
                        switch(op_type)
                        {
                            case OP_STRING:
                                return 78;
                            case OP_TWO_BODY:
                                return 79;
                            default:
                                throw std::runtime_error("this message should not show up.");
                        
                        }
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            default:
                throw std::runtime_error("this message should not show up.");
        
        }
    }

public:

    bosonic_basis_abi(
        const size_t _bits,
        const int _lhss,
        const bool fullspace,
        symmetry_abi& _symmetry,
        const size_t Ns = 0) : lhss(_lhss),basis_switch_code(generate_basis_switch_code(_bits,_lhss,_symmetry.get(),fullspace))
    {
        switch(basis_switch_code)
        {
            case 0:
                {
                    std::shared_ptr<bit_subspace_32> space = std::make_shared<bit_subspace_32>(Ns,_lhss);
                    std::shared_ptr<bit_symmetry<quspin::basis::uint32_t>> symmetry = std::reinterpret_pointer_cast<bit_symmetry<quspin::basis::uint32_t>>(_symmetry.data());
                    std::shared_ptr<symmetric_bitbasis_32> _basis_ptr = std::make_shared<symmetric_bitbasis_32>(*symmetry,space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 1:
                {
                    std::shared_ptr<bit_fullspace_32> space = std::make_shared<bit_fullspace_32>(Ns,_lhss);
                    std::shared_ptr<fullspace_bitbasis_32> _basis_ptr = std::make_shared<fullspace_bitbasis_32>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 2:
                {
                    std::shared_ptr<bit_subspace_32> space = std::make_shared<bit_subspace_32>(Ns,_lhss);
                    std::shared_ptr<subspace_bitbasis_32> _basis_ptr = std::make_shared<subspace_bitbasis_32>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 3:
                {
                    std::shared_ptr<dit_subspace_32> space = std::make_shared<dit_subspace_32>(Ns,_lhss);
                    std::shared_ptr<dit_symmetry<quspin::basis::uint32_t>> symmetry = std::reinterpret_pointer_cast<dit_symmetry<quspin::basis::uint32_t>>(_symmetry.data());
                    std::shared_ptr<symmetric_ditbasis_32> _basis_ptr = std::make_shared<symmetric_ditbasis_32>(*symmetry,space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 4:
                {
                    std::shared_ptr<dit_fullspace_32> space = std::make_shared<dit_fullspace_32>(Ns,_lhss);
                    std::shared_ptr<fullspace_ditbasis_32> _basis_ptr = std::make_shared<fullspace_ditbasis_32>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 5:
                {
                    std::shared_ptr<dit_subspace_32> space = std::make_shared<dit_subspace_32>(Ns,_lhss);
                    std::shared_ptr<subspace_ditbasis_32> _basis_ptr = std::make_shared<subspace_ditbasis_32>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 6:
                {
                    std::shared_ptr<bit_subspace_64> space = std::make_shared<bit_subspace_64>(Ns,_lhss);
                    std::shared_ptr<bit_symmetry<quspin::basis::uint64_t>> symmetry = std::reinterpret_pointer_cast<bit_symmetry<quspin::basis::uint64_t>>(_symmetry.data());
                    std::shared_ptr<symmetric_bitbasis_64> _basis_ptr = std::make_shared<symmetric_bitbasis_64>(*symmetry,space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 7:
                {
                    std::shared_ptr<bit_fullspace_64> space = std::make_shared<bit_fullspace_64>(Ns,_lhss);
                    std::shared_ptr<fullspace_bitbasis_64> _basis_ptr = std::make_shared<fullspace_bitbasis_64>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 8:
                {
                    std::shared_ptr<bit_subspace_64> space = std::make_shared<bit_subspace_64>(Ns,_lhss);
                    std::shared_ptr<subspace_bitbasis_64> _basis_ptr = std::make_shared<subspace_bitbasis_64>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 9:
                {
                    std::shared_ptr<dit_subspace_64> space = std::make_shared<dit_subspace_64>(Ns,_lhss);
                    std::shared_ptr<dit_symmetry<quspin::basis::uint64_t>> symmetry = std::reinterpret_pointer_cast<dit_symmetry<quspin::basis::uint64_t>>(_symmetry.data());
                    std::shared_ptr<symmetric_ditbasis_64> _basis_ptr = std::make_shared<symmetric_ditbasis_64>(*symmetry,space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 10:
                {
                    std::shared_ptr<dit_fullspace_64> space = std::make_shared<dit_fullspace_64>(Ns,_lhss);
                    std::shared_ptr<fullspace_ditbasis_64> _basis_ptr = std::make_shared<fullspace_ditbasis_64>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            case 11:
                {
                    std::shared_ptr<dit_subspace_64> space = std::make_shared<dit_subspace_64>(Ns,_lhss);
                    std::shared_ptr<subspace_ditbasis_64> _basis_ptr = std::make_shared<subspace_ditbasis_64>(space);
                    basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
                    break;
                }
            default:
                throw std::invalid_argument("symmetry reduced basis require a subspace object, please set fullspace=False.");
        
        }
    }

    ~bosonic_basis_abi(){}

    void calc_rowptr(
        operator_abi op,
        PyArrayObject * npy_rowptr) const 
    {
        NPY_TYPES T_typenum = op.get_T_typenum();
        NPY_TYPES J_typenum = npy_typenum(npy_rowptr);
        OPERATOR_TYPES op_type = op.get_op_type();
        const int nterms = op.get_nterms();
        void * terms = op.data();
        void * rowptr = npy_data(npy_rowptr);
        const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,op_type);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 16:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 17:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 18:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 19:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 20:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 21:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 22:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 23:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 24:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 25:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 26:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 27:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 28:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 29:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 30:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 31:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 32:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 33:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 34:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 35:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 36:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 37:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 38:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 39:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 40:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 41:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 42:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 43:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 44:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 45:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 46:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 47:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 48:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 49:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 50:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 51:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 52:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 53:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 54:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 55:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 56:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 57:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 58:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 59:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 60:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 61:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 62:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 63:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 64:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 65:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 66:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 67:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 68:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 69:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 70:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 71:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 72:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 73:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 74:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 75:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 76:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 77:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 78:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 79:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 80:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 81:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 82:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 83:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 84:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 85:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 86:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 87:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 88:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 89:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 90:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 91:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 92:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 93:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 94:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 95:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 96:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 97:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 98:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 99:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 100:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 101:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 102:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 103:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 104:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 105:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 106:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 107:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 108:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 109:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 110:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 111:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 112:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 113:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 114:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 115:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 116:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 117:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 118:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 119:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 120:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 121:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 122:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 123:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 124:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 125:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 126:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 127:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 128:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 129:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 130:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 131:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 132:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 133:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 134:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 135:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 136:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 137:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 138:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 139:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 140:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 141:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 142:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 143:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 144:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 145:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 146:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 147:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 148:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 149:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 150:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 151:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 152:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 153:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 154:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 155:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 156:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 157:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 158:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 159:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 160:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 161:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 162:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 163:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 164:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 165:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 166:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 167:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 168:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 169:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 170:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 171:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 172:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 173:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 174:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 175:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 176:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 177:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 178:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 179:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 180:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 181:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 182:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 183:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 184:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 185:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 186:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 187:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 188:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 189:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 190:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 191:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 192:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 193:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 194:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 195:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 196:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 197:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 198:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 199:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 200:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 201:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 202:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 203:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 204:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 205:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 206:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 207:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 208:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 209:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 210:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 211:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 212:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 213:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 214:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 215:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 216:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 217:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 218:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 219:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 220:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 221:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 222:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 223:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 224:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 225:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 226:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 227:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 228:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 229:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 230:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 231:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 232:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int32*)rowptr); break;
            case 233:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 234:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int32*)rowptr); break;
            case 235:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 236:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr); break;
            case 237:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 238:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr); break;
            case 239:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 240:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 241:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 242:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr); break;
            case 243:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int32*)rowptr); break;
            case 244:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int64*)rowptr); break;
            case 245:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 246:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int64*)rowptr); break;
            case 247:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 248:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr); break;
            case 249:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<float,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 250:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr); break;
            case 251:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<double,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 252:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 253:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            case 254:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr); break;
            case 255:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_int64*)rowptr); break;
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    void calc_matrix(
        operator_abi op,
        PyArrayObject * npy_values,
        PyArrayObject * npy_indices,
        PyArrayObject * npy_rowptr) const 
    {
        NPY_TYPES J_typenum = npy_typenum(npy_indices);
        NPY_TYPES T_typenum = npy_typenum(npy_values);
        OPERATOR_TYPES op_type = op.get_op_type();
        const int nterms = op.get_nterms();
        void * terms = op.data();
        void * values = npy_data(npy_values);
        void * indices = npy_data(npy_indices);
        void * rowptr = npy_data(npy_rowptr);
        const size_t switch_code = generate_term_switch_code(basis_switch_code,J_typenum,T_typenum,op_type);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 16:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 17:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 18:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 19:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 20:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 21:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 22:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 23:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 24:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 25:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 26:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 27:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 28:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 29:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 30:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 31:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 32:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 33:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 34:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 35:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 36:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 37:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 38:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 39:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 40:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 41:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 42:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 43:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 44:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 45:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 46:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 47:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 48:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 49:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 50:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 51:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 52:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 53:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 54:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 55:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 56:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 57:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 58:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 59:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 60:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 61:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 62:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 63:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 64:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 65:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 66:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 67:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 68:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 69:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 70:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 71:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 72:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 73:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 74:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 75:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 76:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 77:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 78:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 79:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 80:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 81:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 82:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 83:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 84:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 85:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 86:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 87:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 88:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 89:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 90:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 91:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 92:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 93:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 94:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 95:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 96:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 97:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 98:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 99:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 100:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 101:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 102:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 103:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 104:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 105:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 106:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 107:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 108:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 109:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 110:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 111:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 112:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 113:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 114:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 115:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 116:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 117:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 118:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 119:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 120:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 121:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 122:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 123:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 124:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 125:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 126:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 127:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 128:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 129:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 130:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 131:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 132:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 133:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 134:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 135:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 136:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 137:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 138:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 139:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 140:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 141:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 142:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 143:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 144:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 145:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 146:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 147:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 148:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 149:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 150:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 151:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 152:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 153:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 154:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 155:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 156:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 157:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 158:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 159:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 160:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 161:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 162:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 163:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 164:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 165:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 166:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 167:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 168:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 169:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 170:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 171:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 172:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 173:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 174:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 175:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 176:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 177:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 178:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 179:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 180:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 181:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 182:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 183:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 184:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 185:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 186:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 187:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 188:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 189:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 190:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 191:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 192:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 193:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 194:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 195:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 196:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 197:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 198:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 199:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 200:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 201:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 202:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 203:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 204:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 205:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 206:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 207:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 208:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 209:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 210:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 211:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 212:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 213:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 214:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 215:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 216:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 217:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 218:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 219:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 220:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 221:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 222:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 223:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 224:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 225:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 226:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 227:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 228:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 229:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 230:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 231:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 232:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 233:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 234:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 235:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 236:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 237:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 238:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 239:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 240:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 241:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 242:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 243:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr); break;
            case 244:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int8>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 245:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, (npy_int8*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 246:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_int16>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 247:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, (npy_int16*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 248:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 249:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<float,2>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 250:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 251:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<double,2>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 252:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 253:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 254:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            case 255:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr); break;
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    void on_the_fly(
        operator_abi op,
        PyArrayObject * npy_a,
        PyArrayObject * npy_input,
        PyArrayObject * npy_b,
        PyArrayObject * npy_output) const 
    {
        NPY_TYPES T_typenum = op.get_T_typenum();
        NPY_TYPES X_typenum = npy_typenum(npy_input);
        NPY_TYPES Y_typenum = npy_typenum(npy_output);
        OPERATOR_TYPES op_type = op.get_op_type();
        const int nterms = op.get_nterms();
        void * terms = op.data();
        void * input = npy_data(npy_input);
        void * output = npy_data(npy_output);
        void * a = npy_data(npy_a);
        void * b = npy_data(npy_b);
        const size_t switch_code = generate_otf_switch_code(basis_switch_code,T_typenum,X_typenum,Y_typenum,op_type);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 16:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 17:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 18:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 19:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 20:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 21:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 22:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 23:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 24:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 25:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 26:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 27:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 28:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 29:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 30:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 31:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 32:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 33:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 34:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 35:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 36:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 37:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 38:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 39:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 40:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 41:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 42:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 43:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 44:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 45:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 46:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 47:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 48:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 49:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 50:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 51:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 52:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 53:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 54:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 55:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 56:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 57:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 58:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 59:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 60:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 61:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 62:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 63:
                std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 64:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 65:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 66:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 67:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 68:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 69:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 70:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 71:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 72:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 73:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 74:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 75:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 76:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 77:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 78:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 79:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 80:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 81:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 82:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 83:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 84:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 85:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 86:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 87:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 88:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 89:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 90:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 91:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 92:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 93:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 94:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 95:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 96:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 97:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 98:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 99:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 100:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 101:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 102:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 103:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 104:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 105:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 106:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 107:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 108:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 109:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 110:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 111:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 112:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 113:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 114:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 115:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 116:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 117:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 118:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 119:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 120:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 121:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 122:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 123:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 124:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 125:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 126:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 127:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 128:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 129:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 130:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 131:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 132:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 133:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 134:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 135:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 136:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 137:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 138:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 139:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 140:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 141:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 142:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 143:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 144:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 145:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 146:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 147:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 148:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 149:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 150:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 151:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 152:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 153:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 154:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 155:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 156:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 157:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 158:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 159:
                std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 160:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 161:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 162:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 163:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 164:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 165:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 166:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 167:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 168:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 169:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 170:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 171:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 172:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 173:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 174:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 175:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 176:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 177:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 178:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 179:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 180:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 181:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 182:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 183:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 184:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 185:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 186:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 187:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 188:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 189:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 190:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 191:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 192:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 193:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 194:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 195:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 196:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 197:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 198:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 199:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 200:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 201:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 202:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 203:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 204:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 205:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 206:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 207:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 208:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 209:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 210:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 211:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 212:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 213:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 214:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 215:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 216:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 217:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 218:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 219:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 220:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 221:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 222:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 223:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 224:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 225:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 226:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 227:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 228:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 229:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 230:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 231:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 232:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 233:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 234:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 235:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 236:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 237:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 238:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 239:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 240:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 241:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 242:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 243:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 244:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 245:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 246:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 247:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 248:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 249:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 250:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 251:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 252:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 253:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 254:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 255:
                std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 256:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 257:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 258:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 259:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 260:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 261:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 262:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 263:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 264:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 265:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 266:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 267:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 268:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 269:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 270:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 271:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 272:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 273:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 274:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 275:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 276:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 277:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 278:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 279:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 280:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 281:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 282:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 283:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 284:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 285:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 286:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 287:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 288:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 289:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 290:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 291:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 292:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 293:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 294:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 295:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 296:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 297:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 298:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 299:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 300:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 301:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 302:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 303:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 304:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 305:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 306:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 307:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 308:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 309:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 310:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 311:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 312:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 313:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 314:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 315:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 316:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 317:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 318:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 319:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 320:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 321:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 322:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 323:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 324:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 325:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 326:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 327:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 328:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 329:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 330:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 331:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 332:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 333:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 334:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 335:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 336:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 337:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 338:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 339:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 340:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 341:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 342:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 343:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 344:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 345:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 346:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 347:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 348:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 349:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 350:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 351:
                std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 352:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 353:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output); break;
            case 354:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 355:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 356:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 357:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 358:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 359:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<float,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 360:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 361:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output); break;
            case 362:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 363:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output); break;
            case 364:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 365:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 366:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 367:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<double,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 368:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 369:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 370:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 371:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 372:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 373:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output); break;
            case 374:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 375:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 376:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 377:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 378:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 379:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 380:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 381:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 382:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            case 383:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output); break;
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    void build_subspace(
        NPY_TYPES T_typenum,
        OPERATOR_TYPES op_type,
        void* terms,
        const int nterms,
        const std::vector<int>& seed_state)
    {
        const size_t switch_code = generate_build_subspace_switch_code(basis_switch_code,T_typenum,op_type);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 8:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int8>*)terms, nterms, seed_state, lhss); break;
            case 9:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, seed_state, lhss); break;
            case 10:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int16>*)terms, nterms, seed_state, lhss); break;
            case 11:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, seed_state, lhss); break;
            case 12:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 13:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 14:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 15:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 16:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 17:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 18:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 19:
                std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 20:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 21:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 22:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 23:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 24:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 25:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 26:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 27:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 28:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int8>*)terms, nterms, seed_state, lhss); break;
            case 29:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, seed_state, lhss); break;
            case 30:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int16>*)terms, nterms, seed_state, lhss); break;
            case 31:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, seed_state, lhss); break;
            case 32:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 33:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 34:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 35:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 36:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 37:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 38:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 39:
                std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 40:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 41:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 42:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 43:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 44:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 45:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 46:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 47:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 48:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int8>*)terms, nterms, seed_state, lhss); break;
            case 49:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_int8,2>*)terms, nterms, seed_state, lhss); break;
            case 50:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int16>*)terms, nterms, seed_state, lhss); break;
            case 51:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_int16,2>*)terms, nterms, seed_state, lhss); break;
            case 52:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 53:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 54:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 55:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 56:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 57:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 58:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 59:
                std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_bit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 60:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 61:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 62:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 63:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 64:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 65:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 66:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 67:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 68:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int8>*)terms, nterms, seed_state, lhss); break;
            case 69:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_int8,2>*)terms, nterms, seed_state, lhss); break;
            case 70:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_int16>*)terms, nterms, seed_state, lhss); break;
            case 71:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_int16,2>*)terms, nterms, seed_state, lhss); break;
            case 72:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss); break;
            case 73:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<float,2>*)terms, nterms, seed_state, lhss); break;
            case 74:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss); break;
            case 75:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<double,2>*)terms, nterms, seed_state, lhss); break;
            case 76:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 77:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cfloat_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            case 78:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss); break;
            case 79:
                std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->build_subspace((const quspin::N_body_dit_op<npy_cdouble_wrapper,2>*)terms, nterms, seed_state, lhss); break;
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    std::vector<quspin::basis::dit_integer_t> get_state(
        const npy_intp state_index) const 
    {
        switch(basis_switch_code)
        {
            case 0:
                return std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 1:
                return std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 2:
                return std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 3:
                return std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 4:
                return std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 5:
                return std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->space->get_state(state_index).to_vector();
            case 6:
                return std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            case 7:
                return std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            case 8:
                return std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            case 9:
                return std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            case 10:
                return std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            case 11:
                return std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->space->get_state(state_index).to_vector();
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }
    npy_intp get_index(
        const std::vector<quspin::basis::dit_integer_t>& state_vector) const 
    {
        switch(basis_switch_code)
        {
            case 0:
                return std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 1:
                return std::reinterpret_pointer_cast<fullspace_bitbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 2:
                return std::reinterpret_pointer_cast<subspace_bitbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 3:
                return std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 4:
                return std::reinterpret_pointer_cast<fullspace_ditbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 5:
                return std::reinterpret_pointer_cast<subspace_ditbasis_32>(basis_ptr)->space->get_index(state_vector);
            case 6:
                return std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->space->get_index(state_vector);
            case 7:
                return std::reinterpret_pointer_cast<fullspace_bitbasis_64>(basis_ptr)->space->get_index(state_vector);
            case 8:
                return std::reinterpret_pointer_cast<subspace_bitbasis_64>(basis_ptr)->space->get_index(state_vector);
            case 9:
                return std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->space->get_index(state_vector);
            case 10:
                return std::reinterpret_pointer_cast<fullspace_ditbasis_64>(basis_ptr)->space->get_index(state_vector);
            case 11:
                return std::reinterpret_pointer_cast<subspace_ditbasis_64>(basis_ptr)->space->get_index(state_vector);
            default:
                throw std::runtime_error("this message should not appear.");
        
        }
    }

};

    

    
    }
    #endif