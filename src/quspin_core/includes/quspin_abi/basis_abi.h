#ifndef __QUSPIN_BASIS_ABI__
#define __QUSPIN_BASIS_ABI__
#define __QUSPIN_VERSION__ "0.0.1"


#include <numpy/ndarrayobject.h>
#include <numpy/ndarraytypes.h>
#include <quspin_abi/complex_ops.h>
#include <memory>

namespace quspin {
    using namespace quspin_abi;
}

#include <quspin/quspin.h>

namespace quspin_abi {

template<typename I>
using bit_perm = quspin::basis::bit_perm<I>;

template<typename I>
using perm_bit = quspin::basis::perm_bit<I>;

template<typename I>
using bit_set = quspin::basis::bit_set<I>;

template<typename I>
using bit_symmetry = quspin::basis::symmetry<bit_perm<I>,perm_bit<I>,bit_set<I>,npy_cdouble_wrapper>;

template<typename I>
using dit_perm = quspin::basis::dit_perm<I>;

template<typename I>
using perm_dit = quspin::basis::perm_dit<I>;

template<typename I>
using dit_set = quspin::basis::dit_set<I>;

template<typename I>
using dit_symmetry = quspin::basis::symmetry<dit_perm<I>,perm_dit<I>,dit_set<I>,npy_cdouble_wrapper>;

// concrete definitions

using bit_subspace_32 = quspin::basis::bit_subspace<quspin::basis::uint32_t,npy_intp,quspin::basis::uint8_t>;

using symmetric_bitbasis_32 = quspin::basis::symmetric_basis<bit_subspace_32,bit_symmetry<quspin::basis::uint32_t>>;

using bit_subspace_64 = quspin::basis::bit_subspace<quspin::basis::uint64_t,npy_intp,quspin::basis::uint8_t>;

using symmetric_bitbasis_64 = quspin::basis::symmetric_basis<bit_subspace_64,bit_symmetry<quspin::basis::uint64_t>>;

using dit_subspace_32 = quspin::basis::dit_subspace<quspin::basis::uint32_t,npy_intp,quspin::basis::uint8_t>;

using symmetric_ditbasis_32 = quspin::basis::symmetric_basis<dit_subspace_32,dit_symmetry<quspin::basis::uint32_t>>;

using dit_subspace_64 = quspin::basis::dit_subspace<quspin::basis::uint64_t,npy_intp,quspin::basis::uint8_t>;

using symmetric_ditbasis_64 = quspin::basis::symmetric_basis<dit_subspace_64,dit_symmetry<quspin::basis::uint64_t>>;

// abi class definitions
class symmetric_bitbasis_abi
{
private:
    std::shared_ptr<void> basis_ptr;

    const size_t bits;

    static const int lhss = 2;

    static size_t term_switch_code_generator(const size_t bits,NPY_TYPES J_typenum,NPY_TYPES T_typenum)
    {
        switch(bits)
        {
            case 32:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 0;
                        case NPY_FLOAT64:
                            return 1;
                        case NPY_COMPLEX64:
                            return 2;
                        case NPY_COMPLEX128:
                            return 3;
                        default:
                            return -1;
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 4;
                        case NPY_FLOAT64:
                            return 5;
                        case NPY_COMPLEX64:
                            return 6;
                        case NPY_COMPLEX128:
                            return 7;
                        default:
                            return -1;
                    
                    }
                }
                else {return -1;}
            case 64:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 8;
                        case NPY_FLOAT64:
                            return 9;
                        case NPY_COMPLEX64:
                            return 10;
                        case NPY_COMPLEX128:
                            return 11;
                        default:
                            return -1;
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 12;
                        case NPY_FLOAT64:
                            return 13;
                        case NPY_COMPLEX64:
                            return 14;
                        case NPY_COMPLEX128:
                            return 15;
                        default:
                            return -1;
                    
                    }
                }
                else {return -1;}
            default:
                return -1;
        
        }
    }

    static size_t on_the_fly_switch_code_generator(const size_t bits,NPY_TYPES T_typenum,NPY_TYPES X_typenum,NPY_TYPES Y_typenum)
    {
        switch(bits)
        {
            case 32:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT32 ? 0 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 1 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 2 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 3 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT64 ? 4 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 5 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 6 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 7 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX64 ? 8 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 9 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 10 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 11 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX128 ? 12 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 13 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 14 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 15 : -1);
                            default:
                                return -1;
                        
                        }
                    default:
                        return -1;
                
                }
            case 64:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT32 ? 16 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 17 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 18 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 19 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT64 ? 20 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 21 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 22 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 23 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX64 ? 24 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 25 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 26 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 27 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX128 ? 28 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 29 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 30 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 31 : -1);
                            default:
                                return -1;
                        
                        }
                    default:
                        return -1;
                
                }
            default:
                return -1;
        
        }
    }
public:

    symmetric_bitbasis_abi(const size_t _bits,std::shared_ptr<void> symmetry,const size_t Ns_est = 0) : bits(_bits)
    {
        if(0){}
        else if(_bits == 32){
            std::shared_ptr<bit_subspace_32> _space = std::make_shared<bit_subspace_32>(Ns_est);
            std::shared_ptr<bit_symmetry<quspin::basis::uint32_t>> _symmetry = std::reinterpret_pointer_cast<bit_symmetry<quspin::basis::uint32_t>>(symmetry);
            std::shared_ptr<symmetric_bitbasis_32> _basis_ptr = std::make_shared<symmetric_bitbasis_32>(*_symmetry,_space);
            basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
        }
        else if(_bits == 64){
            std::shared_ptr<bit_subspace_64> _space = std::make_shared<bit_subspace_64>(Ns_est);
            std::shared_ptr<bit_symmetry<quspin::basis::uint64_t>> _symmetry = std::reinterpret_pointer_cast<bit_symmetry<quspin::basis::uint64_t>>(symmetry);
            std::shared_ptr<symmetric_bitbasis_64> _basis_ptr = std::make_shared<symmetric_bitbasis_64>(*_symmetry,_space);
            basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
        }
        else{throw std::runtime_error("number of bits not supported");}
        
    }

    ~symmetric_bitbasis_abi()
    {
        
    }

    size_t calc_rowptr(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,void* rowptr) const 
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr);
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr);
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr);
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr);
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr);
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr);
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr);
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr);
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            default:
                return -1;
        
        }
    }

    size_t calc_matrix(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,void* values,void* indices,void* rowptr) const 
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            default:
                return -1;
        
        }
    }

    size_t on_the_fly_operator_string(NPY_TYPES T_typenum,NPY_TYPES X_typenum,NPY_TYPES Y_typenum,void* terms,const int nterms,void* a,void* input,void* b,void* output) const 
    {
        const size_t switch_code = on_the_fly_switch_code_generator(bits,T_typenum,X_typenum,Y_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output);
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output);
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 16:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output);
            case 17:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 18:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 19:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 20:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output);
            case 21:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 22:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 23:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 24:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 25:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 26:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 27:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 28:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 29:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 30:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 31:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            default:
                return -1;
        
        }
    }

    size_t calc_rowptr(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,const std::vector<int>& seed_state)
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 1:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 2:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 3:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 4:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 5:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 6:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 7:
                std::reinterpret_pointer_cast<symmetric_bitbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 8:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 9:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 10:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 11:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 12:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 13:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 14:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 15:
                std::reinterpret_pointer_cast<symmetric_bitbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            default:
                return -1;
        
        }
    }
};

class symmetric_ditbasis_abi
{
private:
    std::shared_ptr<void> basis_ptr;

    const size_t bits;

    const int lhss ;

    static size_t term_switch_code_generator(const size_t bits,NPY_TYPES J_typenum,NPY_TYPES T_typenum)
    {
        switch(bits)
        {
            case 32:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 0;
                        case NPY_FLOAT64:
                            return 1;
                        case NPY_COMPLEX64:
                            return 2;
                        case NPY_COMPLEX128:
                            return 3;
                        default:
                            return -1;
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 4;
                        case NPY_FLOAT64:
                            return 5;
                        case NPY_COMPLEX64:
                            return 6;
                        case NPY_COMPLEX128:
                            return 7;
                        default:
                            return -1;
                    
                    }
                }
                else {return -1;}
            case 64:
                if(0){}
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT32))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 8;
                        case NPY_FLOAT64:
                            return 9;
                        case NPY_COMPLEX64:
                            return 10;
                        case NPY_COMPLEX128:
                            return 11;
                        default:
                            return -1;
                    
                    }
                }
                else if(PyArray_EquivTypenums(J_typenum,NPY_INT64))
                {
                    switch(T_typenum)
                    {
                        case NPY_FLOAT32:
                            return 12;
                        case NPY_FLOAT64:
                            return 13;
                        case NPY_COMPLEX64:
                            return 14;
                        case NPY_COMPLEX128:
                            return 15;
                        default:
                            return -1;
                    
                    }
                }
                else {return -1;}
            default:
                return -1;
        
        }
    }

    static size_t on_the_fly_switch_code_generator(const size_t bits,NPY_TYPES T_typenum,NPY_TYPES X_typenum,NPY_TYPES Y_typenum)
    {
        switch(bits)
        {
            case 32:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT32 ? 0 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 1 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 2 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 3 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT64 ? 4 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 5 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 6 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 7 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX64 ? 8 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 9 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 10 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 11 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX128 ? 12 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 13 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 14 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 15 : -1);
                            default:
                                return -1;
                        
                        }
                    default:
                        return -1;
                
                }
            case 64:
                switch(T_typenum)
                {
                    case NPY_FLOAT32:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT32 ? 16 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 17 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 18 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 19 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_FLOAT64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_FLOAT64 ? 20 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_FLOAT64 ? 21 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 22 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 23 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX64:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX64 ? 24 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 25 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX64 ? 26 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 27 : -1);
                            default:
                                return -1;
                        
                        }
                    case NPY_COMPLEX128:
                        switch(X_typenum)
                        {
                            case NPY_FLOAT32:
                                return ( Y_typenum == NPY_COMPLEX128 ? 28 : -1);
                            case NPY_FLOAT64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 29 : -1);
                            case NPY_COMPLEX64:
                                return ( Y_typenum == NPY_COMPLEX128 ? 30 : -1);
                            case NPY_COMPLEX128:
                                return ( Y_typenum == NPY_COMPLEX128 ? 31 : -1);
                            default:
                                return -1;
                        
                        }
                    default:
                        return -1;
                
                }
            default:
                return -1;
        
        }
    }
public:

    symmetric_ditbasis_abi(const size_t _bits,const int _lhss,std::shared_ptr<void> symmetry,const size_t Ns_est = 0) : bits(_bits), lhss(_lhss)
    {
        if(0){}
        else if(_bits == 32){
            auto _space = std::make_shared<dit_subspace_32>(Ns_est,_lhss);
            auto _symmetry = std::reinterpret_pointer_cast<dit_symmetry<quspin::basis::uint32_t>>(symmetry);
            auto _basis_ptr = std::make_shared<symmetric_ditbasis_32>(*_symmetry,_space);
            basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
        }
        else if(_bits == 64){
            auto _space = std::make_shared<dit_subspace_64>(Ns_est,_lhss);
            auto _symmetry = std::reinterpret_pointer_cast<dit_symmetry<quspin::basis::uint64_t>>(symmetry);
            auto _basis_ptr = std::make_shared<symmetric_ditbasis_64>(*_symmetry,_space);
            basis_ptr = std::reinterpret_pointer_cast<void>(_basis_ptr);
        }
        else{throw std::runtime_error("number of bits not supported");}
        
    }

    ~symmetric_ditbasis_abi()
    {
        
    }

    size_t calc_rowptr(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,void* rowptr) const 
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr);
            case 1:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr);
            case 2:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 3:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 4:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr);
            case 5:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr);
            case 6:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 7:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 8:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int32*)rowptr);
            case 9:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int32*)rowptr);
            case 10:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 11:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int32*)rowptr);
            case 12:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<float>*)terms, nterms, (npy_int64*)rowptr);
            case 13:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<double>*)terms, nterms, (npy_int64*)rowptr);
            case 14:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            case 15:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_rowptr((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_int64*)rowptr);
            default:
                return -1;
        
        }
    }

    size_t calc_matrix(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,void* values,void* indices,void* rowptr) const 
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 1:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 2:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 3:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 4:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 5:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 6:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 7:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 8:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 9:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 10:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 11:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int32*)indices, (npy_int32*)rowptr);
            case 12:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<float>*)terms, nterms, (float*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 13:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<double>*)terms, nterms, (double*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 14:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, (npy_cfloat_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            case 15:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->calc_matrix((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, (npy_cdouble_wrapper*)values, (npy_int64*)indices, (npy_int64*)rowptr);
            default:
                return -1;
        
        }
    }

    size_t on_the_fly_operator_string(NPY_TYPES T_typenum,NPY_TYPES X_typenum,NPY_TYPES Y_typenum,void* terms,const int nterms,void* a,void* input,void* b,void* output) const 
    {
        const size_t switch_code = on_the_fly_switch_code_generator(bits,T_typenum,X_typenum,Y_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output);
            case 1:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 2:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 3:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 4:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output);
            case 5:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 6:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 7:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 8:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 9:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 10:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 11:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 12:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 13:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 14:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 15:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 16:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const float*)a, (const float*)input, *(const float*)b, (float*)output);
            case 17:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 18:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 19:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<float>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 20:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const float*)input, *(const double*)b, (double*)output);
            case 21:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const double*)a, (const double*)input, *(const double*)b, (double*)output);
            case 22:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 23:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<double>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 24:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const float*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 25:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 26:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cfloat_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cfloat_wrapper*)b, (npy_cfloat_wrapper*)output);
            case 27:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 28:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const float*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 29:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const double*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 30:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cfloat_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            case 31:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->on_the_fly((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, *(const npy_cdouble_wrapper*)a, (const npy_cdouble_wrapper*)input, *(const npy_cdouble_wrapper*)b, (npy_cdouble_wrapper*)output);
            default:
                return -1;
        
        }
    }

    size_t calc_rowptr(NPY_TYPES J_typenum,NPY_TYPES T_typenum,void* terms,const int nterms,const std::vector<int>& seed_state)
    {
        const size_t switch_code = term_switch_code_generator(bits,J_typenum,T_typenum);
        switch(switch_code)
        {
            case 0:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 1:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 2:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 3:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 4:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 5:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 6:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 7:
                std::reinterpret_pointer_cast<symmetric_ditbasis_32>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 8:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 9:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 10:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 11:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            case 12:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<float>*)terms, nterms, seed_state, lhss);
            case 13:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<double>*)terms, nterms, seed_state, lhss);
            case 14:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cfloat_wrapper>*)terms, nterms, seed_state, lhss);
            case 15:
                std::reinterpret_pointer_cast<symmetric_ditbasis_64>(basis_ptr)->build_subspace((const quspin::operator_string<npy_cdouble_wrapper>*)terms, nterms, seed_state, lhss);
            default:
                return -1;
        
        }
    }
};






}
#endif