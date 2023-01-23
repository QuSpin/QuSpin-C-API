#ifndef __QUSPIN_CORE_OPERATOR_ABI__
    #define __QUSPIN_CORE_OPERATOR_ABI__

    #include <numpy/ndarrayobject.h>
    #include <numpy/ndarraytypes.h>
    #include <quspin_core_abi/complex_ops.h>

    #include <memory>
    #include <vector>

    #include <quspin/quspin.h>


    namespace quspin_core_abi {
    
    
struct operator_string_args
{
    const int nlocs;
    int * locs;
    int * perms;
    void * datas;
    operator_string_args(
        const int _nlocs,
        int * _locs,
        int * _perms,
        void * _datas) : nlocs(_nlocs)
    {
        locs = _locs;
        perms = _perms;
        datas = _datas;
    }
    ~operator_string_args(){}
};

struct N_body_op_args
{
    int * locs;
    void * data;
    N_body_op_args(
        int * _locs,
        void * _data)
    {
        locs = _locs;
        data = _data;
    }
    ~N_body_op_args(){}
};

enum OPERATOR_TYPES {OP_STRING, OP_TWO_BODY};
    // abi class definitions
    class operator_abi
{
private:
    static size_t generate_type_switch_code(
        OPERATOR_TYPES op_type,
        const int lhss,
        NPY_TYPES T_typenum)
    {
        if(lhss < 2){throw std::out_of_range("value of lhss not found in 2 <= lhss <= 255");}
        switch(op_type)
        {
            case OP_STRING:
                switch(T_typenum)
                {
                    case NPY_INT8:
                        return 0;
                    case NPY_INT16:
                        return 1;
                    case NPY_FLOAT32:
                        return 2;
                    case NPY_FLOAT64:
                        return 3;
                    case NPY_COMPLEX64:
                        return 4;
                    case NPY_COMPLEX128:
                        return 5;
                    default:
                        throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                
                }
            case OP_TWO_BODY:
                if(lhss==2)
                {
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            return 6;
                        case NPY_INT16:
                            return 7;
                        case NPY_FLOAT32:
                            return 8;
                        case NPY_FLOAT64:
                            return 9;
                        case NPY_COMPLEX64:
                            return 10;
                        case NPY_COMPLEX128:
                            return 11;
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
                else{
                    switch(T_typenum)
                    {
                        case NPY_INT8:
                            return 12;
                        case NPY_INT16:
                            return 13;
                        case NPY_FLOAT32:
                            return 14;
                        case NPY_FLOAT64:
                            return 15;
                        case NPY_COMPLEX64:
                            return 16;
                        case NPY_COMPLEX128:
                            return 17;
                        default:
                            throw std::invalid_argument("expecting value of Operator dtype to be in: [int8, int16, float32, float64, complex64, or complex128]");
                    
                    }
                }
            default:
                return -1;
        
        }
    }
    NPY_TYPES T_typenum;
    OPERATOR_TYPES op_type;
    const int nterms;
    const size_t type_switch_code;
    std::vector<quspin::operator_string<npy_int8>> operator_string_npy_int8;
    std::vector<quspin::operator_string<npy_int16>> operator_string_npy_int16;
    std::vector<quspin::operator_string<float>> operator_string_float;
    std::vector<quspin::operator_string<double>> operator_string_double;
    std::vector<quspin::operator_string<npy_cfloat_wrapper>> operator_string_npy_cfloat_wrapper;
    std::vector<quspin::operator_string<npy_cdouble_wrapper>> operator_string_npy_cdouble_wrapper;
    std::vector<quspin::N_body_bit_op<npy_int8,2>> two_body_bit_op_npy_int8;
    std::vector<quspin::N_body_bit_op<npy_int16,2>> two_body_bit_op_npy_int16;
    std::vector<quspin::N_body_bit_op<float,2>> two_body_bit_op_float;
    std::vector<quspin::N_body_bit_op<double,2>> two_body_bit_op_double;
    std::vector<quspin::N_body_bit_op<npy_cfloat_wrapper,2>> two_body_bit_op_npy_cfloat_wrapper;
    std::vector<quspin::N_body_bit_op<npy_cdouble_wrapper,2>> two_body_bit_op_npy_cdouble_wrapper;
    std::vector<quspin::N_body_dit_op<npy_int8,2>> two_body_dit_op_npy_int8;
    std::vector<quspin::N_body_dit_op<npy_int16,2>> two_body_dit_op_npy_int16;
    std::vector<quspin::N_body_dit_op<float,2>> two_body_dit_op_float;
    std::vector<quspin::N_body_dit_op<double,2>> two_body_dit_op_double;
    std::vector<quspin::N_body_dit_op<npy_cfloat_wrapper,2>> two_body_dit_op_npy_cfloat_wrapper;
    std::vector<quspin::N_body_dit_op<npy_cdouble_wrapper,2>> two_body_dit_op_npy_cdouble_wrapper;

public:

    operator_abi(
        NPY_TYPES _T_typenum,
        OPERATOR_TYPES _op_type,
        const int lhss,
        std::vector<std::shared_ptr<void>> _op_args) : 
    T_typenum(_T_typenum),
    op_type(_op_type),
    nterms(_op_args.size()),
    type_switch_code(generate_type_switch_code(_op_type,lhss,_T_typenum))
    {
        switch(type_switch_code)
        {
            case 0:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_npy_int8.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(npy_int8 *)op_arg->datas);
                }
                break;
            case 1:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_npy_int16.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(npy_int16 *)op_arg->datas);
                }
                break;
            case 2:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_float.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(float *)op_arg->datas);
                }
                break;
            case 3:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_double.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(double *)op_arg->datas);
                }
                break;
            case 4:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_npy_cfloat_wrapper.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(npy_cfloat_wrapper *)op_arg->datas);
                }
                break;
            case 5:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<operator_string_args> op_arg = std::reinterpret_pointer_cast<operator_string_args>(_op_arg);
                    operator_string_npy_cdouble_wrapper.emplace_back(lhss,op_arg->nlocs,op_arg->locs,op_arg->perms,(npy_cdouble_wrapper *)op_arg->datas);
                }
                break;
            case 6:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_npy_int8.emplace_back(op_arg->locs,(npy_int8 *)op_arg->data);
                }
                break;
            case 7:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_npy_int16.emplace_back(op_arg->locs,(npy_int16 *)op_arg->data);
                }
                break;
            case 8:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_float.emplace_back(op_arg->locs,(float *)op_arg->data);
                }
                break;
            case 9:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_double.emplace_back(op_arg->locs,(double *)op_arg->data);
                }
                break;
            case 10:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_npy_cfloat_wrapper.emplace_back(op_arg->locs,(npy_cfloat_wrapper *)op_arg->data);
                }
                break;
            case 11:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_bit_op_npy_cdouble_wrapper.emplace_back(op_arg->locs,(npy_cdouble_wrapper *)op_arg->data);
                }
                break;
            case 12:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_npy_int8.emplace_back(lhss,op_arg->locs,(npy_int8 *)op_arg->data);
                }
                break;
            case 13:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_npy_int16.emplace_back(lhss,op_arg->locs,(npy_int16 *)op_arg->data);
                }
                break;
            case 14:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_float.emplace_back(lhss,op_arg->locs,(float *)op_arg->data);
                }
                break;
            case 15:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_double.emplace_back(lhss,op_arg->locs,(double *)op_arg->data);
                }
                break;
            case 16:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_npy_cfloat_wrapper.emplace_back(lhss,op_arg->locs,(npy_cfloat_wrapper *)op_arg->data);
                }
                break;
            case 17:
                for(std::shared_ptr<void>  _op_arg : _op_args){
                    std::shared_ptr<N_body_op_args> op_arg = std::reinterpret_pointer_cast<N_body_op_args>(_op_arg);
                    two_body_dit_op_npy_cdouble_wrapper.emplace_back(lhss,op_arg->locs,(npy_cdouble_wrapper *)op_arg->data);
                }
                break;
            default:
                throw std::runtime_error("this message should not show up.");
        
        }
    }

    ~operator_abi(){}

    inline NPY_TYPES get_T_typenum() const 
    {
        return T_typenum;
    }
    inline OPERATOR_TYPES get_op_type() const 
    {
        return op_type;
    }
    inline int get_nterms() const 
    {
        return nterms;
    }
    void * data()
    {
        switch(type_switch_code)
        {
            case 0:
                return (void *)operator_string_npy_int8.data();
            case 1:
                return (void *)operator_string_npy_int16.data();
            case 2:
                return (void *)operator_string_float.data();
            case 3:
                return (void *)operator_string_double.data();
            case 4:
                return (void *)operator_string_npy_cfloat_wrapper.data();
            case 5:
                return (void *)operator_string_npy_cdouble_wrapper.data();
            case 6:
                return (void *)two_body_bit_op_npy_int8.data();
            case 7:
                return (void *)two_body_bit_op_npy_int16.data();
            case 8:
                return (void *)two_body_bit_op_float.data();
            case 9:
                return (void *)two_body_bit_op_double.data();
            case 10:
                return (void *)two_body_bit_op_npy_cfloat_wrapper.data();
            case 11:
                return (void *)two_body_bit_op_npy_cdouble_wrapper.data();
            case 12:
                return (void *)two_body_dit_op_npy_int8.data();
            case 13:
                return (void *)two_body_dit_op_npy_int16.data();
            case 14:
                return (void *)two_body_dit_op_float.data();
            case 15:
                return (void *)two_body_dit_op_double.data();
            case 16:
                return (void *)two_body_dit_op_npy_cfloat_wrapper.data();
            case 17:
                return (void *)two_body_dit_op_npy_cdouble_wrapper.data();
            default:
                return nullptr;
        
        }
    }

};
    
    }
    #endif