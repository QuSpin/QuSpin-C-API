#ifndef __QUSPIN_CORE_SYMMETRY_ABI__
    #define __QUSPIN_CORE_SYMMETRY_ABI__

    #include <numpy/ndarrayobject.h>
    #include <numpy/ndarraytypes.h>
    #include <complex>
    #include <memory>
    #include <vector>
    
    #include <quspin/quspin.h>


    namespace quspin_core_abi {
    
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

struct bit_perm_args
{
    std::vector<int> perm;
    bit_perm_args(
        std::vector<int> _perm)
    {
        perm = _perm;
    }
    ~bit_perm_args(){}
};

struct dit_perm_args
{
    std::vector<int> perm;
    dit_perm_args(
        std::vector<int> _perm)
    {
        perm = _perm;
    }
    ~dit_perm_args(){}
};

struct perm_bit_args
{
    std::vector<int> mask;
    perm_bit_args(
        std::vector<int> _mask)
    {
        mask = _mask;
    }
    ~perm_bit_args(){}
};

struct perm_dit_args
{
    std::vector<std::vector<int>> perms;
    std::vector<int> locs;
    perm_dit_args(
        std::vector<std::vector<int>> _perms,
        std::vector<int> _locs)
    {
        perms = _perms;
        locs = _locs;
    }
    ~perm_dit_args(){}
};

    // abi class definitions
    class symmetry_abi
{
private:
    std::shared_ptr<void> symmetry_ptr;
    static size_t generate_type_switch_code(
        const int lhss,
        const size_t bits)
    {
        if(lhss<2){throw std::domain_error("expecting value of lhss to be in range: lhss >= 2");}
        else if(lhss==2)
        {
            switch(bits)
            {
                case 32:
                    return 0;
                case 64:
                    return 1;
                default:
                    throw std::invalid_argument("expecting value of bits to be in: [32, 64]");
            
            }
        }
        else
        {
            switch(bits)
            {
                case 32:
                    return 2;
                case 64:
                    return 3;
                default:
                    throw std::invalid_argument("expecting value of bits to be in: [32, 64]");
            
            }
        }
    }

public:

    symmetry_abi(
        const int lhss,
        const size_t bits,
        std::vector<std::shared_ptr<void>> _lat_args,
        std::complex<double> * _lat_char,
        std::vector<std::shared_ptr<void>> _loc_args,
        std::complex<double> * _loc_char)
    {
        const size_t type_switch_code = generate_type_switch_code(lhss,bits);
        if(_lat_args.size() == 0 || _loc_args.size() == 0){return;}
        switch(type_switch_code)
        {
            case 0:
                {
                    std::vector<bit_perm<quspin::basis::uint32_t>> lat_symm;
                    std::vector<perm_bit<quspin::basis::uint32_t>> loc_symm;
                    std::vector<npy_cdouble_wrapper> lat_char((npy_cdouble_wrapper*)_lat_char, ((npy_cdouble_wrapper*)_lat_char) + lat_symm.size());
                    std::vector<npy_cdouble_wrapper> loc_char((npy_cdouble_wrapper*)_loc_char, ((npy_cdouble_wrapper*)_loc_char) + loc_symm.size());
                    for(std::shared_ptr<void> _lat_arg : _lat_args){
                        std::shared_ptr<bit_perm_args> lat_arg =  std::reinterpret_pointer_cast<bit_perm_args>(_lat_arg);
                        lat_symm.emplace_back(lat_arg->perm);
                    }
                    for(std::shared_ptr<void> _loc_arg : _loc_args){
                        std::shared_ptr<perm_bit_args> loc_arg =  std::reinterpret_pointer_cast<perm_bit_args>(_loc_arg);
                        loc_symm.emplace_back(loc_arg->mask);
                    }
                    std::shared_ptr<bit_symmetry<quspin::basis::uint32_t>> symmetry = std::make_shared<bit_symmetry<quspin::basis::uint32_t>>(lat_symm,lat_char,loc_symm,loc_char);
                    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);
                    break;
                }
            case 1:
                {
                    std::vector<bit_perm<quspin::basis::uint64_t>> lat_symm;
                    std::vector<perm_bit<quspin::basis::uint64_t>> loc_symm;
                    std::vector<npy_cdouble_wrapper> lat_char((npy_cdouble_wrapper*)_lat_char, ((npy_cdouble_wrapper*)_lat_char) + lat_symm.size());
                    std::vector<npy_cdouble_wrapper> loc_char((npy_cdouble_wrapper*)_loc_char, ((npy_cdouble_wrapper*)_loc_char) + loc_symm.size());
                    for(std::shared_ptr<void> _lat_arg : _lat_args){
                        std::shared_ptr<bit_perm_args> lat_arg =  std::reinterpret_pointer_cast<bit_perm_args>(_lat_arg);
                        lat_symm.emplace_back(lat_arg->perm);
                    }
                    for(std::shared_ptr<void> _loc_arg : _loc_args){
                        std::shared_ptr<perm_bit_args> loc_arg =  std::reinterpret_pointer_cast<perm_bit_args>(_loc_arg);
                        loc_symm.emplace_back(loc_arg->mask);
                    }
                    std::shared_ptr<bit_symmetry<quspin::basis::uint64_t>> symmetry = std::make_shared<bit_symmetry<quspin::basis::uint64_t>>(lat_symm,lat_char,loc_symm,loc_char);
                    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);
                    break;
                }
            case 2:
                {
                    std::vector<dit_perm<quspin::basis::uint32_t>> lat_symm;
                    std::vector<perm_dit<quspin::basis::uint32_t>> loc_symm;
                    std::vector<npy_cdouble_wrapper> lat_char((npy_cdouble_wrapper*)_lat_char, ((npy_cdouble_wrapper*)_lat_char) + lat_symm.size());
                    std::vector<npy_cdouble_wrapper> loc_char((npy_cdouble_wrapper*)_loc_char, ((npy_cdouble_wrapper*)_loc_char) + loc_symm.size());
                    for(std::shared_ptr<void> _lat_arg : _lat_args){
                        std::shared_ptr<dit_perm_args> lat_arg =  std::reinterpret_pointer_cast<dit_perm_args>(_lat_arg);
                        lat_symm.emplace_back(lhss,lat_arg->perm);
                    }
                    for(std::shared_ptr<void> _loc_arg : _loc_args){
                        std::shared_ptr<perm_dit_args> loc_arg = std::reinterpret_pointer_cast<perm_dit_args>(_loc_arg);
                        loc_symm.emplace_back(lhss,loc_arg->perms,loc_arg->locs);
                    }
                    std::shared_ptr<dit_symmetry<quspin::basis::uint32_t>> symmetry = std::make_shared<dit_symmetry<quspin::basis::uint32_t>>(lat_symm,lat_char,loc_symm,loc_char);
                    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);
                    break;
                }
            case 3:
                {
                    std::vector<dit_perm<quspin::basis::uint64_t>> lat_symm;
                    std::vector<perm_dit<quspin::basis::uint64_t>> loc_symm;
                    std::vector<npy_cdouble_wrapper> lat_char((npy_cdouble_wrapper*)_lat_char, ((npy_cdouble_wrapper*)_lat_char) + lat_symm.size());
                    std::vector<npy_cdouble_wrapper> loc_char((npy_cdouble_wrapper*)_loc_char, ((npy_cdouble_wrapper*)_loc_char) + loc_symm.size());
                    for(std::shared_ptr<void> _lat_arg : _lat_args){
                        std::shared_ptr<dit_perm_args> lat_arg =  std::reinterpret_pointer_cast<dit_perm_args>(_lat_arg);
                        lat_symm.emplace_back(lhss,lat_arg->perm);
                    }
                    for(std::shared_ptr<void> _loc_arg : _loc_args){
                        std::shared_ptr<perm_dit_args> loc_arg = std::reinterpret_pointer_cast<perm_dit_args>(_loc_arg);
                        loc_symm.emplace_back(lhss,loc_arg->perms,loc_arg->locs);
                    }
                    std::shared_ptr<dit_symmetry<quspin::basis::uint64_t>> symmetry = std::make_shared<dit_symmetry<quspin::basis::uint64_t>>(lat_symm,lat_char,loc_symm,loc_char);
                    symmetry_ptr = std::reinterpret_pointer_cast<void>(symmetry);
                    break;
                }
            default:
                throw std::runtime_error("cannot parse arguments.");
        
        }
    }

    ~symmetry_abi(){}

    std::shared_ptr<void> data()
    {
        return symmetry_ptr;
    }
    void* get()
    {
        return symmetry_ptr.get();
    }

};

    
    }
    #endif