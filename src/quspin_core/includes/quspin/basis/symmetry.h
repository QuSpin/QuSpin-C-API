#ifndef __QUSPIN_BASIS_SYMMETRY_H__
#define __QUSPIN_BASIS_SYMMETRY_H__


#include <vector>
#include <memory>
#include "quspin/basis/types.h"
#include "quspin/basis/bitbasis/dits.h"
#include "quspin/basis/bitbasis/bits.h"
#include "quspin/basis/bitbasis/benes.h"


namespace quspin::basis {

template <typename I>
class dit_perm // permutation of dit locations
{
private:
    const int lhss;
    const int length;
    benes::tr_benes<I> benes;

public:
    dit_perm() : lhss(0), length(0)  {}// identity constructor
    dit_perm(const int _lhss,const int * _perm,const int _length) : 
    lhss(_lhss), length(_length)
    {
        benes::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<_length;i++){
            const int dst = _perm[i];
            const int src = i;
            for(int j=0;j<constants::bits[lhss];j++){
                const int srcbit = constants::bits[_lhss] * src + j;
                const int dstbit = constants::bits[_lhss] * dst + j;
                index[srcbit] = dstbit;
            }   
        }

        benes::gen_benes(&benes,index);

    }
    ~dit_perm() {}

    template<typename T>
    inline const dit_set<I> app(const dit_set<I>& s, T& coeff) const { 
        return (length != 0 ? dit_set<I>(benes::benes_bwd(benes,s.content),s.lhss,s.mask,s.bits) : s); 
    }

    template<typename T>
    inline const dit_set<I> inv(const dit_set<I>& s, T& coeff) const { 
        return (length != 0 ? dit_set<I>(benes::benes_fwd(benes,s.content),s.lhss,s.mask,s.bits) : s);  
    }
};

template <typename I>
class perm_dit // permutations of the dit states locally
{
private:
    const int lhss;
    const int length;
    std::vector<dit_integer_t> perm;
    std::vector<dit_integer_t> inv_perm;
    std::vector<int> locs;

public:
    perm_dit() : lhss(0), length(0) {} // identity
    perm_dit(
        const int _lhss, const dit_integer_t * _perm, const int * _locs, const int _length
    ) : lhss(_lhss), length(_length)
    { 
        perm.resize(_length * _lhss);
        inv_perm.resize(_length * _lhss);
        locs.resize(_length);

        dit_integer_t * loc_perm = perm.data();
        dit_integer_t * loc_inv_perm = inv_perm.data();

        for(int i=0;i<_length;++i){
            locs[i] = _locs[i];
            for(int j=0;j<_lhss;++j){
                loc_perm[j] = _perm[j];
                loc_inv_perm[_perm[j]] = j;
            }
            loc_perm += _lhss;
            loc_inv_perm += _lhss;
            _perm += _lhss;
        }

    }
    ~perm_dit() {}

    template<typename T>
    dit_set<I> app(const dit_set<I>& s, T& coeff) const {
        
        dit_set<I> r(s);
        
        dit_integer_t * perm_ptr = perm.data();

        for(int i=0;i<length;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,perm_ptr[bits],locs[i]);
            perm_ptr += lhss;
        }

        return r;
    }

    template<typename T>
    dit_set<I> inv(const dit_set<I>& s, T& coeff) const {
        dit_set<I> r(s);
        
        dit_integer_t * perm_ptr = inv_perm.data();

        for(int i=0;i<length;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,perm_ptr[bits],locs[i]);
            perm_ptr += lhss;
        }
        
        return r;
    }
};

template <typename I>
class bit_perm // permutation of dit locations
{
private:
    const int length;
    benes::tr_benes<I> benes;

public:
    bit_perm() : length(0) {} // identity constructor
    bit_perm(const int * _perm,const int _length) : length(_length)
    {
        benes::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<_length;i++){
            index[i] = _perm[i];  
        }

        benes::gen_benes(&benes,index);

    }
    ~bit_perm() {}

    template<typename T>
    inline const bit_set<I> app(const bit_set<I>& s, T& coeff) const { 
        return (length != 0 ? bit_set<I>(benes::benes_bwd(benes,s.content)) : s); 
    }

    template<typename T>
    inline const bit_set<I> inv(const bit_set<I>& s, T& coeff) const { 
        return (length != 0 ?bit_set<I>(benes::benes_fwd(benes,s.content)) : s); 
    }
};

template <typename I>
class perm_bit // permutations of the dit states locally
{
private:
    const I mask; // which bits to flip

public:
    perm_bit() : mask(0) {} // identity (don't flip any bits)
    perm_bit(const I _mask) : mask(_mask) { }
    ~perm_bit() {}

    template<typename T>
    inline bit_set<I> app(const bit_set<I> s, T& coeff) const {

        return  bit_set<I>( s.content^mask );
    }

    template<typename T>
    inline bit_set<I> inv(const I s, T& coeff) const {
        return bit_set<I>( s.content^mask );
    }
};


template<typename lat_perm_t,typename loc_perm_t,typename dits_or_bits,typename T>
class symmetry
{
private:
    std::vector<lat_perm_t> lat_symm;
    std::vector<loc_perm_t> loc_symm;
    std::vector<T> lat_char;
    std::vector<T> loc_char;



public:
    symmetry(std::vector<lat_perm_t> &_lat_symm,std::vector<T> &_lat_char,
    std::vector<loc_perm_t> &_loc_symm,std::vector<T> &_loc_char)
    {
        assert((lat_symm.size() == lat_char.size()) && (loc_symm.size() == loc_char.size()));
        lat_symm = _lat_symm;
        loc_symm = _loc_symm;
        lat_char = _lat_char;
        loc_char = _loc_char;
    }

    symmetry(symmetry<lat_perm_t,loc_perm_t,dits_or_bits,T>& other){
        lat_symm = other.lat_symm;
        loc_symm = other.loc_symm;
        lat_char = other.lat_char;
        loc_char = other.loc_char;
    }
    ~symmetry() {}

    size_t size() const {return lat_symm.size() * loc_symm.size();}

    std::pair<dits_or_bits,T> get_refstate(const dits_or_bits &s) const {

        dits_or_bits ss(s);
        T sign = T(1);
        T coeff = T(1);
        
        for(int i=0;i<loc_symm.size();++i) 
        {
            const auto r = loc_symm[i].app(s,sign);
            for(int j=0;j<lat_symm.size();++j)
            {
                const auto rr = lat_symm[i].app(r,sign);

                if(rr > ss){
                    ss = rr;
                    coeff = sign * lat_char[j] * loc_char[i];
                }
            }
        }

        return std::make_pair(ss,coeff);
    }

    double calc_norm(const dits_or_bits &s) const {
        double norm = 0.0;
        T sign = T(1);

        for(int i=0;i<loc_symm.size();++i) 
        {
            const auto r = loc_symm[i].app(s,sign);
            for(int j=0;j<lat_symm.size();++j)
            {
                const auto rr = lat_symm[i].app(r,sign);
                if(rr == s) norm += real_value(sign * lat_char[j] * loc_char[i]);
            }
        }
        return norm;
    }

    double check_refstate(const dits_or_bits &s) const {
        double norm=0.0;
        T sign = T(1);

        for(int i=0;i<loc_symm.size();++i) 
        {
            const auto r = loc_symm[i].app(s,sign);
            for(int j=0;j<lat_symm.size();++j)
            {
                const auto rr = lat_symm[i].app(r,sign);

                if(rr >  s){return std::numeric_limits<double>::quiet_NaN();};
                if(rr == s) norm += real_value(sign * lat_char[j] * loc_char[i]);
            }
        }

        return norm;
    }

};

}


#ifdef QUSPIN_UNIT_TESTS


namespace quspin::basis {

template class dit_perm<uint8_t>;
template class perm_dit<uint8_t>;
template class bit_perm<uint8_t>;
template class perm_bit<uint8_t>;
template class symmetry<bit_perm<uint8_t>,perm_bit<uint8_t>,bit_set<uint8_t>,double>;

}



#endif

#endif