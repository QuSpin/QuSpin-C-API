#ifndef __QUSPIN_BASIS_SYMMETRY_H__
#define __QUSPIN_BASIS_SYMMETRY_H__


#include <vector>
#include <algorithm>
#include <cassert>
#include "quspin/basis/types.h"
#include "quspin/basis/bitbasis/utils.h"
#include "quspin/basis/bitbasis/dits.h"
#include "quspin/basis/bitbasis/bits.h"
#include "quspin/basis/bitbasis/benes.h"
#include "quspin/utils/functions.h"


namespace quspin::basis {

template <typename I>
class dit_perm // permutation of dit locations
{
private:
    benes::tr_benes<I> benes;

public:
    dit_perm(const int lhss,const int * perm,const int length) {
        benes::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes::no_index;}

        
        // number of bits to store lhss worth of data:
        const dit_integer_t bits = constants::bits[lhss];

        // permute chucks of bits of length 'bits' 
        for(int i=0;i<length;i++){
            const int dst = perm[i];
            const int src = i;
            for(int j=0;j<bits;j++){
                const int srcbit = bits * src + j;
                const int dstbit = bits * dst + j;
                index[srcbit] = dstbit;
            }   
        }

        benes::gen_benes(&benes,index);

    }

    ~dit_perm() {}

    template<typename T>
    inline dit_set<I> app(const dit_set<I>& s, T& coeff) const { 
        return dit_set<I>(benes::benes_bwd(&benes,s.content),s.lhss,s.mask,s.bits); 
    }

    template<typename T>
    inline dit_set<I> inv(const dit_set<I>& s, T& coeff) const { 
        return dit_set<I>(benes::benes_fwd(&benes,s.content),s.lhss,s.mask,s.bits);  
    }
};

template <typename I>
class perm_dit // permutations of the dit states locally
{
private:
    int lhss;
    std::vector<int> perm;
    std::vector<int> inv_perm;
    std::vector<int> locs;

public:
    perm_dit(const int _lhss,const std::vector<std::vector<int>>& _perm,const std::vector<int> _locs) : lhss(_lhss)
    { 

        assert(_perm.size() == _locs.size());

        locs.insert(locs.end(),_locs.begin(),_locs.end());

        for(int i=0;i<locs.size();++i){
            const std::vector<int> p = _perm.at(i);
            std::vector<int> ip(p.size());
            int j = 0;
            for(const int ele : p){ip[ele] = j++;}
            perm.insert(perm.end(),p.begin(),p.end());
            inv_perm.insert(inv_perm.end(),ip.begin(),ip.end());
        }

    }
    ~perm_dit() {}

    template<typename T>
    dit_set<I> app(const dit_set<I>& s, T& coeff) const {
        
        dit_set<I> r(s);
        
        const int * perm_ptr = perm.data();

        for(const int loc : locs){
            const int bits = get_sub_bitstring(s,loc);
            r = set_sub_bitstring(r,perm_ptr[bits],loc);
            perm_ptr += lhss;
        }

        return r;
    }

    template<typename T>
    dit_set<I> inv(const dit_set<I>& s, T& coeff) const {
        dit_set<I> r(s);
        
        const int * perm_ptr = inv_perm.data();

        for(const int loc : locs){
            const int bits = get_sub_bitstring(s,loc);
            r = set_sub_bitstring(r,perm_ptr[bits],loc);
            perm_ptr += lhss;
        }
        
        return r;
    }
};

template <typename I>
class bit_perm // permutation of bit locations
{
private:
    benes::tr_benes<I> benes;

public:
    bit_perm(const int * perm,const int length) 
    {
        benes::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<length;i++){
            index[i] = perm[i];  
        }

        benes::gen_benes(&benes,index);

    }
    ~bit_perm() {}

    template<typename T>
    inline bit_set<I> app(const bit_set<I>& s, T& coeff) const { 
        return  bit_set<I>(benes::benes_bwd(&benes,s.content)); 
    }

    template<typename T>
    inline bit_set<I> inv(const bit_set<I>& s, T& coeff) const { 
        return bit_set<I>(benes::benes_fwd(&benes,s.content)); 
    }
};

template <typename I>
class perm_bit // permutations of the bit states locally
{
private:
    I mask; // which bits to flip

public:
    perm_bit(const I _mask) : mask(_mask) { }
    ~perm_bit() {}

    template<typename T>
    inline bit_set<I> app(const bit_set<I>& s, T& coeff) const {

        return  bit_set<I>( s.content^mask );
    }

    template<typename T>
    inline bit_set<I> inv(const bit_set<I>& s, T& coeff) const {
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
        lat_symm.insert(lat_symm.end(),_lat_symm.begin(),_lat_symm.end());
        loc_symm.insert(loc_symm.end(),_loc_symm.begin(),_loc_symm.end());
        lat_char.insert(lat_char.end(),_lat_char.begin(),_lat_char.end());
        loc_char.insert(loc_char.end(),_loc_char.begin(),_loc_char.end());
    }

    symmetry(symmetry<lat_perm_t,loc_perm_t,dits_or_bits,T>& other){
        lat_symm.insert(lat_symm.end(),other.lat_symm.begin(),other.lat_symm.end());
        loc_symm.insert(loc_symm.end(),other.loc_symm.begin(),other.loc_symm.end());
        lat_char.insert(lat_char.end(),other.lat_char.begin(),other.lat_char.end());
        loc_char.insert(loc_char.end(),other.loc_char.begin(),other.loc_char.end());
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
                if(rr == s) norm += real(sign * lat_char[j] * loc_char[i]);
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
                if(rr == s) norm += real(sign * lat_char[j] * loc_char[i]);
            }
        }

        return norm;
    }

};

}



#ifdef QUSPIN_UNIT_TESTS


namespace quspin::basis { // test cases

template class dit_perm<uint8_t>;
template dit_set<uint8_t> dit_perm<uint8_t>::app<double>(const dit_set<uint8_t>&, double&) const;
template dit_set<uint8_t> dit_perm<uint8_t>::inv<double>(const dit_set<uint8_t>&, double&) const;

template class perm_dit<uint8_t>;
template dit_set<uint8_t> perm_dit<uint8_t>::app<double>(const dit_set<uint8_t>&, double&) const;
template dit_set<uint8_t> perm_dit<uint8_t>::inv<double>(const dit_set<uint8_t>&, double&) const;

template class bit_perm<uint8_t>;
template bit_set<uint8_t> bit_perm<uint8_t>::app<double>(const bit_set<uint8_t>&, double&) const;
template bit_set<uint8_t> bit_perm<uint8_t>::inv<double>(const bit_set<uint8_t>&, double&) const;

template class perm_bit<uint8_t>;
template bit_set<uint8_t> perm_bit<uint8_t>::app<double>(const bit_set<uint8_t>&, double&) const;
template bit_set<uint8_t> perm_bit<uint8_t>::inv<double>(const bit_set<uint8_t>&, double&) const;

template class symmetry<bit_perm<uint8_t>,perm_bit<uint8_t>,bit_set<uint8_t>,double>;


}

#endif

#endif