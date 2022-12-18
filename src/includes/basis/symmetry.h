#ifndef __QUSPIN_BASIS_SYMMETRY_H__
#define __QUSPIN_BASIS_SYMMETRY_H__


#include <vector>

#include "basis/bitbasis/dits.h"
#include "basis/bitbasis/bits.h"
#include "basis/bitbasis/types.h"

namespace quspin::basis {

template <typename I>
class dit_perm // permutation of dit locations
{
private:
    const int lhss;
    const int length;
    bit_basis::benes_perm::tr_benes<I> benes;

public:
    dit_perm() : lhss(0), length(0) // identity constructor
    dit_perm(const int _lhss,const int * _perm,const int _length) : 
    lhss(_lhss), length(_length)
    {
        bit_basis::benes_perm::ta_index<I> index;
        for(int i=0;i<bit_basis::bit_info<I>::bits;i++){index[i] = bit_basis::benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<_length;i++){
            const int dst = _perm[i];
            const int src = i;
            for(int j=0;j<bits[lhss];j++){
                const int srcbit = bit_basis::constants::bits[_lhss] * src + j;
                const int dstbit = bit_basis::constants::bits[_lhss] * dst + j;
                index[srcbit] = dstbit;
            }   
        }

        bit_basis::benes_perm::gen_benes(&benes,index);

    }
    ~dit_perm() {}

    template<typename T>
    inline const bit_basis::dit_set<I> app(const bit_basis::dit_set<I>& s, T& coeff) const { 
        return (length != 0 ? bit_basis::dit_set<I>(benes_bwd(benes,s.content),s.lhss,s.mask,s.bits) : s); 
    }

    template<typename T>
    inline const bit_basis::dit_set<I> inv(const bit_basis::dit_set<I>& s, T& coeff) const { 
        return (length != 0 ? bit_basis::dit_set<I>(benes_fwd(benes,s.content),s.lhss,s.mask,s.bits) : s);  
    }
};

template <typename I>
class perm_dit<I> // permutations of the dit states locally
{
private:
    const int lhss;
    const int length;
    unique_ptr<bit_basis::types::dit_integer_t> * perm;
    unique_ptr<bit_basis::types::dit_integer_t> * inv_perm;
    unique_ptr<int> * locs;

public:
    perm_dit() : lhss(0), perm(NULL), inv_perm(NULL), locs(NULL), length(0) // identity
    perm_dit(
        const int _lhss, const dit_integer_t * _perm, const int * _locs, const int _length
    ) : lhss(_lhss), length(_length)
    { 
        perm(new bit_basis::types::dit_integer_t[_length * _lhss]);
        inv_perm(new bit_basis::types::dit_integer_t[_length * _lhss]);
        locs(new int[_length]);

        int * loc_perm = perm.get();
        int * loc_inv_perm = inv_perm.get();

        for(int i=0;i<_length;++i){
            locs[i] = _locs[i];
            for(int j=0;j<_lhss;++j){
                loc_perm[j] = _perm[j];
                loc_inv_perm[_perm[j]] = j;
            }
            local_perm += _lhss;
            local_inv_perm += _lhss;
            _perm += _lhss;
        }

    }
    ~dit_perm() {}

    template<typename T>
    bit_basis::dit_set<I> app(const bit_basis::dit_set<I>& s, T& coeff) const {
        
        bit_basis::dit_set<I> r(s);
        
        const dit_integer_t * perm_ptr = perm.get();

        for(int i=0;i<length;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,perm_ptr[bits],locs[i]);
            perm_ptr += lhss;
        }

        return r;
    }

    template<typename T>
    bit_basis::dit_set<I> inv(const bit_basis::dit_set<I>& s, T& coeff) const {
        bit_basis::dit_set<I> r(s);
        
        const dit_integer_t * perm_ptr = perm_inv.get();

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
    bit_basis::benes_perm::tr_benes<I> benes;

public:
    dit_perm() : length(0) // identity constructor
    bit_perm(const int * _perm,const int _length) : 
    lhss(_lhss), length(_length)
    {
        bit_basis::benes_perm::ta_index<I> index;
        for(int i=0;i<bit_basis::bit_info<I>::bits;i++){index[i] = bit_basis::benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<_length;i++){
            index[i] = _perm[i];  
        }

        bit_basis::benes_perm::gen_benes(&benes,index);

    }
    ~bit_perm() {}

    template<typename T>
    inline const bit_basis::bit_set<I> app(const bit_basis::bit_set<I>& s, T& coeff) const { 
        return (length != 0 ? bit_basis::bit_set<I>(benes_bwd(benes,s.content)) : s); 
    }

    template<typename T>
    inline const bit_basis::bit_set<I> inv(const bit_basis::bit_set<I>& s, T& coeff) const { 
        return (length != 0 ?bit_basis::bit_set<I>(benes_fwd(benes,s.content)) : s); 
    }
};

template <typename I>
class perm_bit // permutations of the dit states locally
{
private:
    const I mask; // which bits to flip

public:
    perm_bit() : mask(0) // identity (don't flip any bits)
    perm_bit(const I _mask) : mask(_mask) { }
    ~bit_perm() {}

    template<typename T>
    inline bit_basis::bit_set<I> app(const bit_basis::bit_set<I> s, T& coeff) const {

        return  bit_basis::bit_set<I>( s.content^mask );
    }

    template<typename T>
    inline bit_basis::bit_set<I> inv(const I s, T& coeff) const {
        return bit_basis::bit_set<I>( s.content^mask );
    }
};


template<typename lat_perm_t,typename loc_perm_t,typename dits_or_bits,typename T>
class symmetry
{
private:
    std::vector<lat_perm_t> lat_symm;
    std::vector<loc_perm_t> loc_symm;
    std::vector<T> lat_chars;
    std::vector<T> loc_chars;



public:
    lattice_symmetry(std::vector<lat_perm_t> &_lat_symm,std::vector<T> &_lat_chars,
    std::vector<loc_perm_t> &_loc_symm,std::vector<T> &_loc_chars) : 
    {
        lat_symm = _lat_symm;
        loc_symm = _loc_symm;
        lat_char = _lat_char;
        loc_char = _loc_char;

        assert((lat_symm.size() == lat_char.size()) && (loc_symm.size() == loc_char.size()));

    }

    lattice_symmetry(lattice_symmetry<lat_perm_t,loc_perm_t,dits_or_bits,T>& other){
        lat_symm = other.lat_symm;
        loc_symm = other.loc_symm;
        lat_char = other.lat_char;
        loc_char = other.loc_char;
    }
    ~lattice_symmetry() {}

    size_t size() const {return lattice_perms.size() * local_perms.size();}

    std::pair<dits_or_bits,T> get_refstate(const dits_or_bits &s){

        dits_or_bits ss(s);
        T sign = T(1);
        
        for(int i=0;i<loc_symm.size();++i) 
        {
            const auto r = loc_symm[i].app(s,sign);
            for(int j=0;j<lat_symm.size();++j)
            {
                const auto rr = lat_symm[i].app(r,sign);

                if(rr > ss){
                    ss = rr;
                    coeff = sign * lat_chars[p] * loc_chars[i];
                }
            }
        }

        return make_pair(ss,coeff);
    }

    double check_refstate(const dits_or_bits &s){
        double norm=0.0;
        

        for(int i=0;i<loc_symm.size();++i) 
        {
            const auto r = loc_symm[i].app(s,norm);
            for(int j=0;j<lat_symm.size();++j)
            {
                const auto rr = lat_symm[i].app(r,norm);

                if(rr >  s){return std::numeric_limits<double>::quiet_NaN()};
                if(rr == s) norm += std::real(lat_char[j] * loc_char[i]);
            }
        }

        return norm;
    }

};

}

#endif