#ifndef __QUSPIN_BASIS_SYMMETRY_H__
#define __QUSPIN_BASIS_SYMMETRY_H__


#include <vector>

#include "basis/bitbasis/dits.h"
#include "basis/bitbasis/bits.h"

namespace quspin::basis {

template <typename I>
class dit_perm // permutation of dit locations
{
private:
    const int lhss;
    const int * perm;
    const int length;
    bit_basis::benes_perm::tr_benes<I> benes;

public:
    dit_perm(const int _lhss,const int * _perm,const int _length) : 
    lhss(_lhss), perm(_perm), length(_length)
    {
        bit_basis::benes_perm::ta_index<I> index;
        for(int i=0;i<bit_basis::bit_info<I>::bits;i++){index[i] = bit_basis::benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<length;i++){
            const int dst = perm[i];
            const int src = i;
            for(int j=0;j<bits[lhss];j++){
                const int srcbit = bits[lhss] * src + j;
                const int dstbit = bits[lhss] * dst + j;
                index[srcbit] = dstbit;
            }   
        }

        bit_basis::benes_perm::gen_benes(&benes,index);

    }
    ~dit_perm() {}

    inline bit_basis::dit_set<I> app(const bit_basis::dit_set<I>& s) const { 
        return bit_basis::dit_set<I>(benes_bwd(benes,s.content),s.lhss,s.mask,s.bits); 
    }

    inline bit_basis::dit_set<I> inv(const bit_basis::dit_set<I>& s) const { 
        return bit_basis::dit_set<I>(benes_fwd(benes,s.content),s.lhss,s.mask,s.bits); 
    }
};

template <typename I>
class perm_dit<I> // permutations of the dit states locally
{
private:
    const int lhss;
    const dit_integer_t * perm;
    const dit_integer_t * inv_perm;
    const int * locs;
    const int nlocs;
    const int length;

public:
    perm_dit(const int _lhss,const int * _perm,const int _length) : perm(_perm), length(_length), lhss(_lhss) { }
    ~dit_perm() {}

    bit_basis::dit_set<I> app(const bit_basis::dit_set<I>& s) const {
        bit_basis::dit_set<I> r(s);
        for(int i=0;i<nlocs;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,perm[bits],locs[i]);
        }
        return r;
    }
   bit_basis::dit_set<I> inv(const bit_basis::dit_set<I>& s) const {
        bit_basis::dit_set<I> r(s);
        for(int i=0;i<nlocs;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,inv_perm[bits],locs[i]);
        }
        return r;
    }
};

template <typename I>
class bit_perm // permutation of dit locations
{
private:
    const int * perm;
    const int length;
    bit_basis::benes_perm::tr_benes<I> benes;

public:
    bit_perm(const int * _perm,const int _length) : 
    lhss(_lhss), perm(_perm), length(_length)
    {
        bit_basis::benes_perm::ta_index<I> index;
        for(int i=0;i<bit_basis::bit_info<I>::bits;i++){index[i] = bit_basis::benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<length;i++){
            index[i] = perm[i];  
        }

        bit_basis::benes_perm::gen_benes(&benes,index);

    }
    ~bit_perm() {}

    inline bit_basis::bit_set<I> app(const bit_basis::bit_set<I>& s) const { 
        return bit_basis::bit_set<I>(benes_bwd(benes,s.content)); 
    }

    inline bit_basis::bit_set<I> inv(const bit_basis::bit_set<I>& s) const { 
        return bit_basis::bit_set<I>(benes_fwd(benes,s.content)); 
    }
};

template <typename I>
class perm_bit // permutations of the dit states locally
{
private:
    const I mask; // which bits to flip

public:
    perm_bit(const I _mask) : mask(_mask) { }
    ~bit_perm() {}

    inline bit_basis::bit_set<I> app(const bit_basis::bit_set<I> s) const {

        return  bit_basis::bit_set<I>( s.content^mask );
    }

    inline I inv(const I s) const {
        return bit_basis::bit_set<I>( s.content^mask );
    }
};


template<typename perm,typename T>
class symmetry
{
private:
    std::vector<perm> perms;
    std::vector<T> chars;

public:
    lattice_symmetry(std::vector<perm> &_perms,std::vector<T> &_chars) : perms(_perms), chars(_chars) {
        assert(_perms.size() == _chars.size());
    }
    ~lattice_symmetry() {}

    size_t size() const {return perms.size();}
    T character(const size_t i) const {return chars[i];}
    const perm& operator[](const size_t i) const {return perms[i];}

};

}

#endif