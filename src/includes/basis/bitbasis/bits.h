#ifndef __QUSPIN_BASIS_BITBASIS_BITS_H__
#define __QUSPIN_BASIS_BITBASIS_BITS_H__

#include "basis/bitbasis/types.h"

namespace quspin::basis::bit_basis {

// local degrees of freedom stored in contiguous chunks of bits
template<typename I>
typedef struct bit_set { // thin wrapper used for convience

    I content;

    bit_set(I _content) : content(_content) {}
    bit_set(const bit_set<I>& other) : content(other.content) {}

} bit_set;

template<typename I>
int get_sub_bitstring(bit_set<I> s,const int i){
    return integer_cast<int,I>( (s.content >> i) & I(1))
}

template<typename I>
int get_sub_bitstring(const bit_set<I> s,const int * locs,const int nlocs){
    int out = 0;
    for(int i=0;i<nlocs;++i){
        out |= get_sub_bitstring(s,locs[i]);
        out <<= 1;
    }
    return out;
}

template<typename I>
bit_set<I> set_sub_bitstring(const bit_set<I> s,const I in,const int i){
    return  bit_set<I>(s.content ^ (in << i));
}

template<typename I>
bit_set<I> set_sub_bitstring(const bit_set<I> s,int in,const int * locs,const int nlocs){
    bit_set<I> r(s);
    for(int i=0;i<nlocs;++i){
        r = set_sub_bitstring(r,  I(in & 1), locs[i]);
        in >>= 1;
    }
    return  r;
}


// permutations

template <typename I>
class bit_perm // permutation of dit locations
{
private:
    const int * perm;
    const int length;
    benes_perm::tr_benes<I> benes;

public:
    bit_perm(const int * _perm,const int _length) : 
    lhss(_lhss), perm(_perm), length(_length)
    {
        benes_perm::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<length;i++){
            index[i] = perm[i];  
        }

        benes_perm::gen_benes(&benes,index);

    }
    ~bit_perm() {}

    inline bit_set<I> app(const bit_set<I> s) const { 
        return bit_set<I>(benes_bwd(benes,s.content)); 
    }

    inline bit_set<I> inv(const bit_set<I> s) const { 
        return bit_set<I>(benes_fwd(benes,s.content)); 
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

    inline bit_set<I> app(const bit_set<I> s) const {

        return  bit_set<I>( s.content^mask );
    }

    inline I inv(const I s) const {
        return bit_set<I>( s.content^mask );
    }
};


}

#endif