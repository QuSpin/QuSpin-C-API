#ifndef __QUSPIN_BASIS_BITBASIS_DITS_H__
#define __QUSPIN_BASIS_BITBASIS_DITS_H__

#include "basis/bitbasis/types.h"


namespace bit_basis {


typedef uint_fast8_t dit_integer_t; 

namespace constants {

static const dit_integer_t nbits[255] = {
    1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};

static const dit_integer_t mask[255] = {
    1,   1,   1,   3,   3,   7,   7,   7,   7,  15,  15,  15,  15,
    15,  15,  15,  15,  31,  31,  31,  31,  31,  31,  31,  31,  31,
    31,  31,  31,  31,  31,  31,  31,  63,  63,  63,  63,  63,  63,
    63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,
    63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,  63,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255
};

}

// local degrees of freedom stored in contiguous chunks of bits
template<typename I>
typedef struct dit_set { // thin wrapper used for convience

    I content;
    const int lhss;
    const I mask;
    const dit_integer_t nbits;

    dit_set(I _content,const int _lhss, const dit_integer_t nbits, const I _mask) : 
    content(_content), 
    lhss(_lhss), 
    mask(_mask), 
    nbits(nbits) {}

    dit_set(const dit_set<I> other) : 
    content(other.content), 
    lhss(other.lhss), 
    mask(other.mask), 
    nbits(other.nbits) {}


} dit_set;




template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int i){
    return integer_cast<int,I>( (s.content >> (i * s.nbits)) & s.mask)
}

template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int * locs,const int nlocs){
    int out = 0;
    for(int i=0;i<nlocs;++i){
        out += get_sub_bitstring(s,locs[i])
        out *= s.lhss;
        /* // implementation with padding for congituous blocks 
        out |= integer_cast<int,I>(get_sub_bitstring(s,locs[i]));
        out <<= s.nbits;
        */
    }
    return out;
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,const I in,const int i){
    const int shift =  i * s.nbits;
    const I r = s.content ^ ( -( in << shift ) ^ s.content) & (s.mask << shift)
    return  dit_set<I>(r,s.lhss,s.mask,s.nbits);
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,int in,const int * locs,const int nlocs){
    dit_set<I> r(s);
    for(int i=0;i<nlocs;++i){
        r = set_sub_bitstring(r,  I(in % s.lhss), locs[i]);
        in /= s.lhss;
        /* // implementation with padding for congituous blocks 
        r = set_sub_bitstring(r,  in & s.mask, locs[i]);
        in >>= s.nbits;
        */
    }

    return  r;
}

// permutations

template <typename I>
class dit_perm // permutation of dit locations
{
private:
    const int lhss;
    const int * perm;
    const int length;
    benes_perm::tr_benes<I> benes;

public:
    dit_perm(const int _lhss,const int * _perm,const int _length) : 
    lhss(_lhss), perm(_perm), length(_length)
    {
        benes_perm::ta_index<I> index;
        for(int i=0;i<bit_info<I>::bits;i++){index[i] = benes_perm::no_index;}

        // benes permutation is agnostic to lhss
        for(int i=0;i<length;i++){
            const int dst = perm[i];
            const int src = i;
            for(int j=0;j<nbits[lhss];j++){
                const int srcbit = nbits[lhss] * src + j;
                const int dstbit = nbits[lhss] * dst + j;
                index[srcbit] = dstbit;
            }   
        }

        benes_perm::gen_benes(&benes,index);

    }
    ~dit_perm() {}

    inline dit_set<I> app(const dit_set<I> s) const { 
        return dit_set<I>(benes_bwd(benes,s.content),s.lhss,s.mask,s.nbits); 
    }

    inline dit_set<I> inv(const dit_set<I> s) const { 
        return dit_set<I>(benes_fwd(benes,s.content),s.lhss,s.mask,s.nbits); 
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

    dit_set<I> app(const dit_set<I> s) const {
        dit_set<I> r(s);
        for(int i=0;i<nlocs;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,perm[bits],locs[i]);
        }
        return r;
    }
    dit_set<I> inv(const dit_set<I> s) const {
        dit_set<I> r(s);
        for(int i=0;i<nlocs;++i){
            const int bits = get_sub_bitstring(s,locs[i]);
            r = get_sub_bitstring(r,inv_perm[bits],locs[i]);
        }
        return r;
    }
};
}


#endif