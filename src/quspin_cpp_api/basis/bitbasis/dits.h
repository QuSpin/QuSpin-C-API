#ifndef __QUSPIN_BASIS_BITBASIS_DITS_H__
#define __QUSPIN_BASIS_BITBASIS_DITS_H__

#include "basis/bitbasis/types.h"

namespace quspin::basis {

namespace constants {

static const quspin::basis::types::dit_integer_t bits[255] = {
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

static const quspin::basis::types::dit_integer_t mask[255] = {
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
    
    typedef I bitset_t;

    I content;
    const int lhss;
    const I mask;
    const types::dit_integer_t bits;

    dit_set(I _content,const int _lhss, const I _mask, const dit_integer_t bits) : 
    content(_content), 
    lhss(_lhss), 
    mask(_mask), 
    bits(bits) {}

    dit_set(const dit_set<I> other) : 
    content(other.content), 
    lhss(other.lhss), 
    mask(other.mask), 
    bits(other.bits) {}


} dit_set;

template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int i){
    return integer<int,I>::cast( (s.content >> (i * s.bits)) & s.mask)
}

template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int * locs,const int nlocs){
    int out = 0;
    for(int i=0;i<nlocs;++i){
        out += get_sub_bitstring(s,locs[i])
        out *= s.lhss;
        /* // implementation with padding for congituous blocks 
        out |= bit_basis<int,I>(get_sub_bitstring(s,locs[i]));
        out <<= s.bits;
        */
    }
    return out;
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,const I in,const int i){
    const int shift =  i * s.bits;
    const I r = s.content ^ ( ( in << shift ) ^ s.content) & (s.mask << shift)
    return  dit_set<I>(r,s.lhss,s.mask,s.bits);
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,int in,const int * locs,const int nlocs){
    dit_set<I> r(s);
    for(int i=0;i<nlocs;++i){
        r = set_sub_bitstring(r,  I(in % s.lhss), locs[i]);
        in /= s.lhss;
        /* // implementation with padding for congituous blocks 
        r = set_sub_bitstring(r,  in & s.mask, locs[i]);
        in >>= s.bits;
        */
    }

    return  r;
}

} // end namespace quspin::basis


#ifdef QUSPIN_UNIT_TESTS

TEST_CASE("get_bit_substring") {

}

TEST_CASE("set_sub_bitstring") {

}

TEST_CASE("operators") {

}

TEST_CASE("to_/from_vector") {
    
}

#endif


#endif