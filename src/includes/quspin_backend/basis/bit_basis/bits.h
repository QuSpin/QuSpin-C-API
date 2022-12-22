#ifndef __QUSPIN_BASIS_BITBASIS_BITS_H__
#define __QUSPIN_BASIS_BITBASIS_BITS_H__

#include "quspin_backend/basis/bit_basis/types.h"
#include "quspin_backend/basis/bit_basis/info.h"

namespace quspin::basis::bit_basis {

// local degrees of freedom stored in contiguous chunks of bits
template<typename I>
typedef struct bit_set { // thin wrapper used for convience

    // doesn't allocate any data
    static const int lhss = 2;
    static const types::dit_integer_t mask = 1;
    static const types::dit_integer_t bits = 1;

    typedef I bitset_t;
    I content;

    bit_set(I _content) : content(_content) {}
    bit_set(const bit_set<I>& other) : content(other.content) {}

};

template<typename I>
int get_sub_bitstring(bit_set<I> s,const int i){
    return integer_cast<int,I>( (s.content >> i) & I(1));
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

template<typename I>
inline bool operator<(const bit_set<I>& lhs, const bit_set<I>& rhs){return lhs.content < rhs.content;}

template<typename I>
inline bool operator>(const bit_set<I>& lhs, const bit_set<I>& rhs){return lhs.content > rhs.content;}

template<typename I>
inline bool operator==(const bit_set<I>& lhs, const bit_set<I>& rhs){return lhs.content == rhs.content;}

template<typename I>
std::vector<bool> to_vector(const bit_set<I>& s,const int length=0){
    std::vector<bool> out;

    const int niter = (length>0 ? length : bit_info<I>::bits/s.bits);

    for(int i=0;i<niter;++i){
        out.push_back(static_cast<bool>(get_sub_bitstring(s,i)));
    }
    return out;
}

template<typename I>
bit_set<I> from_vector(const std::vector<bool>& input){
    bit_set<I>  state;

    for(int i=0;i<input.size();i++){
        set_sub_bitstring(state,static_cast<int>(input[i]),i);
    }

    return state;
}


} // end namespace quspin::basis


#ifdef __UNIT_TESTS__

TEST_CASE("get_bit_substring") {

}

TEST_CASE("set_sub_bitstring") {

}

TEST_CAST("operators") {

}

TEST_CAST("to_/from_vector") {
    
}

#endif



#endif