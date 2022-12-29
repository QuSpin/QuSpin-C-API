#ifndef __QUSPIN_BASIS_BITBASIS_DITS_H__
#define __QUSPIN_BASIS_BITBASIS_DITS_H__

#include <vector>
#include <string>
#include <sstream>

#include "quspin/basis/types.h"
#include "quspin/basis/bitbasis/info.h"

namespace quspin::basis {

namespace constants {

static const quspin::basis::dit_integer_t bits[255] = {
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

static const quspin::basis::dit_integer_t mask[255] = {
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
struct dit_set { // thin wrapper used for convience
    
    int lhss;
    I mask;
    dit_integer_t bits;

    typedef I bitset_t;
    I content;

    dit_set(I _content,const int _lhss, const I _mask, const dit_integer_t bits) : 
    content(_content), 
    lhss(_lhss), 
    mask(_mask), 
    bits(bits) {}

    dit_set(I _content,const int _lhss) : 
    content(_content), 
    lhss(_lhss), 
    mask(constants::mask[_lhss]), 
    bits(constants::bits[_lhss]) {}

    dit_set(const dit_set<I>& other) : 
    content(other.content), 
    lhss(other.lhss), 
    mask(other.mask), 
    bits(other.bits) {}

    dit_set(const std::vector<dit_integer_t>& dits,const int _lhss):
    lhss(_lhss),
    mask(constants::mask[_lhss]),
    bits(constants::bits[_lhss])
    {
        content = 0;
        for(int i=0;i<dits.size();i++){
            content |= (I(dits[i]) << i*bits);
        }
    }

    std::vector<dit_integer_t> to_vector(const int length=0){
        const int niter = (length>0 ? length : bit_info<I>::bits/bits);

        std::vector<dit_integer_t> out(niter);
        for(int i=0;i<niter;++i){
            int shift = i*bits;
            out[i] = integer<dit_integer_t,I>::cast((content >> shift) & mask);
        }
        return out;
    }

    std::string to_string(const int length=0){
        auto dit_vec = to_vector(length);
        std::stringstream out;
        for(auto ele : dit_vec){
            out << (int)ele; 
        }
        return out.str();
    }

};

template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int i){
    return integer<int,I>::cast( (s.content >> (i * s.bits)) & s.mask);
}

template<typename I>
int get_sub_bitstring(const dit_set<I>& s,const int * locs,const int nlocs){
    int out = 0;
    for(int i=nlocs-1;i>0;--i){
        out += get_sub_bitstring(s,locs[i]) % s.lhss;
        out *= s.lhss;
        /* // implementation with padding for congituous blocks 
        out |= bit_basis<int,I>(get_sub_bitstring(s,locs[i]));
        out <<= s.bits;
        */
    }
    out += get_sub_bitstring(s,locs[0]);

    return out;
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,const int in,const int i){
    const int shift =  i * s.bits;
    const I r = s.content ^ ( ( I(in) << shift ) ^ s.content) & (s.mask << shift);
    return  dit_set<I>(r,s.lhss,s.mask,s.bits);
}

template<typename I>
dit_set<I> set_sub_bitstring(const dit_set<I>& s,int in,const int * locs,const int nlocs){
    I out = s.content;
    I in_I = I(in);
    for(int i=0;i<nlocs;++i){
        const int shift = s.bits * locs[i];
        out ^= (((in_I % s.lhss) << shift ) ^ s.content)  &  ( s.mask << shift );
        in_I /= s.lhss;
    }

    return dit_set<I>(out,s.lhss,s.mask,s.bits);
}

template<typename I>
inline bool operator<(const dit_set<I>& lhs, const dit_set<I>& rhs){return lhs.content < rhs.content;}

template<typename I>
inline bool operator>(const dit_set<I>& lhs, const dit_set<I>& rhs){return lhs.content > rhs.content;}

template<typename I>
inline bool operator==(const dit_set<I>& lhs, const dit_set<I>& rhs){return lhs.content == rhs.content;}


} // end namespace quspin::basis


#ifdef QUSPIN_UNIT_TESTS

namespace quspin::basis { // explicit instantiation for code coverage

template struct dit_set<uint8_t>; 
template int get_sub_bitstring<uint8_t>(const dit_set<uint8_t>&,const int);
template int get_sub_bitstring<uint8_t>(const dit_set<uint8_t>&,const int*, const int);
template dit_set<uint8_t> set_sub_bitstring<uint8_t>(const dit_set<uint8_t>&,const int, const int);
template dit_set<uint8_t> set_sub_bitstring<uint8_t>(const dit_set<uint8_t>&,const int,const int *,const int);
template bool operator< <uint8_t>(const dit_set<uint8_t>&, const dit_set<uint8_t>&);
template bool operator> <uint8_t>(const dit_set<uint8_t>&, const dit_set<uint8_t>&);
template bool operator== <uint8_t>(const dit_set<uint8_t>&, const dit_set<uint8_t>&);

}

TEST_CASE("get_bit_substring") {
    using namespace quspin::basis;

    dit_set<uint8_t> state(0b01100100,3); // possible states: 00 01 10

    CHECK(get_sub_bitstring(state,0) == 0);
    CHECK(get_sub_bitstring(state,1) == 1);
    CHECK(get_sub_bitstring(state,2) == 2);
    CHECK(get_sub_bitstring(state,3) == 1);

    int l1[2] = {0,1};
    int l2[3] = {0,2,1};

    CHECK(get_sub_bitstring(state,l1,2) == 0 + 3 * 1);
    CHECK(get_sub_bitstring(state,l2,3) == 0 + 3 * 2 + 9 * 1);

}

TEST_CASE("set_sub_bitstring") {
    using namespace quspin::basis;

    dit_set<uint8_t> state(0b01100100,3); // possible states: 00 01 10

    dit_set<uint8_t> result = set_sub_bitstring(state,1,0);
    CHECK(result.content == 0b01100101);

    result = set_sub_bitstring(state,1,2);
    CHECK(result.content == 0b01010100);

    int l1[2] = {0,1};
    int in1 = 1 + (3*2);
    result = set_sub_bitstring(state,in1,l1,2);
    CHECK(result.content == 0b01101001);

    int l2[3] = {0,3,1};
    int in2 = 2 + (3*2) + (9*0);
    result = set_sub_bitstring(state,in2,l2,3);
    CHECK(result.content == 0b10100010);
}

TEST_CASE("operators") {
    using namespace quspin::basis;

    dit_set<uint8_t> s1(0b1010111,3);
    dit_set<uint8_t> s2(0b1010111,3);
    dit_set<uint8_t> s3(0b1011111,3);
    dit_set<uint8_t> s4(0b1010011,3);


    CHECK(s1==s2);
    CHECK(s3>s1);
    CHECK(s4<s1);
}

TEST_CASE("from_/to_vector") {
    using namespace quspin::basis;

    dit_set<uint8_t> s(0b01100100,3);
    std::vector<dit_integer_t> dits = {0,1,2,1};
    
    CHECK(dit_set<uint8_t>(dits,3) == s);
    CHECK(s.to_vector() == dits);
}

TEST_CASE("to_string") {
    using namespace quspin::basis;

    dit_set<uint8_t> s(0b01100100,3);
    std::string dits = "0121"; // note reverse order
    CHECK(s.to_string() == dits);
}

#endif


#endif