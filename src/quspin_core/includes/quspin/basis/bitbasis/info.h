#ifndef __QUSPIN_BASIS_BITBASIS_INFO_H__
#define __QUSPIN_BASIS_BITBASIS_INFO_H__

#include "quspin/basis/types.h"

namespace quspin::basis {

template<typename I>
struct bit_info{};

#ifdef USE_BOOST
template<>
struct bit_info<uint16384_t>
{ enum {ld_bits=14,bits=16384,bytes=2048};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint8192_t>
{ enum {ld_bits=13,bits=4096,bytes=1024};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint4096_t>
{ enum {ld_bits=12,bits=4096,bytes=512};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint2048_t>
{ enum {ld_bits=11,bits=2048,bytes=256};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint1024_t>
{ enum {ld_bits=10,bits=1024,bytes=128};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint512_t>
{ enum {ld_bits=9,bits=512,bytes=64};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint256_t>
{ enum {ld_bits=8,bits=256,bytes=32};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint128_t>
{ enum {ld_bits=7,bits=128,bytes=16};
  typedef int bit_index_type;
};


#endif

template<>
struct bit_info<uint64_t>
{  enum {ld_bits=6,bits=64,bytes=8};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint32_t>
{  enum {ld_bits=5,bits=32,bytes=4};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint16_t>
{  enum {ld_bits=4,bits=16,bytes=2};
  typedef int bit_index_type;
};

template<>
struct bit_info<uint8_t>
{ enum {ld_bits=3,bits=8,bytes=1};
  typedef int bit_index_type;
};

#ifdef USE_BOOST

template<class J,class I>
struct integer {
  static J cast(const I s){
    try 
    {
      return boost::numeric_cast<J>(s);; // This conversion succeeds (is in range)
    }
    catch(boost::numeric::positive_overflow& e) {
      return -1;
    }
  }
};

template<class J>
struct integer<J,uint64_t> {
  static J cast(const uint64_t s){
    return static_cast<J>(s);
  }
};

template<class J>
struct integer<J,uint32_t> {
  static J cast(const uint32_t s){
    return static_cast<J>(s);
  }
};

template<class J>
struct integer<J,uint16_t> {
  static J cast(const uint16_t s){
    return static_cast<J>(s);
  }
};

template<class J>
struct integer<J,uint8_t> {
  static J cast(const uint8_t s){
    return static_cast<J>(s);
  }
};

#else

template<class J,class I>
struct integer {
  static J cast(const I s){
    return static_cast<J>(s);
  }
};

#endif

template<class T>
typename bit_info<T>::bit_index_type bit_pos(T x, typename bit_info<T>::bit_index_type *idx)
{
  typename bit_info<T>::bit_index_type n = 0;
  do {
    if (x & 1) *(idx++) = n;
    n++;
  } while (x >>= 1); 
  return n;
}

template<class T>
int inline bit_count(T v,const int l){
  v = v & (((~(T)0) >> 1) >> (bit_info<T>::bits - 1 - l));
  v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
  v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
  v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
  T res = (T)(v * ((T)~(T)0/255)) >> ((bit_info<T>::bytes - 1) * 8); // count
  return (int)res;

}

}

#ifdef QUSPIN_UNIT_TESTS

namespace quspin::basis {

template int bit_pos<uint8_t>(uint8_t,int *);
template int bit_count<uint8_t>(uint8_t,const int);

}

#endif

#endif