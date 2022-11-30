#ifndef __QUSPIN_BASIS_BITBASIS_TYPES_H__
#define __QUSPIN_BASIS_BITBASIS_TYPES_H__


#include "boost/multiprecision/cpp_int.hpp"
#include "boost/numeric/conversion/cast.hpp"


namespace BitBasis {

typedef boost::multiprecision::uint32_t uint32_t;
typedef boost::multiprecision::uint64_t uint64_t;
typedef boost::multiprecision::uint128_t uint128_t;
typedef boost::multiprecision::uint256_t uint256_t;
typedef boost::multiprecision::uint512_t uint512_t;
typedef boost::multiprecision::uint1024_t uint1024_t;

typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<2048, 2048, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void> > uint2048_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<4096, 4096, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void> > uint4096_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<8192, 8192, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void> > uint8192_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<16384, 16384, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void> > uint16384_t;

}

#endif