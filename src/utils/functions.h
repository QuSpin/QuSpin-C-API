#ifndef __QUSPIN_UTILS_FUNCTIONS_H__
#define __QUSPIN_UTILS_FUNCTIONS_H__

namespace quspin {

template<typename dits_or_bits>
inline
bool operator>(dits_or_bits& lhs,dits_or_bits& rhs){return lhs.content < rhs.content;}

template<typename dits_or_bits>
inline
bool operator<(dits_or_bits& lhs,dits_or_bits& rhs){return lhs.content < rhs.content;}

template<typename dits_or_bits>
inline
bool operator==(dits_or_bits& lhs,dits_or_bits& rhs){return lhs.content < rhs.content;}

template<typename T>
T real_value(const T v){return v;}

template<typename T>
T real_value(const std::complex<T>& v){return v.real();}

}



#endif