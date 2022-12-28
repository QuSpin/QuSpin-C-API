#ifndef __QUSPIN_UTILS_FUNCTIONS_H__
#define __QUSPIN_UTILS_FUNCTIONS_H__

#include <vector>
#include <complex>

namespace quspin {


template<typename T> T real_value(const T v){return v;}

template<typename T> T real_value(const std::complex<T>& v){return v.real();}

template<typename T> T conj(const T v){return v;}

template<typename T> T conj(const std::complex<T>& v){return std::conj(v);}

}

#endif