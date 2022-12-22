#ifndef __QUSPIN_UTILS_FUNCTIONS_H__
#define __QUSPIN_UTILS_FUNCTIONS_H__

#include <vector>

namespace quspin {


template<typename T>
T real_value(const T v){return v;}

template<typename T>
T real_value(const std::complex<T>& v){return v.real();}

}


template<typename I>
std::hash(quspin::basis::bitbasis::bit_set<I> &state){
    return std::hash(state.content);
}


template<typename I>
std::hash(quspin::basis::bitbasis::dit_set<I> &state){
    return std::hash(state.content);
}

#endif