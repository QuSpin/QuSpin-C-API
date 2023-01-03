#ifndef __QUSPIN_UTILS_FUNCTIONS_H__
#define __QUSPIN_UTILS_FUNCTIONS_H__

#include <vector>
#include <complex>

namespace quspin {


template<typename T> T real(const T v){return v;}
template<typename T> T real(const std::complex<T>& v){return v.real();}
template<typename T> T conj(const T v){return v;}
template<typename T> T conj(const std::complex<T>& v){return std::conj(v);}

template<int base,std::size_t N>
struct integer_pow {
    enum {value = base * integer_pow<base,N-1>::value};
};

template<int base>
struct integer_pow<base,1>{
    enum {value = base};
};

template<int base>
struct integer_pow<base,0>{
    enum {value = 1};
};

}

#ifdef QUSPIN_UNIT_TESTS

TEST_SUITE("quspin/utils/functions.h"){
    TEST_CASE("integer_pow"){
        using namespace quspin;

        CHECK(integer_pow<2,4>::value == 16);
        CHECK(integer_pow<3,3>::value == 27);
        CHECK(integer_pow<10,3>::value == 1000);
        CHECK(integer_pow<5,5>::value == 3125);
        
    }
}



#endif

#endif