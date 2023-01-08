#ifndef __QUSPIN_UTILS_FUNCTIONS_H__
#define __QUSPIN_UTILS_FUNCTIONS_H__

#include <vector>
#include <complex>

namespace quspin {

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


#ifdef USE_STD_COMPLEX

#include <complex>

namespace quspin {
    template<class T>
    std::complex<T> conj(const std::complex<T>&A){return std::conj(A);}

    template<class T>
    T conj(const T& A){return A;}

    template<class T>
    T real(const std::complex<T>&A){return A.real();}

    template<class T>
    T real(const T& A){return A;}

}
#endif

#ifdef QUSPIN_UNIT_TESTS

TEST_SUITE("quspin/utils/functions.h"){
    using namespace quspin;

    TEST_CASE("integer_pow"){
        CHECK(integer_pow<2,4>::value == 16);
        CHECK(integer_pow<3,3>::value == 27);
        CHECK(integer_pow<10,3>::value == 1000);
        CHECK(integer_pow<5,5>::value == 3125);
    }
}



#endif

#endif