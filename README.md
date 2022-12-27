![C++ CI](https://github.com/QuSpin/QuSpin-Core/actions/workflows/cpp_api_ci.yml/badge.svg)
![C++ CI](https://github.com/QuSpin/QuSpin-Core/actions/workflows/cython_api_ci.yml/badge.svg)

# QuSpin-Core (work in progress)
Low-level C++ API for QuSpin 

## TODO:

1. generate code coverage results for both Cython and C++ tests. 
    a. Only runs on linux github action runners
    c. for C++ API use gcov (comes standard with gcc) integrate into CMAKE
        i. https://jhbell.com/using-cmake-and-gcov
    d. for Cython API use coverage.py, see (NOTE: This is a bit lower on priorities since cython code hasn't even been designed.)
        i. https://cython.readthedocs.io/en/latest/src/tutorial/profiling_tutorial.html#enabling-coverage-analysis
        ii. https://coverage.readthedocs.io/en/7.0.1/
        iii. needs to be integrated into install/build commands in setup.py

2. create unit tests for C++ API. 
    a. This is done inside the header files, e.g. here is an example from `src/quspin_core/includes/quspin/basis/bitbasis/bits.h`
        ```
        #ifdef QUSPIN_UNIT_TESTS

        TEST_CASE("get_bit_substring") {

            using namespace quspin::basis;

            bit_set<uint8_t> state(0b1010111);

            CHECK(get_sub_bitstring(state,0) == 1);
            CHECK(get_sub_bitstring(state,3) == 0);
            CHECK(get_sub_bitstring(state,7) == 0);
            CHECK(get_sub_bitstring(state,9) == 0);

            int l1[2] = {0,1};
            int l2[3] = {0,3,6};

            CHECK(get_sub_bitstring(state,l1,2) == 0b11);
            CHECK(get_sub_bitstring(state,l2,3) == 0b101);

        }

        #endif
        ```

3. design C++ ABI (Application backend interface) 
4. Cython API:
    a. design
    b. test
5. fermion core API:
    a. design
    b. test
6. Start to plan integration for QuSpin-Core into QuSpin
    a. move relavant C++ code from QuSpin to QuSpin-Core
        i. redesign
        ii. test
    b. replace old backend with new backend
        i. removing old cython code for 1D basis sets
        ii. redesign backend interface 
        iii. potentially remove old arguments from functions


