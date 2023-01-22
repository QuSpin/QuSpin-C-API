
![C++ CI](https://github.com/QuSpin/QuSpin-Core/actions/workflows/cpp_api_ci.yml/badge.svg)
![C++ CI](https://github.com/QuSpin/QuSpin-Core/actions/workflows/cython_api_ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/QuSpin/QuSpin-Core/branch/main/graph/badge.svg?token=RMRQBFPFT6)](https://codecov.io/gh/QuSpin/QuSpin-Core)

# QuSpin-Core (work in progress)
Low-level C++ and Cython API for QuSpin. Note currently the Cython CI just tests if the API builds. 

## TODO:

~~1. generate code coverage results for both Cython and C++ tests.~~

2. create unit tests for C++ API. 

    a. This is done inside the header files. Here is an example from [`src/quspin_core/includes/quspin/basis/bitbasis/bits.h`](https://github.com/QuSpin/QuSpin-Core/blob/83e273776a6421ca58b5a20302e8a1bdd5950163/src/quspin_core/includes/quspin/basis/bitbasis/bits.h#L94)
    
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
	All one has to do is define the macro boundary `#ifdef QUSPIN_UNIT_TESTS ... #endif` and write tests inside the body of the if statement. Following the [doctest](https://github.com/doctest/doctest) instructions to set up tests. To run the tests you need CMAKE and you can run the following commands in the terminal from QuSpin-Core root directory:
	
	i. `cmake -B build .`
	
	ii. `cmake --build build/ `
	
	iii `cd build/`
	
	iv. `ctest `
	
	it will output something like:
	
	```
	Test project .../QuSpin-Core/build
	Start 1: test_cpp
	1/1 Test #1: test_cpp .........................   Passed    0.07 sec

	100% tests passed, 0 tests failed out of 1

	Total Test time (real) =   0.11 sec
	```
	
	These tests will automatically get run every time the code is updated on github so make these test as light as possible. The idea is to test the most basic functionality with simple cases that can be worked out by hand. 

3. design C++ ABI (Application backend interface) 

    ~~a. preliminary design concept and implementation.~~
    
    b. revise during development of Cython API

4. Cython API:

    a. design
    
    b. test
    
5. fermion core ABI/API:
    
    a. design
    
    b. test

6. Start to plan integration for QuSpin-Core into QuSpin

    a. move relevant C++ code from QuSpin to QuSpin-Core

		 i. redesign

	    ii. test
        
    b. replace old backend withthe  new backend
    
	      i. removing old cython code for 1D basis sets
	        
	      ii.redesign thee  backend interface 
	        
	      iii. potentially remove old arguments from functions


