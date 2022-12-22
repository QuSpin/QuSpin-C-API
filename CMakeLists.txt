cmake_minimum_required(VERSION 3.21)
project(quspin_c_api VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(doctest 2.4.9 REQUIRED)
add_subdirectory(test/test_cpp)
