cmake_minimum_required(VERSION 3.24)
project(quspin_c_api VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

enable_testing()
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    doctest
    PREFIX ${CMAKE_BINARY_DIR}/doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)

# Expose required variable (DOCTEST_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(doctest source_dir)
set(DOCTEST_INCLUDE_DIR ${source_dir}/doctest CACHE INTERNAL "Path to include folder for doctest")
set(QUSPIN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/quspin_cpp_api)

# add it globally
include_directories(${DOCTEST_INCLUDE_DIR})
include_directories(${QUSPIN_INCLUDE_DIR})


add_subdirectory(test/test_cpp)

