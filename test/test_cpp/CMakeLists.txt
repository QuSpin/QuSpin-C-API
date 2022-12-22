
add_executable(test
    test.cpp
)

target_include_directories(test
    ${PROJECT_SOURCE_DIR}/src/quspin_cpp_api
    ${${doctest_INCLUDE_DIRS}}
)
