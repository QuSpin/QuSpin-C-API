
# Make test executable
add_executable(test_cpp test.cpp)

sif(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(test_cpp PRIVATE --coverage -fno-inline -fno-inline-small-functions -fno-default-inline)
    target_link_libraries(test_cpp gcov)
endif()

add_test(NAME test_cpp
        COMMAND test_cpp)

