
# Make test executable
add_executable(test_cpp test.cpp)

add_test(NAME test_cpp
         COMMAND test_cpp)