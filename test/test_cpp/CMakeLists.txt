

# Prepare doctest for other targets to use
find_package(doctest REQUIRED)

# Make test executable
add_executable(test_cpp test.cpp)
target_link_libraries(test_cpp PRIVATE doctest::doctest)

add_test(NAME test_cpp
         COMMAND test_cpp)