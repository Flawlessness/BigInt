cmake_minimum_required(VERSION 3.21)
project(lab1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g -fprofile-arcs -ftest-coverage")

add_subdirectory(src)
add_subdirectory(tst)
add_subdirectory(lib/googletest)

