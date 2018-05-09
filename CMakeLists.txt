project(LP_MP-Cell-tracking)
cmake_minimum_required(VERSION 2.8.12)

# C++17
add_compile_options(-std=c++17)

# compiler options
add_definitions(-march=native)

include_directories(include)
add_subdirectory(external/LP_MP)
include_directories(external/PEGTL)
include_directories(external/MCF-SSP)

#add_subdirectory(test)
add_subdirectory(src)
