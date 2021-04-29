cmake_minimum_required(VERSION 2.8)
#set(CMAKE_CXX_COMPILER "icpc")
project(cpu_vah)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require GSL for this project:
find_package(GSL REQUIRED)

# gcc flags
#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -march=native -std=c++11")

# intel flags
set(CMAKE_CXX_FLAGS "-O3 -qopenmp -std=c++11")

include_directories(${GSL_INCLUDE_DIR})

# header files
include_directories(rhic/include)

# cornelius external package
include_directories(cornelius-c++-1.3)

add_subdirectory(rhic/src)