project("Critical Cells")
cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 23)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_library(readInput STATIC deps/readInput.cpp deps/readInput.hpp)
find_package(OpenMP)
include_directories(${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}) 

add_executable(CriticalCells criticalCells.cpp criticalCells.hpp)

target_include_directories(CriticalCells PRIVATE deps)
target_link_libraries(CriticalCells PRIVATE readInput OpenMP::OpenMP_CXX)