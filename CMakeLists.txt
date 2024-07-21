project("Critical Cells")
cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(hopcroft_karp STATIC deps/hopcroft_karp.cpp deps/hopcroft_karp.hpp)
add_library(readInput STATIC deps/readInput.cpp deps/readInput.hpp)
add_library(distMat STATIC deps/distMat.cpp)

# Create library for CriticalCells
add_library(CriticalCells SHARED criticalCells.cpp criticalCells.hpp)
target_link_libraries(CriticalCells PUBLIC Eigen3::Eigen hopcroft_karp readInput distMat OpenMP::OpenMP_CXX)

# Create executable
add_executable(CritCells main.cpp)
target_link_libraries(CritCells PRIVATE CriticalCells)
