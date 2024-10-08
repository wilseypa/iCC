cmake_minimum_required(VERSION 3.22)

project("Critical Cells")

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

#Add Subdirectory
add_subdirectory(deps)

# Create library for CriticalCells
add_library(CriticalCells SHARED criticalCells.cpp criticalCells.hpp)
target_link_libraries(CriticalCells PUBLIC Eigen3::Eigen hopcroft_karp readInput distMat OpenMP::OpenMP_CXX)

# Create executable
add_executable(CritCells main.cpp)
target_link_libraries(CritCells PRIVATE CriticalCells)
