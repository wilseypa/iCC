project("Parallel_Match")
cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(bi_graph STATIC bi_graph.cpp)
add_library(cycle_rm STATIC cycle_rm.cpp)

# Create executable
add_executable(Parallel_Match main.cpp)
target_link_libraries(Parallel_Match PRIVATE bi_graph cycle_rm OpenMP::OpenMP_CXX)
