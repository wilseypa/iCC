cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)

add_library(hopcroft_karp STATIC hopcroft_karp.cpp)
target_include_directories(hopcroft_karp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(biGraphMorseMatch STATIC biGraphMorseMatch.cpp)
target_include_directories(biGraphMorseMatch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(biGraphMorseMatch PUBLIC OpenMP::OpenMP_CXX)

add_library(readInput STATIC readInput.cpp)
target_include_directories(readInput PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(distMat STATIC distMat.cpp)