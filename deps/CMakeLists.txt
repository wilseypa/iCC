cmake_minimum_required(VERSION 3.22)

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)

if (BUILD_LEGACY_ICC)
    add_library(hopcroft_karp STATIC hopcroft_karp.cpp)
    target_include_directories(hopcroft_karp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_library(readInput STATIC readInput.cpp)
target_include_directories(readInput PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Library for the DistanceMatrix class and its constructors
add_library(DistanceMatrix STATIC DistanceMatrix.cpp)
target_include_directories(DistanceMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Library for the Approximate Morse Matching strategy
add_library(ApproximateMorseMatching STATIC ApproximateMorseMatching.cpp)
target_include_directories(ApproximateMorseMatching PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Library for the Maximum Morse Matching strategy
add_library(MaximumMorseMatching STATIC MaximumMorseMatching.cpp)
target_include_directories(MaximumMorseMatching PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Library for the Quotient and Expand workflow
add_library(QuotientAndExpand STATIC QuotientAndExpand.cpp)
target_include_directories(QuotientAndExpand PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})