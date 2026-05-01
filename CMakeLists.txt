cmake_minimum_required(VERSION 3.22)

project("Critical Cells")

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)

include(FetchContent)
find_package(CLI11 CONFIG QUIET)
if (NOT CLI11_FOUND)
    FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.6.2
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(CLI11)
endif()


# Option for alpha-complex support
option(BUILD_ALPHA_COMPLEX "Build alpha-complex support (requires CGAL and Eigen3)" OFF)
# Option for legacy ICC
option(BUILD_LEGACY_ICC "Build with legacy ICC code" OFF)

if (BUILD_ALPHA_COMPLEX)
    add_compile_definitions(BUILD_ALPHA_COMPLEX)
endif()

if (BUILD_LEGACY_ICC)
    add_compile_definitions(BUILD_LEGACY_ICC)
endif()

#Subdirectory
add_subdirectory(deps)

# Create library for CriticalCells
add_library(CriticalCells STATIC criticalCells.cpp)

target_include_directories(CriticalCells PUBLIC ${CMAKE_SOURCE_DIR}/deps)

target_link_libraries(CriticalCells PUBLIC
    readInput
    DistanceMatrix
    SimplexEnumerator
    ApproximateMorseMatching
    MaximumMorseMatching
    QuotientAndExpand
    OpenMP::OpenMP_CXX
)

# Create executable
add_executable(iCC main.cpp)
target_link_libraries(iCC PRIVATE CriticalCells CLI11::CLI11)

# Link legacy ICC library if option is enabled
if (BUILD_LEGACY_ICC)
    target_link_libraries(CriticalCells PUBLIC hopcroft_karp)
endif()


# if(MPI_FOUND)
#     target_link_libraries(CriticalCells PUBLIC MPI::MPI_CXX)
#     add_compile_definitions(MPI_ENABLED=1)
# endif()
