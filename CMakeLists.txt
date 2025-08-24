cmake_minimum_required(VERSION 3.22)

project("Critical Cells")

set (CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
# for later use
find_package(Eigen3 REQUIRED)
# find_package(MPI)
find_package(CGAL REQUIRED)

# Option for legacy ICC
option(BUILD_LEGACY_ICC "Build with legacy ICC code" OFF)

#Subdirectory
add_subdirectory(deps)

# Create library for CriticalCells
add_library(CriticalCells SHARED criticalCells.cpp)

target_include_directories(CriticalCells PUBLIC ${CMAKE_SOURCE_DIR}/deps)

# target_link_libraries(CriticalCells PUBLIC Eigen3::Eigen hopcroft_karp readInput distMat OpenMP::OpenMP_CXX biGraphMorseMatch CGAL::CGAL)
target_link_libraries(CriticalCells PUBLIC
    readInput
    DistanceMatrix
    SimplexEnumerator
    ApproximateMorseMatching
    MaximumMorseMatching
    QuotientAndExpand
    # External libraries
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    CGAL::CGAL
)

# Link legacy ICC library if option is enabled
if (BUILD_LEGACY_ICC)
    target_compile_definitions(CriticalCells PUBLIC BUILD_ICC_LEGACY)
    target_link_libraries(CriticalCells PUBLIC hopcroft_karp)
endif()


# if(MPI_FOUND)
#     target_link_libraries(CriticalCells PUBLIC MPI::MPI_CXX)
#     add_compile_definitions(MPI_ENABLED=1)
# endif()

# Create executable
add_executable(CritCells main.cpp)
target_link_libraries(CritCells PRIVATE CriticalCells)
