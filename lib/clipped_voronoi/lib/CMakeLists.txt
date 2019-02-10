cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(lib)


include_directories(${Meschach_SOURCE_DIR})
include_directories(${GPC_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)

SET(GCC_COVERAGE_COMPILE_FLAGS " -O3 -std=c99 ")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/library)

add_library(CVLIB STATIC src/generate_voronoi.c src/callVoronoi.c src/clipVoronoi.c)

