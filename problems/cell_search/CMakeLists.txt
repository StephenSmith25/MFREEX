cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(cell_search)


include_directories(${triangulate_SOURCE_DIR}/include)
include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${util_SOURCE_DIR}/include)
include_directories(${EIGEN_SOURCE_DIR}/include)

include_directories(${GPC_SOURCE_DIR})
include_directories(${MFREEX_SOURCE_DIR}/MFREEX/include)
include_directories(${thread_pool_SOURCE_DIR})


link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})



set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/cell_search)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads REQUIRED)


SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -g -std=c99 -pthread -lpthread")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

add_executable(cell_search cell_search.c)
target_link_libraries(cell_search  CVLIB thread_pool triangulate MFREEX Eigen EFSFLIB meschach m gpc util )


