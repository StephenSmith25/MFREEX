cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test)


include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${util_SOURCE_DIR}/include)

include_directories(${GPC_SOURCE_DIR})
include_directories(${MFREEX_SOURCE_DIR}/MFREEX/include)
include_directories(${triangulate_SOURCE_DIR}/include)


link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/Test\ problem)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads)


add_executable(test_program test.c)
target_link_libraries(test_program CVLIB MFREEX EFSFLIB triangulate meschach m gpc util)