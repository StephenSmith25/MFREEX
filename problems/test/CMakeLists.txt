cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test)


include_directories(${Meschach_SOURCE_DIR})
link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})



set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/test)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads)

SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -std=c99")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})


add_executable(test test.c)
target_link_libraries(test meschach m)


