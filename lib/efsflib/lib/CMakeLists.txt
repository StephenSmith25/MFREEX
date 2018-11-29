cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(lib)


include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package (Threads)

set(LIBRARY_OUTPUT_PATH ${LIBS_OUT_PATH})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


add_library(EFSFLIB STATIC src/polynomial_basis.c src/mls_shapefunction.c src/neighbours.c src/setDomain.c src/weight_function.c
	src/v_outer_product.c src/quartic_spline.c src/cubic_spline.c)

#target_link_libraries(Shapefunction meschach m ${CMAKE_THREAD_LIBS_INIT} )