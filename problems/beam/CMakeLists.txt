cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(beam)


include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${util_SOURCE_DIR}/include)

include_directories(${GPC_SOURCE_DIR})
include_directories(${MFREEX_SOURCE_DIR}/MFREEX/include)
include_directories(${triangulate_SOURCE_DIR}/include)


link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


configure_file(square.segs square.segs COPYONLY)
configure_file(square.nodes square.nodes COPYONLY)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/beam)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads)


add_executable(beam beam.c)
target_link_libraries(beam CVLIB MFREEX EFSFLIB triangulate meschach m gpc util)


# copy input files to location
add_custom_command(
        TARGET beam POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/square.segs
                ${EXECUTABLE_OUTPUT_PATH}/square.segs)
add_custom_command(
        TARGET beam POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/square.nodes
                ${EXECUTABLE_OUTPUT_PATH}/square.nodes)