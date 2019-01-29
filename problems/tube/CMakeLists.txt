cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(tube)


include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${util_SOURCE_DIR}/include)

include_directories(${GPC_SOURCE_DIR})
include_directories(${MFREEX_SOURCE_DIR}/MFREEX/include)
include_directories(${triangulate_SOURCE_DIR}/include)


link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})



set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/tube)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads)


add_executable(tube tube.c)
target_link_libraries(tube CVLIB MFREEX EFSFLIB triangulate meschach m gpc util)


# copy input files to location
# add_custom_command(
#         TARGET tube POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy
#                 ${CMAKE_CURRENT_SOURCE_DIR}/tube.segs
#                 ${EXECUTABLE_OUTPUT_PATH}/tube.segs)
# add_custom_command(
#         TARGET tube POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy
#                 ${CMAKE_CURRENT_SOURCE_DIR}/tube.nodes
#                 ${EXECUTABLE_OUTPUT_PATH}/tube.nodes)
