project(util)


file(GLOB_RECURSE SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${Meschach_SOURCE_DIR})

SET(GCC_COVERAGE_COMPILE_FLAGS " -O3 -std=c99 ")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

add_library(util STATIC ${SRCS} ${HEADERS})