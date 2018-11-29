project(triangulate)


file(GLOB_RECURSE SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${Meschach_SOURCE_DIR})


add_library(triangulate STATIC ${SRCS} ${HEADERS})