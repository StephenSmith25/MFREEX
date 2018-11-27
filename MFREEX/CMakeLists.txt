file(GLOB_RECURSE SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${GPC_SOURCE_DIR})

add_library(MFREEX STATIC ${SRCS} ${HEADERS})