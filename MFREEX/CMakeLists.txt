file(GLOB_RECURSE SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)


SET(EIGEN_SRC_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib/Eigen/src)
SET(EIGEN_HEADER_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../lib/Eigen/include)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/lib)
include_directories(${EIGEN_HEADER_DIRECTORY})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${CVLIB_SOURCE_DIR}/test/include)
include_directories(${GPC_SOURCE_DIR})
include_directories(${util_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Threads)

SET(GCC_COVERAGE_COMPILE_FLAGS " -O3 -fPIC -std=c99 ")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})


add_library(MFREEX STATIC ${SRCS} 


	${EIGEN_SRC_DIRECTORY}/dsytrd3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevc3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevj3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevq3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevv3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevh3.c 
	${EIGEN_SRC_DIRECTORY}/dsyevd3.c 

	${HEADERS}

)

