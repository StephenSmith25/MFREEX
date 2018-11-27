cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(MFREEX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Whether to compile test program
OPTION(PROBLEM_TEST "Compile test program" ON)




add_subdirectory(lib)
add_subdirectory(MFREEX)


if(PROBLEM_TEST)
	add_subdirectory(problems/test)
endif(PROBLEM_TEST)	
