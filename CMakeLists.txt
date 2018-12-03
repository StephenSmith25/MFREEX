cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(MFREEX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Whether to compile test program
OPTION(PROBLEM_BEAM "Compile Beam program" ON)
OPTION(PROBLEM_BEAM_UL "Compile updated lagrangian version of Beam program" ON)




add_subdirectory(lib)
add_subdirectory(MFREEX)


if(PROBLEM_BEAM)
	add_subdirectory(problems/beam)
endif(PROBLEM_BEAM)	

if(PROBLEM_BEAM_UL)
	add_subdirectory(problems/beamUL)
endif(PROBLEM_BEAM_UL)	