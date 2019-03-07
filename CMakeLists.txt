cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(MFREEX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Whether to compile test program
OPTION(PROBLEM_BEAM "Compile Beam program" OFF)
OPTION(PROBLEM_BEAM_UL "Compile updated lagrangian version of Beam program" OFF)
OPTION(PROBLEM_PREFORM "Compile preform program" ON)
OPTION(PROBLEM_BIAXIAL "Compile biaxial program" ON)
OPTION(PROBLEM_TUBE "Compile tube program" OFF)
OPTION(PROBLEM_CYLINDER "Compile cylinder program" ON)
OPTION(PROBLEM_TAYLOR_BAR "Compile cylinder program" ON)
OPTION(PROBLEM_TEST "Compile TEST program" OFF)


add_subdirectory(lib)
add_subdirectory(MFREEX)


if(PROBLEM_BEAM)
	add_subdirectory(problems/beam)
endif(PROBLEM_BEAM)	

if(PROBLEM_BEAM_UL)
	add_subdirectory(problems/beamUL)
endif(PROBLEM_BEAM_UL)	

if(PROBLEM_PREFORM)
	add_subdirectory(problems/preform)
endif(PROBLEM_PREFORM)	

if(PROBLEM_BIAXIAL)
	add_subdirectory(problems/Biaxial)
endif(PROBLEM_BIAXIAL)	

if(PROBLEM_TUBE)
	add_subdirectory(problems/tube)
endif(PROBLEM_TUBE)	


if(PROBLEM_CYLINDER)
	add_subdirectory(problems/Cylinder)
endif(PROBLEM_CYLINDER)	

if(PROBLEM_TAYLOR_BAR)
	add_subdirectory(problems/taylor_bar)
endif(PROBLEM_TAYLOR_BAR)	


if(PROBLEM_TEST)
	add_subdirectory(problems/test)
endif(PROBLEM_TEST)	