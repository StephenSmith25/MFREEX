cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(MFREEX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Whether to compile test program
OPTION(PROBLEM_BEAM "Compile Beam program" ON)
OPTION(PROBLEM_BEAM_UL "Compile updated lagrangian version of Beam program" ON)
OPTION(PROBLEM_PREFORM "Compile preform program" ON)
OPTION(PROBLEM_BIAXIAL "Compile biaxial program" ON)
OPTION(PROBLEM_TUBE "Compile tube program" OFF)
OPTION(PROBLEM_CYLINDER "Compile cylinder program" ON)
OPTION(PROBLEM_TAYLOR_BAR "Compile cylinder program" ON)
OPTION(PROBLEM_TEST "Compile TEST program" OFF)
OPTION(PROBLEM_PREFORM_UPDATED "Compile UPDATED PREFORM program" ON)
OPTION(PROBLEM_RUBBER "Compile RUBBER program" ON)
OPTION(PROBLEM_RUBBER_COMP "Compile RUBBER_COMP program" ON)
OPTION(PROBLEM_PUNCH "Compile PUNCH program" ON)
OPTION(PROBLEM_CELL_SEARCH "Compile CELL_SEARCH program" ON)


add_subdirectory(lib)
add_subdirectory(MFREEX)


if(PROBLEM_BEAM)
	add_subdirectory(problems/beam)
endif(PROBLEM_BEAM)	

if(PROBLEM_BEAM_UL)
	add_subdirectory(problems/beam_updated)
endif(PROBLEM_BEAM_UL)	

if(PROBLEM_PREFORM)
	add_subdirectory(problems/preform)
endif(PROBLEM_PREFORM)	

if(PROBLEM_PREFORM_UPDATED)
	add_subdirectory(problems/preform_updated)
endif(PROBLEM_PREFORM_UPDATED)	
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

if(PROBLEM_RUBBER)
	add_subdirectory(problems/rubber)
endif(PROBLEM_RUBBER)	


if(PROBLEM_RUBBER_COMP)
	add_subdirectory(problems/rubber_comp)
endif(PROBLEM_RUBBER_COMP)	


if(PROBLEM_TEST)
	add_subdirectory(problems/test)
endif(PROBLEM_TEST)	

if(PROBLEM_PUNCH)
	add_subdirectory(problems/rigid_punch)
endif(PROBLEM_PUNCH)	


if(PROBLEM_CELL_SEARCH)
	add_subdirectory(problems/cell_search)
endif(PROBLEM_CELL_SEARCH)	
