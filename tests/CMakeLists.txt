# MFREEX regression tests.  Run with:   ctest --test-dir build --output-on-failure
#
#   unit_*      property-based unit tests (patch tests, reproduction
#               conditions, energy consistency, analytic solutions)
#   golden_*    rerun example problems and compare with reference output
#   memcheck_*  unit tests under valgrind (if installed)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/unit)
include_directories(${MFREEX_SOURCE_DIR}/MFREEX/include)
include_directories(${EFSFLIB_SOURCE_DIR}/lib/include)
include_directories(${Meschach_SOURCE_DIR})
include_directories(${CVLIB_SOURCE_DIR}/lib/include)
include_directories(${GPC_SOURCE_DIR})
include_directories(${util_SOURCE_DIR}/include)
include_directories(${Eigen_SOURCE_DIR}/include)
include_directories(${triangulate_SOURCE_DIR}/include)
include_directories(${thread_pool_SOURCE_DIR})

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(Threads)

set(MFREEX_TEST_LIBS CVLIB MFREEX EFSFLIB triangulate util Eigen gpc thread_pool meschach m
    ${CMAKE_THREAD_LIBS_INIT})

set(UNIT_TESTS
  test_mls_shapefunctions
  test_nodal_integration
  test_support_domain
  test_weight_functions
  test_linear_algebra
  test_kinematics_materials
  test_hyperelastic_energy)

foreach(t IN LISTS UNIT_TESTS)
  add_executable(${t} unit/${t}.c)
  target_link_libraries(${t} ${MFREEX_TEST_LIBS})
  set_target_properties(${t} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
  string(REPLACE "test_" "unit_" name ${t})
  add_test(NAME ${name} COMMAND ${t})
  set_tests_properties(${name} PROPERTIES LABELS unit TIMEOUT 120 ENVIRONMENT OMP_NUM_THREADS=2)
endforeach()

# ---------------------------------------------------------------- golden --
add_executable(numdiff tools/numdiff.c)
target_link_libraries(numdiff m)
set_target_properties(numdiff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

set(GOLDEN_REL_TOL 1e-6)
set(GOLDEN_ABS_TOL 1e-9)

function(add_golden_test name exe_target inputs files)
  if(NOT TARGET ${exe_target})
    message(STATUS "golden_${name}: problem '${exe_target}' not built, test skipped")
    return()
  endif()
  add_test(NAME golden_${name}
    COMMAND ${CMAKE_COMMAND}
      -DEXE=$<TARGET_FILE:${exe_target}>
      -DWORKDIR=${CMAKE_BINARY_DIR}/tests/golden_${name}
      "-DINPUTS=${inputs}"
      -DREF_DIR=${CMAKE_CURRENT_SOURCE_DIR}/golden/${name}
      "-DFILES=${files}"
      -DNUMDIFF=$<TARGET_FILE:numdiff>
      -DREL_TOL=${GOLDEN_REL_TOL}
      -DABS_TOL=${GOLDEN_ABS_TOL}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/RunGolden.cmake)
  set_tests_properties(golden_${name} PROPERTIES LABELS golden TIMEOUT 300
    ENVIRONMENT OMP_NUM_THREADS=1)
endfunction()

add_golden_test(biaxial biaxial "" "Stress_11_nominal.txt")
add_golden_test(beam beam
  "${MFREEX_SOURCE_DIR}/problems/beam/square.nodes;${MFREEX_SOURCE_DIR}/problems/beam/square.segs"
  "nodes.csv;boundary.txt;search_cells.csv")

# -------------------------------------------------------------- memcheck --
find_program(VALGRIND valgrind)
if(VALGRIND)
  foreach(t test_support_domain test_linear_algebra test_hyperelastic_energy)
    string(REPLACE "test_" "memcheck_" name ${t})
    add_test(NAME ${name}
      COMMAND ${VALGRIND} -q --error-exitcode=1 --track-origins=yes $<TARGET_FILE:${t}>)
    set_tests_properties(${name} PROPERTIES LABELS memcheck TIMEOUT 600
      ENVIRONMENT OMP_NUM_THREADS=1)
  endforeach()
endif()
