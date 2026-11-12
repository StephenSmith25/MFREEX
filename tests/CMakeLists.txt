# MFREEX regression tests.  Run with:   ctest --test-dir <build> --output-on-failure
#
#   unit_*      property-based unit tests (patch tests, reproduction
#               conditions, energy consistency, analytic solutions)
#   golden_*    rerun example problems and compare with reference output
#   memcheck_*  unit tests under valgrind (if installed)

set(MFREEX_TEST_OUT "${CMAKE_BINARY_DIR}/tests")

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
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/unit)
  target_link_libraries(${t} PRIVATE MFREEX::MFREEX)
  set_target_properties(${t} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${MFREEX_TEST_OUT}")
  string(REPLACE "test_" "unit_" name ${t})
  add_test(NAME ${name} COMMAND ${t})
  set_tests_properties(${name} PROPERTIES LABELS unit TIMEOUT 120 ENVIRONMENT OMP_NUM_THREADS=2)
endforeach()

# ---------------------------------------------------------------- golden --
add_executable(numdiff tools/numdiff.c)
if(NOT MSVC)
  target_link_libraries(numdiff PRIVATE m)
endif()
set_target_properties(numdiff PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${MFREEX_TEST_OUT}")

set(GOLDEN_REL_TOL 1e-6)
set(GOLDEN_ABS_TOL 1e-9)

# add_golden_test(<name> <problem target> <input files> <output files>)
function(add_golden_test name exe_target inputs files)
  if(NOT TARGET ${exe_target})
    message(STATUS "golden_${name}: problem '${exe_target}' not built, test skipped")
    return()
  endif()
  add_test(NAME golden_${name}
    COMMAND ${CMAKE_COMMAND}
      -DEXE=$<TARGET_FILE:${exe_target}>
      -DWORKDIR=${MFREEX_TEST_OUT}/golden_${name}
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

set(P ${PROJECT_SOURCE_DIR}/problems)
add_golden_test(biaxial biaxial "" "Stress_11_nominal.txt")
add_golden_test(beam beam "${P}/beam/square.nodes;${P}/beam/square.segs"
  "nodes.csv;boundary.txt;search_cells.csv")

# -------------------------------------------------------------- memcheck --
find_program(VALGRIND valgrind)
if(VALGRIND AND NOT MFREEX_SANITIZE)
  foreach(t test_support_domain test_linear_algebra test_hyperelastic_energy)
    string(REPLACE "test_" "memcheck_" name ${t})
    add_test(NAME ${name}
      COMMAND ${VALGRIND} -q --error-exitcode=1 --track-origins=yes $<TARGET_FILE:${t}>)
    set_tests_properties(${name} PROPERTIES LABELS memcheck TIMEOUT 600
      ENVIRONMENT OMP_NUM_THREADS=1)
  endforeach()
endif()
