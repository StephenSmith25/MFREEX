cmake_minimum_required(VERSION 3.16)

project(MFREEX
  VERSION 0.2.0
  DESCRIPTION "Meshfree (element-free Galerkin) library for non-linear explicit dynamics"
  LANGUAGES C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)               # gnu99: M_PI, strdup, drand48 ...
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)    # for clangd / IDEs (replaces .clang_complete)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(MFREEX_USE_OPENMP     "Build with OpenMP threading (if available)"        ON)
option(MFREEX_ENABLE_WARNINGS "Compile MFREEX sources with -Wall -Wextra"         OFF)
option(MFREEX_WARNINGS_AS_ERRORS "Treat warnings in MFREEX sources as errors"    OFF)
option(MFREEX_SANITIZE       "Build with AddressSanitizer + UBSan (debugging)"   OFF)
option(MFREEX_BUILD_PROBLEMS "Build the example problems in problems/"           ON)
option(MFREEX_BUILD_TESTS    "Build the regression test suite (ctest)"           ON)

# Individual example problems (names kept from the original build)
option(PROBLEM_BEAM            "Build beam problem"                         ON)
option(PROBLEM_BEAM_UL         "Build updated-Lagrangian beam problem"      ON)
option(PROBLEM_BEAM_UL_LEGACY  "Build legacy beamUL problem"                OFF)
option(PROBLEM_PREFORM         "Build preform problem"                      ON)
option(PROBLEM_PREFORM_UPDATED "Build updated preform problem"              ON)
option(PROBLEM_BIAXIAL         "Build biaxial problem"                      ON)
option(PROBLEM_TUBE            "Build tube problem"                         OFF)
option(PROBLEM_CYLINDER        "Build cylinder problem"                     ON)
option(PROBLEM_TAYLOR_BAR      "Build Taylor bar problem"                   ON)
option(PROBLEM_TEST            "Build test problem"                         OFF)
option(PROBLEM_RUBBER          "Build rubber (tension) problem"             ON)
option(PROBLEM_RUBBER_COMP     "Build rubber (compression) problem"         ON)
option(PROBLEM_PUNCH           "Build rigid punch problem"                  ON)
option(PROBLEM_CELL_SEARCH     "Build cell search problem"                  ON)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(Threads REQUIRED)

if(MFREEX_USE_OPENMP)
  find_package(OpenMP COMPONENTS C)
  if(NOT OpenMP_C_FOUND)
    message(WARNING
      "OpenMP not found - building single-threaded.\n"
      "On macOS: `brew install libomp` and configure with "
      "-DOpenMP_ROOT=$(brew --prefix libomp) (or use the 'macos' preset).")
  endif()
endif()

# Interface target carrying project-wide usage requirements (OpenMP shim, libm)
add_library(mfreex_compat INTERFACE)
target_include_directories(mfreex_compat INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/compat/include")
if(MFREEX_USE_OPENMP AND OpenMP_C_FOUND)
  target_link_libraries(mfreex_compat INTERFACE OpenMP::OpenMP_C)
endif()
if(NOT MSVC)
  target_link_libraries(mfreex_compat INTERFACE m)
endif()
if(MFREEX_SANITIZE)
  target_compile_options(mfreex_compat INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(mfreex_compat INTERFACE -fsanitize=address,undefined)
endif()

include(MfreexProblem)

# ---------------------------------------------------------------------------
# Sub-projects
# ---------------------------------------------------------------------------
add_subdirectory(lib)
add_subdirectory(MFREEX)

if(MFREEX_BUILD_PROBLEMS)
  add_subdirectory(problems)
endif()

if(MFREEX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "MFREEX ${PROJECT_VERSION}")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C compiler : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
if(MFREEX_USE_OPENMP AND OpenMP_C_FOUND)
  message(STATUS "  OpenMP     : ON (${OpenMP_C_VERSION})")
else()
  message(STATUS "  OpenMP     : OFF (serial build)")
endif()
message(STATUS "")
