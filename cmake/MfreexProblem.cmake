# mfreex_add_problem(<name> SOURCES <src>... [OUTPUT_DIR <dir>])
#
# Builds an example problem executable into ${CMAKE_BINARY_DIR}/bin/<dir>
# (default <dir> = <name>) and copies any input data files found next to the
# sources (*.nodes, *.segs, *.poly, *.coords, *.boundary) alongside it, since
# the problems open their inputs relative to the working directory.
function(mfreex_add_problem name)
  cmake_parse_arguments(ARG "" "OUTPUT_DIR" "SOURCES" ${ARGN})
  if(NOT ARG_SOURCES)
    message(FATAL_ERROR "mfreex_add_problem(${name}): SOURCES is required")
  endif()
  if(NOT ARG_OUTPUT_DIR)
    set(ARG_OUTPUT_DIR ${name})
  endif()

  set(out_dir "${CMAKE_BINARY_DIR}/bin/${ARG_OUTPUT_DIR}")

  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE MFREEX::MFREEX)
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${out_dir}")

  file(GLOB data_files CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.nodes"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.segs"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.poly"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.coords"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.boundary")
  if(data_files)
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory "${out_dir}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${data_files} "${out_dir}"
      COMMENT "Copying input files for ${name}"
      VERBATIM)
  endif()
endfunction()
