# Runs an example problem in a scratch directory and compares its output
# files against stored reference output with numdiff.
#
# Variables (passed with -D):
#   EXE       executable to run
#   WORKDIR   scratch directory (recreated)
#   INPUTS    ;-separated input files to copy into WORKDIR
#   REF_DIR   directory with reference output files
#   FILES     ;-separated output files to compare
#   NUMDIFF   numdiff executable
#   REL_TOL / ABS_TOL  tolerances

foreach(v EXE WORKDIR REF_DIR FILES NUMDIFF REL_TOL ABS_TOL)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "RunGolden: ${v} not set")
  endif()
endforeach()

file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
foreach(f IN LISTS INPUTS)
  file(COPY "${f}" DESTINATION "${WORKDIR}")
endforeach()
# Output folders some problems expect to exist
foreach(d Displacement Stress Strain Domain Mesh srcVolume History)
  file(MAKE_DIRECTORY "${WORKDIR}/${d}")
endforeach()

execute_process(COMMAND "${EXE}"
  WORKING_DIRECTORY "${WORKDIR}"
  RESULT_VARIABLE rc
  OUTPUT_FILE "${WORKDIR}/stdout.log"
  ERROR_FILE "${WORKDIR}/stderr.log")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${EXE} exited with '${rc}' (logs in ${WORKDIR})")
endif()

set(failed "")
foreach(f IN LISTS FILES)
  execute_process(COMMAND "${NUMDIFF}" "${REF_DIR}/${f}" "${WORKDIR}/${f}" ${REL_TOL} ${ABS_TOL}
    RESULT_VARIABLE drc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  message("${f}: ${out}${err}")
  if(NOT drc EQUAL 0)
    list(APPEND failed ${f})
  endif()
endforeach()
if(failed)
  message(FATAL_ERROR "Output differs from reference: ${failed}")
endif()
