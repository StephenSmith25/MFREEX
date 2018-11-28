cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(CVLIB VERSION 0.1.0 LANGUAGES C)


SET(CVLIB_LIB_OUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(CVLIB_BIN_OUT_PATH ${CVLIB_SOURCE_DIR}/test/bin)



add_subdirectory(lib)
