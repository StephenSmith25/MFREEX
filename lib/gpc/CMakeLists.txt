cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GPC)

SET(GPC_LIB_OUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


add_library(gpc STATIC gpc.c)