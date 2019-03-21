cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)


project(libs)


SET(MFREEX_LIBS_OUT_PATH ${MFREEX_SOURCE_DIR}/lib)

add_subdirectory(gpc)
add_subdirectory(Meschach)
add_subdirectory(clipped_voronoi)
add_subdirectory(efsflib)
add_subdirectory(util)
add_subdirectory(triangulate)
add_subdirectory(Eigen)
add_subdirectory(thread_pool)
