add_library(CVLIB STATIC
  src/generate_voronoi.c
  src/callVoronoi.c
  src/clipVoronoi.c)
target_include_directories(CVLIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(CVLIB PUBLIC meschach gpc)
