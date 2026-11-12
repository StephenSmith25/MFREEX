add_library(triangulate STATIC
  src/triangle.c
  src/trigen.c
  src/connect_segments.c)
target_include_directories(triangulate
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(triangulate PUBLIC meschach)
