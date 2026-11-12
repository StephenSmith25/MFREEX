add_library(EFSFLIB STATIC
  src/polynomial_basis.c
  src/mls_shapefunction.c
  src/neighbours.c
  src/setDomain.c
  src/weight_function.c
  src/v_outer_product.c
  src/quartic_spline.c
  src/cubic_spline.c)
target_include_directories(EFSFLIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(EFSFLIB PUBLIC meschach)
