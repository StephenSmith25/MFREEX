# Joachim Kopp's analytical/iterative eigen-solvers for 3x3 symmetric matrices.
# (Directory is historically named "Eigen"; this is unrelated to the C++ Eigen library.)
add_library(eig3 STATIC
  src/dsyev2.c
  src/dsytrd3.c
  src/dsyevc3.c
  src/dsyevj3.c
  src/dsyevq3.c
  src/dsyevv3.c
  src/dsyevh3.c
  src/dsyevd3.c)
target_include_directories(eig3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(eig3 PUBLIC mfreex_compat)
