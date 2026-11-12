# Vendored third-party / support libraries.
# Each target exports its own include directories, so consumers only need
# target_link_libraries().

add_subdirectory(Meschach)      # meschach    - dense/sparse matrix library
add_subdirectory(gpc)           # gpc         - General Polygon Clipper
add_subdirectory(Eigen)         # eig3        - Kopp's 3x3 eigen-solvers (NOT the Eigen C++ library)
add_subdirectory(util)          # util        - small matrix helpers
add_subdirectory(triangulate)   # triangulate - Shewchuk's Triangle + wrappers
add_subdirectory(efsflib)       # EFSFLIB     - EFG shape functions
add_subdirectory(clipped_voronoi) # CVLIB     - clipped Voronoi diagrams
add_subdirectory(thread_pool)   # thread_pool - C-Thread-Pool
