project(Eigen)

SET(EIGEN_INTERNAL_SRC
	slvsec3.c
	dsyev2.c 
	dsytrd3.c
	dsyevc3.c
	dsyevj3.c
	dsyevq3.c
	dsyevv3.c 
	dsyevh3.c 
	dsyevd3.c 
	zhetrd3.c 
	zheevc3.c 
	zheevj3.c 
	zheevq3.c 
	zheevv3.c 
	zheevh3.c 
	zheevd3.c
)


FILE(GLOB EIGEN_HEADERS *.h)

SET(EIGEN_LIBRARY_SOURCE
	${EIGEN_INTERNAL_SRC}
	${EIGEN_HEADERS})

add_definitions("-O3 -fPIC -Wall -std=c99 -c ")



add_library(Eigen STATIC ${EIGEN_LIBRARY_SOURCE})