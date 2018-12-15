#ifndef NEW_PRESSURE_LOAD_H_
#define NEW_PRESSURE_LOAD_H_

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"



typedef struct pressure_boundary
{
	MAT * coords;
	IVEC * points;
	double * weights;
	double magnitude;
	shape_function_container * sf_traction;

} pressure_boundary;

pressure_boundary * new_pressure_load(IVEC * points, meshfreeDomain * mFree);



#endif