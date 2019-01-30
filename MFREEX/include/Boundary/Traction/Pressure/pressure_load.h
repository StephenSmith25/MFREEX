#ifndef PRESSURE_LOAD_H_
#define PRESSURE_LOAD_H_

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


typedef struct pressure_boundary
{
	MAT * coords;
	IVEC * points;
	double * segment_weights;
	MAT * segment_normals;
	double magnitude;
	shape_function_container * sf_traction;
	int is_axi;

} pressure_boundary;

pressure_boundary * new_pressure_boundary(IVEC * points, meshfreeDomain * mFree);
int assemble_pressure_load(VEC * Fext, double Pressure, pressure_boundary * pB);
int update_pressure_boundary(pressure_boundary *pB, MAT * coords);
int free_pressure_boundary(pressure_boundary * pB);



#endif