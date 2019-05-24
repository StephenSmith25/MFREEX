#ifndef STRESS_POINTS_H_
#define STRESS_POINTS_H_

#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"
#include "mls_shapefunction.h"
#include "generate_voronoi.h"
#include "Integration/SCNI/generate_scni.h"
#include "mat2csv.h"
#include "trigen.h"
#include "Integration/Triangle/triangle_quadrature.h"
#include "determinant.h"
#define PI 3.14159265358979323846




typedef struct STRESS_POINT
{
	VEC * phi;
	IVEC * neighbours;
	VEC * fInt;
	MAT * B;
	double volume;
	double r;
	double coords[3];

}STRESS_POINT;

STRESS_POINT ** create_stress_points(MAT * stress_point_coords,voronoi_diagram * voronoi, meshfreeDomain * mfree);





#endif