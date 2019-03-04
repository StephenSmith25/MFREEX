#ifndef TRIANGLE_QUADRATURE_H_
#define TRIANGLE_QUADRATURE_H_
#include "matrix.h"
#include "matrix2.h"





typedef struct QUAD_TRIANGLE
{
	VEC * VOLUMES;
	MAT * QUAD_POINTS;

}QUAD_TRIANGLE;


QUAD_TRIANGLE * create_triangle_quadrature_points(double * points, int * triangles,
	int number_of_triangles, int * quad_orders );

#endif