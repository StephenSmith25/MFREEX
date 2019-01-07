#ifndef MESHFREE_STRUCTS_H_
#define MESHFREE_STRUCTS_H_

#include "matrix.h"
#include "matrix2.h"

typedef struct meshfreeDomain
{
	MAT * nodes;
	VEC * di;
	int num_nodes;
	int dim;
	int IS_AXI;
} meshfreeDomain;

typedef struct shape_function{


	// construction point
	double * x;
	IVEC * neighbours;
	VEC * phi;


	MAT * dphi;
	MAT * d2phi;


} shape_function;

typedef struct shape_function_container{

	int number_of_points;
	int compute;
	shape_function ** sf_list;

} shape_function_container;
#endif