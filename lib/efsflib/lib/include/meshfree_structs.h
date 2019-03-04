#ifndef MESHFREE_STRUCTS_H_
#define MESHFREE_STRUCTS_H_

#include "matrix.h"
#include "matrix2.h"


enum SUPPORT_TYPE {RADIAL,RECTANGULAR,ELLIPTICAL};



typedef struct meshfreeDomain
{
	MAT * nodes;
	VEC * di;
	MAT * di_tensor;

	MAT ** MI;

	double dmax_radial ;
	double * dmax_tensor;
	double beta ;

	char * kernel_shape;
	char * basis_type;
	char * weight_function;


	enum SUPPORT_TYPE kernel_support;
	int num_nodes;
	int dim;
	int IS_AXI;
	int is_constant_support_size;
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