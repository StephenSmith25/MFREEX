#ifndef MATERIAL_POINT_H_
#define MATERIAL_POINT_H_

#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"
#include "mls_shapefunction.h"
#include "generate_voronoi.h"
#include "Integration/SCNI/generate_scni.h"
#include "mat2csv.h"
#include "trigen.h"
#include "Integration/Triangle/triangle_quadrature.h"


typedef struct MATERIAL_POINT
{


	// Geometric information
	double volume;
	double * coords;

	// shape functions
	VEC * phi;
	IVEC * neighbours;
	MAT * B;
	MAT * F_r;
	VEC * fInt;


	// density
	double rho;

	// material state 
	state_variables * state_n;
	state_variables * state_n_1;


}MATERIAL_POINT;

typedef struct MATERIAL_POINTS
{

	// Material points
	MATERIAL_POINT ** MP;
	
	// information about material points
	int num_material_points;
	int dim;
	int IS_AXI;

	// Raw coordinates of material points
	MAT * MP_COORDS;


	// Shape functions for each material points
	shape_function_container * sf_material_points;

}MATERIAL_POINTS;

typedef struct BOUNDARY_POINT
{




}BOUNDARY_POINT;



MATERIAL_POINTS * create_material_points(void * cells, int IS_AXI, int dim, 
	char * integration_type, char * material, double rho, meshfreeDomain * mfree);

int write_material_points(char * filename, MATERIAL_POINTS * MPS);





#endif