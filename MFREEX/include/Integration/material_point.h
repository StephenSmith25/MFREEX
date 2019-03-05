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
#include "determinant.h"
#define PI 3.14159265358979323846


typedef struct MATERIAL_POINT_TRACTION
{



	double area;
	double * normal;

	
}MATERIAL_POINT_TRACTION;


typedef struct MATERIAL_POINTS_TRACTION
{
	// Material points
	MATERIAL_POINT_TRACTION ** MP;
	
	// information about material points
	int num_material_points;
	int dim;
	int IS_AXI;

	// Raw coordinates of material points
	MAT * MP_COORDS;


}MATERIAL_POINTS_TRACTION;


typedef struct MATERIAL_POINT
{

	// Geometric information
	double volume;
	double * coords;

	// shape functions
	VEC * phi;
	IVEC * neighbours;
	MAT * B;
	MAT * F_n;
	MAT * inc_F;
	VEC * fInt;
	VEC * stressVoigt;

	// density
	double rho;
	double temperature;
	double mass;

	// material state 
	state_variables * stateNew;
	state_variables * stateOld;


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

	// Material
	MATERIAL * material;




}MATERIAL_POINTS;

typedef struct BOUNDARY_POINT
{




}BOUNDARY_POINT;


MATERIAL_POINTS * create_material_points(void * cells, int IS_AXI, int dim, 
	char * integration_type, MATERIAL_TYPE material_type, double rho,  meshfreeDomain * mfree);

int write_material_points(char * filename, MATERIAL_POINTS * MPS);





#endif