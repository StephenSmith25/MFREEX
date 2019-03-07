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

	double * coords;
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
	double INTEGRATION_FACTOR;
	double * coords;

	// shape functions
	MAT * B;
	MAT * F_n;
	MAT * inc_F;
	VEC * fInt;
	VEC * stressVoigt;

	// density
	double rho;
	double temperature;
	double mass;
	double Jn ;


	shape_function * shape_function;

	// material state 
	state_variables * stateNew;
	state_variables * stateOld;

	// Support
	double beta;
	MAT * MI;
	MAT * invMI;
	enum SUPPORT_TYPE kernel_support;

}MATERIAL_POINT;

typedef struct MATERIAL_POINTS
{

	// Material points
	MATERIAL_POINT ** MP;
	
	// information about material points
	int num_material_points;
	int dim;
	int IS_AXI;


	// Material
	MATERIAL * material;




}MATERIAL_POINTS;

typedef struct BOUNDARY_POINT
{




}BOUNDARY_POINT;


MATERIAL_POINTS * create_material_points(void * cells, int IS_AXI, int dim, 
	char * integration_type, MATERIAL_TYPE material_type, double rho, double beta, meshfreeDomain * mfree);

int write_material_points(char * filename, MATERIAL_POINTS * MPS);

MATERIAL_POINT * update_material_point(MATERIAL_POINT * MP, MAT * NODES, VEC * nodal_mass);

void set_material_point_temperature(double temperature , MATERIAL_POINT * MP);




#endif