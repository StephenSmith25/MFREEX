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
#include "cellSearch.h"

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
	double * coords_n_1;
	double * coords_n;

	// shape functions
	MAT * B;
	MAT * F_n;
	MAT * inc_F;
	VEC * fInt;
	VEC * stressVoigt;
	MAT * temp;
	MAT * temp_1;

	// density
	double rho;
	double temperature;
	double mass;
	double Jn ;


	shape_function * shape_function;

	int num_neighbours;
	double r_cutoff ;
	IVEC * neighbours;

	// material state 
	state_variables * stateNew;
	state_variables * stateOld;

	// Support
	double beta;
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


typedef struct UPDATE_STRUCT
{
	MAT * nodes;
	MATERIAL_POINT * MP;
	CELLS * grid;
}UPDATE_STRUCT;

MATERIAL_POINTS * create_material_points(void * cealls, int IS_AXI, int dim, 
	char * integration_type, MATERIAL_TYPE material_type, CELLS * cells, double rho, double beta, meshfreeDomain * mfree);

int write_material_points(char * filename, MATERIAL_POINTS * MPS);

MATERIAL_POINT * update_material_point(MATERIAL_POINT * MP, CELLS * grid, MAT * NODES, VEC * nodal_mass);
void * update_material_point_a(void * threadarg);

void set_material_point_temperature(double temperature , MATERIAL_POINT * MP);
void update_material_point_coords(MATERIAL_POINT * MP, MAT * NODES);




#endif