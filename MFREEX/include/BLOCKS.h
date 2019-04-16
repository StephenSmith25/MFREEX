#ifndef BLOCKS_H_
#define BLOCKS_H_



#include <stdlib.h>
#include "matrix.h"
#include <stdbool.h>
#include "Mesh/Elements.h"

typedef enum PHYSICAL_GROUP
{
	PRESSURE_TYPE=1,
	DISPLACEMENT_TYPE=2,
	PHYSICAL_MATERIAL_TYPE=3,

}PHYSICAL_GROUP;


typedef enum BC_TYPE
{
	DOF_FIXED = 1,
	DOF_PRESCRIBED=2,
}BC_TYPE;

typedef enum DOF_CONSTRAINT
{
	X = 1,
	Y = 2,
	Z = 3,
	ALL=4
}DOF_CONSTRAINT;

typedef struct BLOCKSET 
{

	MAT * nodes;
	MAT * ELEMENTS;
	ELEMENT_TYPE element_type;
	char * material;
	int ID;

}BLOCKSET;

typedef struct NODESET
{

	DOF_CONSTRAINT dof_constrained;
	BC_TYPE bc_type;

	int ID;
}NODESET;

typedef struct SIDESET
{
	double pressure_magnitude;
	MAT * ELEMENTS;
	ELEMENT_TYPE element_type;

	double t_start;
	double t_end;
	bool RAMP; 

}SIDESET;


#endif

typedef struct DOMAIN
{


	MAT * NODES;

	// PHYSCIAL_GROUPS;
	PHYSICAL_GROUP * physical_groups;
	int number_of_physical_groups;
	// BLOCKSETS
	int NUM_BLOCK_SETS;
	BLOCKSET ** blocksets;
	
	//NODESETS
	int NUM_NODE_SETS;
	NODESET ** nodesets;

	//SIDESETS
	int NUM_SIDE_SETS;
	SIDESET ** sidesets;




}DOMAIN;