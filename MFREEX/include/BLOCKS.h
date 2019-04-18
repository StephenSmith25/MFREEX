#ifndef BLOCKS_H_
#define BLOCKS_H_



#include <stdlib.h>
#include "matrix.h"
#include <stdbool.h>
#include "Mesh/Elements.h"
#include "Boundary/Displacement/DisplacementBC.h"



typedef enum PHYSICAL_TYPE
{
	PRESSURE_TYPE=1,
	DISPLACEMENT_TYPE=2,
	PHYSICAL_MATERIAL_TYPE=3,

}PHYSICAL_TYPE;

typedef struct _PHYSICAL_GROUPS
{
	int ID;
	PHYSICAL_TYPE type;
	struct _PHYSICAL_GROUPS * next;
}PHYSICAL_NAME;


typedef struct _BLOCKSET 
{

	MAT * nodes;
	ELEMENT * elements;
	char * material;
	int ID; 

	struct _BLOCKSET * next;

}BLOCKSET;

typedef struct _NODESET
{
	ELEMENT * elements;
	DOF_CONSTRAINT * dof_constrained;
	BC_TYPE bc_type;

	// nodes which belong to nodeset
	IVEC * nodes;
	int ID; 

	struct _NODESET * next;

}NODESET;

typedef struct _SIDESET
{
	double pressure_magnitude;

	ELEMENT * elements;

	ELEMENT_TYPE element_type;

	double t_start;
	double t_end;
	bool RAMP; 
	int ID; 

	struct _SIDESET * next;

}SIDESET;



typedef struct DOMAIN
{


	MAT * NODES;

	// PHYSCIAL_GROUPS;
	PHYSICAL_NAME * physical_names;
	int number_of_physical_groups;
	// BLOCKSETS
	int NUM_BLOCK_SETS;
	BLOCKSET * blocksets;
	
	//NODESETS
	int NUM_NODE_SETS;
	NODESET * nodesets;

	//SIDESETS
	int NUM_SIDE_SETS;
	SIDESET * sidesets;




}DOMAIN;




// Domain constructors  
int AddBlockSetToDomain(DOMAIN * domain, int ID );
int AddSideSetToDomain(DOMAIN * domain, int ID );
int AddNodeSetToDomain(DOMAIN * domain, int ID );
int AddPhysicalNameToDomain(DOMAIN * domain, PHYSICAL_TYPE type, int ID);

PHYSICAL_NAME * FindPhysicalNameByID(DOMAIN * domain, int ID);

// Find sets
BLOCKSET * FindBlockSetByID(DOMAIN * domain, int ID);
SIDESET * FindSideSetByID(DOMAIN * domain, int ID);
NODESET * FindNodeSetByID(DOMAIN * domain, int ID);

// Add elements to blocksets
ELEMENT * AddElementToBlockSet(BLOCKSET * blockset, ELEMENT_TYPE etype, int * verticies);
ELEMENT * AddElementToSideSet(SIDESET * sideset, ELEMENT_TYPE etype, int * verticies);
ELEMENT * AddElementToNodeSet(NODESET * nodeset, ELEMENT_TYPE etype, int * verticies);

// Outputs
int WriteBlockSetElementsToFile(BLOCKSET * blockset, char * FILENAME);
int WriteSideSetElementsToFile(SIDESET * sideset, char * FILENAME);
int WriteNodeSetElementsToFile(NODESET * nodeset, char * FILENAME);

//DisplacmentBCS;
int AddDOFConstraintToNodeSet(NODESET * nodeset,DOF_CONSTRAINT * dof_constraint );
DOF_CONSTRAINT *  FindDOFConstraintInNodeSet(NODESET * nodeset, DIRECTION dir);

#endif 
