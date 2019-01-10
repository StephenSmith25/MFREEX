#ifndef ESSENTIAL_BOUNDARY_H_
#define ESSENTIAL_BOUNDARY_H_ 


#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"


typedef enum{
	PRESCRIBED_BOUNDARY,
	FIXED_BOUNDARY
}BOUNDARY_TYPE;

typedef struct ESSENTIAL_BOUNDARY 
{

	int * NB;
	int * NI;
	BOUNDARY_TYPE b;


}ESSENTIAL_BOUNDARY;




ESSENTIAL_BOUNDARY * new_essential_boundary(IVEC * boundary_nodes, BOUNDARY_TYPE b);




#endif