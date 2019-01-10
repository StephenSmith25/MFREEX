
#include "Boundary/Displacement/essential_boundary.h"

ESSENTIAL_BOUNDARY * new_essential_boundary(IVEC * boundary_nodes, BOUNDARY_TYPE b)
{

	ESSENTIAL_BOUNDARY *  eb = malloc(1*sizeof(ESSENTIAL_BOUNDARY));

	// Find NB and NI
	// NB 

	MAT * PHI_B = m_get(boundary_nodes->max_dim,boundary_nodes->max_dim);
	



	return NULL;


}
