#include "Integration/integration_structure.h"


integration_structure * create_integration_points(int num_integration_points, void * integration_points)
{

	integration_structure * int_struct = malloc(1*sizeof(integration_structure));

	int_struct->num_integration_points = num_integration_points;
	int_struct->integration_points = integration_points;


	return int_struct;

}





