#ifndef INTEGRATION_STRUCTURE_H_
#define INTEGRATION_STRUCTURE_H_
#include <stdlib.h>
#include <stdio.h>


typedef struct {

	void * integration_points;
	int num_integration_points;

} integration_structure;



integration_structure * create_integration_points(int num_integration_points, void * integration_points);





#endif