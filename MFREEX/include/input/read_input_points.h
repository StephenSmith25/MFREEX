#ifndef READ_INPUT_POINTS_H_
#define READ_INPUT_POINTS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int read_input_points(char * fileName, double ** points, int ** boundary, 
	int ** boundary_markers, double ** attribute, int * num_nodes, int * num_boundary_nodes);



#endif