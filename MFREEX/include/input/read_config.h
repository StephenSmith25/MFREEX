#ifndef READ_CONFIG_H_
#define READ_CONFIG_H_
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "BLOCKS.h"
#include "Boundary/Displacement/DisplacementBC.h"
#include "Boundary/Traction/Pressure/pressure_load.h"
#ifndef DIM
#define DIM 2
#endif


// read inputs into domain configuration
int read_config_file(DOMAIN * domain, char * filename);


// meshfree parameters are read into meshfree domain 
int read_meshfree_parameters();

// solution parameters kept in domain ( not sure if we'll use these in the end )
int read_solution_parameters();



#endif 