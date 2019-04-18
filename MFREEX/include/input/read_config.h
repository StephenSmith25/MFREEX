#ifndef READ_CONFIG_H_
#define READ_CONFIG_H_
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "BLOCKS.h"
#include "Boundary/Displacement/DisplacementBC.h"
#ifndef DIM
#define DIM 2
#endif


// read inputs into domain configuration
int read_config_file(DOMAIN * domain, char * filename);

// pressure conditions are read into sidesets
int read_pressure_loads();

// timestep parameters are read into the timestep structure
int read_timestep_parameters();

// material parameters are read into blockset 
int read_material_parameters();

// meshfree parameters are read into meshfree domain 
int read_meshfree_parameters();

// solution parameters kept in domain ( not sure if we'll use these in the end )
int read_solution_parameters();



#endif 