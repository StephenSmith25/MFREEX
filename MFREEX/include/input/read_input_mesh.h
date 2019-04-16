#ifndef READ_INPUT_MESH_H
#define READ_INPUT_MESH_H


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "BLOCKS.h"

#ifndef DIM
#define DIM 2
#endif

DOMAIN * read_intput_mesh(char * input_file_name);

#endif 