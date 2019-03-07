#ifndef NEIGHBOURS_MATERIALPOINT_H_
#define NEIGHBOURS_MATERIALPOINT_H_

#include "matrix.h"
#include "matrix2.h"
#include <math.h>
#include "Integration/material_point.h"


IVEC * get_materialpoint_neighbours(IVEC * neighbours, MATERIAL_POINT * MP, MAT * nodes);



#endif