#ifndef WEIGHT_FUNCTION_MATERIALPOINT_H_
#define WEIGHT_FUNCTION_MATERIALPOINT_H_

#include "matrix.h"
#include "matrix2.h"
#include <string.h>
#include "quartic_spline.h"
#include <math.h>
#include "cubic_spline.h"
#include "meshfree_structs.h"
#include "Integration/material_point.h"

int weight_function_materialpoint(VEC * weights, MATERIAL_POINT * MP, double  * xS, int dim);


#endif 