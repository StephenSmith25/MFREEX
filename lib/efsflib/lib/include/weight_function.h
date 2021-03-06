#ifndef WEIGHT_FUNCTION_H_
#define WEIGHT_FUNCTION_H_

#include "matrix.h"
#include "matrix2.h"
#include <string.h>
#include "quartic_spline.h"
#include <math.h>
#include "cubic_spline.h"
#include "meshfree_structs.h"

int weight_function (VEC * weights, double  * xS, int I, meshfreeDomain * mfree,  int compute);


#endif 