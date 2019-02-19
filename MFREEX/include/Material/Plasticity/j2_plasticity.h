#ifndef J2_PLASTICITY_H_
#define J2_PLASTICITY_H_


#include "Material/material.h"
#include "matrix.h"
#include "matrix2.h"
#include "contraction.h"
#include <math.h>

int j2_plasticity(state_variables * stateNew, state_variables * stateOld, VEC * params, double dt);


#endif
