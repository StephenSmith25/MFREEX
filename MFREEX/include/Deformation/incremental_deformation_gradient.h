#ifndef INCREMENTAL_DEFORMATION_GRADIENT_H_
#define INCREMENTAL_DEFORMATION_GRADIENT_H_


#include "Material/material.h"
#include "matrix.h"
#include "m_inverse_small.h"
#include "matrix2.h"

MAT *  
incremental_deformation_gradient(state_variables * stateNew, state_variables * stateOld,
double alpha);


#endif