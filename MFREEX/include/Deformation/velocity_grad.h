#ifndef VELOCITY_GRAD_h_
#define VELOCITY_GRAD_h_


#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"
#include "Deformation/incremental_deformation_gradient.h"



int velocity_grad(state_variables * stateNew, state_variables * stateOld, 
double delta_t,double alpha);



#endif