#ifndef INTERNALFORCE_INELASTIC_H_
#define INTERNALFORCE_INELASTIC_H_

#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "Integration/SCNI/generate_scni.h"
#include "Material/material.h"
#include "Integration/defgrad.h"
#include "mat2csv.h"
#include "m_inverse_small.h"
#include "Integration/gMat.h"
#include "Deformation/velocity_grad.h"
#define PI 3.14159265359



double 
internalForce_Inelastic(VEC * Fint, SCNI_OBJ * scni_obj,
	VEC * disp, VEC * velocity,
	VEC * matParams, 
	state_variables ** stateNew, 
	state_variables ** stateOld,
	int is_axi, int dim, 
	double DT, double t_n_1, 
	char * Material);



#endif