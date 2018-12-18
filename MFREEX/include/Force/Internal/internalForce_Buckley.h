#ifndef INTERNALFORCE_BUCKLEY_H_
#define INTERNALFORCE_BUCKLEY_H_

#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "Integration/SCNI/generate_scni.h"
#include "Material/material.h"
#include "Integration/defgrad.h"
#include "Material/Buckley/buckleyStress.h"


#include "Integration/gMat.h"
#define PI 3.14159265359



double internalForce_hyperelastic(VEC * Fint, SCNI_OBJ * scni_obj, VEC * disp, VEC * velocity,
	VEC * matParams, VEC * critLambdaParams, state_Buckley ** stateNew, state_Buckley ** stateOld,
	 int is_axi, int dim, double deltat);



#endif