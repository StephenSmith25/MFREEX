#ifndef BUCKLEY_STRESS_H_
#define BUCKLEY_STRESS_H_

#include "Integration/SCNI/generate_scni.h"
#include "Deformation/velocity_grad.h"
#include "Material/material.h"
#include "Integration/defgrad.h"
#include "Deformation/poldec.h"
#include "Material/Buckley/buckleyBond.h"
#include "Material/Buckley/buckleyConf.h"
#include "Material/Buckley/lambdaCrit.h"
#include "matrix.h"
#include "matrix2.h"
#include <math.h>
#include "m_inverse_small.h"


int buckleyStress(state_Buckley * stateNew, state_Buckley * stateOld, 
	VEC * matParams, VEC * critLambdaParams, double dt, int i, int IS_AXI );



#endif