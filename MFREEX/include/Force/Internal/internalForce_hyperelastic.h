#ifndef INTERNALFORCE_HYPERELASTIC_H_
#define INTERNALFORCE_HYPERELASTIC_H_

#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "determinant.h"
#include "Integration/SCNI/generate_scni.h"
#include "Material/Hyperelastic/hyperelastic_materials.h"
#include "Integration/defgrad.h"
#include "m_inverse_small.h"
#include "mat2csv.h"
#include "Integration/material_point.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

double internalForce_hyperelastic(VEC * Fint, MATERIAL_POINT * MP, VEC * disp, VEC * velocity,
 	VEC * matParams, int (*mat_func_ptr)(VEC *, MAT*,VEC*), double t_n_1);










#endif

