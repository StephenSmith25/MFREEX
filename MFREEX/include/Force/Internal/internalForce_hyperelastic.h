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

#define PI 3.14159265359


double internalForce_hyperelastic(VEC * Fint, MATERIAL_POINTS * material_points, VEC * disp, VEC * velocity, VEC * matParams, 
	char * Material, int is_axi, int dim, double t_n_1);









#endif

