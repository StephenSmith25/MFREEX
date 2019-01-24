#ifndef INTERNALFORCE_HYPERELASTIC_S_H_
#define INTERNALFORCE_HYPERELASTIC_S_H_
#include "m_inverse_small.h"

#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "Integration/SCNI/generate_scni.h"
#include "Integration/SCNI/generate_mscni.h"

#include "Material/Hyperelastic/hyperelastic_materials.h"
#include "Integration/defgrad.h"

#define PI 3.14159265359


double internalForce_hyperelastic_S(VEC * Fint, MSCNI_OBJ * mscni_obj, VEC * disp, VEC * velocity,VEC * matParams, char * material, int is_axi, int dim);










#endif

