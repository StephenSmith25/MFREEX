#ifndef MASS_MATRIX_H_
#define MASS_MATRIX_H_

#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include "Integration/material_point.h"

#define PI 3.14159265359


VEC *  mass_vector(MATERIAL_POINTS * MPS, meshfreeDomain * mfree);


#endif