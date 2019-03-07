#ifndef GENERATE_BMAT_H_
#define GENERATE_BMAT_H_

#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"


MAT * generate_Bmat(MAT * phi_der, int dim, int is_axi, double r );


MAT * BMAT(MAT * Bmat, shape_function * basis_functions, int dim, int is_axi, double r );


#endif