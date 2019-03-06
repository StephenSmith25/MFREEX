#ifndef MLS_SHAPEFUNCTION_H_
#define MLS_SHAPEFUNCTION_H_

#include "matrix.h"
#include "matrix2.h"
#include "meshfree_structs.h"
#include "neighbours.h"
#include "polynomial_basis.h"
#include "weight_function.h"
#include "v_outer_product.h"
#include <assert.h>
#include <omp.h>

shape_function_container * mls_shapefunction(MAT * compute_points, int compute, meshfreeDomain * mfree);
int free_shapefunction_container(shape_function_container * sf);

shape_function * MLS_SHAPEFUNCTION(shape_function * sf_point, double compute_point[3], int compute, meshfreeDomain * mfree);

int update_mls_shapefunction(shape_function * sf_point, double compute_point[3], int compute, meshfreeDomain * mfree);

#endif