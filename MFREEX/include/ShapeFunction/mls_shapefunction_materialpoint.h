#ifndef MLS_SHAPEFUNCTION_MATERIALPOINT_H_
#define MLS_SHAPEFUNCTION_MATERIALPOINT_H_

#include "matrix.h"
#include "matrix2.h"
#include "meshfree_structs.h"
#include "polynomial_basis.h"
#include "v_outer_product.h"
#include <assert.h>
#include <omp.h>
#include "Integration/material_point.h"
#include "ShapeFunction/neighbours_materialpoint.h"
#include "ShapeFunction/weight_function_materialpoint.h"


shape_function * mls_shapefunction_materialpoint(MATERIAL_POINT * MP, int compute, MAT * nodes);


#endif