#ifndef GRADIENT_STABILISATION_H_
#define GRADIENT_STABILISATION_H_

#include <stdlib.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include <math.h>
#include <omp.h>

typedef struct stabalised_gradient
{
	MAT * g ;

}stabalised_gradient;

stabalised_gradient **  gradient_stabilisation( shape_function_container * sf_point, meshfreeDomain * mfree  );


#endif
