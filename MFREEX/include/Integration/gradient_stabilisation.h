#ifndef GRADIENT_STABILISATION_H_
#define GRADIENT_STABILISATION_H_

#include <stdlib.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"
#include <math.h>
#include <omp.h>
#include "Integration/SCNI/generate_scni.h"



int gradient_stabilisation(SCNI_OBJ * scni_, shape_function_container * sf_point, meshfreeDomain * mfree  );


#endif
