#ifndef INTERNALFORCE_HYPERELASTIC_H_
#define INTERNALFORCE_HYPERELASTIC_H_

#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "smoothedDefGrad.h"
#define PI 3.14159265359


int internalForce_hyperlastic(VEC * Fint, SCNI * scni,VEC * disp,VEC * matParams, int numnode );










#endif

