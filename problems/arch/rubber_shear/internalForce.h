#ifndef INTERNALFORCE_H_
#define INTERNALFORCE_H_

#include "structure.h"
#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "smoothedDefGrad.h"
#define PI 3.14159265359


int internalForce(VEC * Fint, SCNI * scni,VEC * disp,VEC * matParams, int numnode );










#endif

