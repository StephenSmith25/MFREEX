#ifndef COMPRESSNEO_H_
#define COMPRESSNEO_H_




#include "Material/material.h"
#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <stdio.h>
#include "trace.h"
#include "determinant.h"
#include <math.h>
#include "m_inverse_small.h"

int compress_neo(VEC * stressVoigt, state_variables * states, VEC * params);




#endif
