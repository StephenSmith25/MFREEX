
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  yeoh
 *  Description:  Finds the stress using the hyperelastic Yeoh material model 
 * =====================================================================================
 */
// 22/11/17 09:52:32
//
//
#ifndef YEOH_H_
#define YEOH_H_

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include <math.h>
#include "determinant.h"
#include "trace.h"

int yeoh(VEC * stressVoigt, MAT * defGrad, VEC * params); 




#endif




