/*
 * =====================================================================================
 *
 *       Filename:  SVK.h
 *
 *    Description:  St Venant Kirchhoff Material
 *
 *        Version:  1.0
 *        Created:  19/02/18 14:13:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Stephen Smith
 *   Organization:  Queen's University Belfast
 *
 * =====================================================================================
 */


#ifndef SVK_H_
#define SVK_H_
#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <stdio.h>
#include "trace.h"


int SVK(VEC * stressVoigt, MAT * defGrad, VEC * params); 






#endif 
