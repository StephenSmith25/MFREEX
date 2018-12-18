/*
 * =====================================================================================
 *
 *       Filename:  gMat.c
 *
 *    Description:  Computs the G matrix such that {P} = J[G]{sigma} where sigma is the 
 *    		    Cauchy stress, {.} denotes the Voigt form of a second order tensor, 
 *    		    [.] a matrix and J the Jacobian, defined as the determiant of the 
 *    		    deformation gradient det[F] 
 *
 *        Version:  1.0
 *        Created:  04/05/18 16:30:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef GMAT_H_
#define GMAT_H_
#include <stdlib.h>

#include "matrix.h"
#include "matrix2.h"
#include <math.h>


int gMat(MAT * G, MAT * invDefGrad, int is_AXI);

#endif 
