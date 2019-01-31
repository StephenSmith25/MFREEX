/*
 * =====================================================================================
 *
 *       Filename:  lambdaCrit.h
 *
 *    Description:  Calculate the critical network stretch 
 *
 *        Version:  1.0
 *        Created:  08/05/18 15:07:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef LAMBDACRIT_H_
#define LAMBDACRIT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"


double 
lambdaCrit(double critLambda_n, state_variables * state, VEC * para, double temperature, double dt);




#endif
