/*
 * =====================================================================================
 *
 *       Filename:  buckleyBond.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/05/18 12:14:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef BUCKLEYBOND_H_
#define BUCKLEYBOND_H_

#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "matrix2.h"
#include "Material/material.h"
#include "contraction.h"

int buckleyBond(state_variables * stateNew, state_variables * stateOld , VEC * para, double dt);


#endif
