/*
 * =====================================================================================
 *
 *       Filename:  buckleyConf.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/05/18 14:24:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#ifndef BUCKLEYCONF_H_
#define BUCKLEYCONF_H_
#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <stdio.h>
#include "Material/material.h"
#include "gamma.h"
#include <math.h>
#include "Material/Buckley/edwardsVilgis.h"

int buckleyConf(state_variables * stateNew, state_variables * stateOld, 
	VEC * para, double deltaT);








#endif 

