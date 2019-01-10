/*
 * =====================================================================================
 *
 *       Filename:  enforceBC.h
 *
 *    Description:  Enforces the boundary conditions using Joldes (2017 )
 *
 *        Version:  1.0
 *        Created:  14/02/18 16:09:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Stephen Smith	 
 *   Organization:  Queen's University Belfast
 *
 * =====================================================================================
 */



#include "matrix.h"
#include "matrix2.h"
#include <stdlib.h>
#include <stdio.h>
#include "Boundary/Displacement/setUpBC.h"


int enforceBC(EBC * ebc_, VEC * disp);
