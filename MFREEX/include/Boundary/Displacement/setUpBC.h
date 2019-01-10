/*
 * =====================================================================================
 *
 *       Filename:  setUpBC.h
 *
 *    Description:  Enforces essential boundary conditions using Joldes (2017)
 *
 *        Version:  1.0
 *        Created:  14/02/18 09:54:29
 *       Revision:  1
 *       Compiler:  gcc
 *
 *         Author:  Stephen Smith 
 *   Organization:  Queen's University Belfast
 *
 * =====================================================================================
 */


#ifndef SETUPBC_H_
#define SETUPBC_H_

#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix2.h"
#include "mls_shapefunction.h"


typedef struct EBC {
	MAT * V;
	MAT * phi;
	IVEC * nodes;
	MAT * coords; 
	MAT * P; 
	int dofFixed;
	VEC * uBar1;
	VEC * uBar2;
	VEC * uCorrect1;
	VEC * uCorrect2;
}EBC;

int setUpBC(EBC * ebc_, VEC * invMass, meshfreeDomain * mfree );








#endif
