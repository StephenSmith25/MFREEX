/*
 * =====================================================================================
 *
 *       Filename:  gMat.c
 *
 *    Description:  Computes the G matrix such that {P} = J[G]{sigma}
 *    		    Where sigma is the Cauchy stress, and {.} denotes Voigt form of a second
 *    		    order tensor, J the Jacobian of the deformation gradient and [.] a matrix.
 *
 *        Version:  1.0
 *        Created:  04/05/18 15:58:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Stephen Smith 
 *   Organization:  Queen's University Belfast 
 *
 * =====================================================================================
 */

#include "Integration/gMat.h"
int gMat(MAT * G, MAT * invDefGrad, int is_AXI){

	m_zero(G);


	// row 1
	G->me[0][0] = invDefGrad->me[0][0];
	G->me[0][2] = invDefGrad->me[0][1];
	// row 2
	
	G->me[1][1] = invDefGrad->me[1][1];
	G->me[1][2] = invDefGrad->me[1][0];
	// row 3
	G->me[2][0] = invDefGrad->me[1][0];
	G->me[2][2] = invDefGrad->me[1][1];
	// row 4
	G->me[3][1] = invDefGrad->me[0][1];
	G->me[3][2] = invDefGrad->me[0][0];
	// row 5

	if ( is_AXI == 1){

	
	G->me[4][3] = invDefGrad->me[2][2];
	}






	return 0;
}
