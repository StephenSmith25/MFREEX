/*
 * =====================================================================================
 *
 *       Filename:  contraction.c
 *
 *    Description:  Performs contraction of two second order tensors
 *
 *        Version:  1.0
 *        Created:  04/05/18 11:29:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "contraction.h"


double contraction(MAT * A, MAT * B){

	double returnValue = 0;
	for ( int i = 0 ; i < A->m ; i++){
		for ( int j = 0 ; j < A->m; j++){
			returnValue += A->me[i][j]*B->me[i][j];
		}
	}
	return returnValue;
}
