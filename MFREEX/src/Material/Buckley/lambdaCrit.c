/*
 * =====================================================================================
 *
 *       Filename:  lambdaCrit.c
 *
 *    Description:  Calculates the critical network stretch in the conformational part
 *
 *        Version:  1.0
 *        Created:  08/05/18 11:12:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMTH 
 *   Organization:  QUEEN'S UNIVERSITY BELFAST 
 *
 * =====================================================================================
 */

#include "Material/Buckley/lambdaCrit.h"

double lambdaCrit(double critLambda_n, VEC * lambda, VEC * eigD, VEC * para){

	double critLambda = 1;
	VEC * strainRate = v_get(3);
	int  index = 0;
	double theta = 0;
	// find the strain rate

	for ( int i = 0 ; i < 3 ; i++){
		strainRate->ve[i] = eigD->ve[i]*lambda->ve[i] ;
	}

	double maxSr = v_max(strainRate,&index);


	if ( maxSr == 0){
		maxSr = 0.00000001; 
	}

	// find theta, the ratio of in plane natural strain rate
	if ( eigD->ve[1] == 0){
		theta = 0;
	}else{
		theta = eigD->ve[2]/eigD->ve[1];
	}

	// find xi 
	double xi = ( 2.00 * theta + 1.00)/(theta +2.00);
	if ( xi > 1){
		xi = 1;
	}else if( xi < -1 ) {
		xi = -1;
	}else{
		// do nothing
	}


	// update critical network stretch 
	double C1 = para->ve[0];
	double C2 = para->ve[1];
	double beta = para->ve[2];
	double k =  para->ve[3];
	double b = para->ve[4];
	double temperature = para->ve[5];


	double shift_factor = pow(10,(C1*(maxSr-1)/(C2 + maxSr -1))*pow(beta,2-2*xi));
	double shifted_temperature = temperature*shift_factor; 

	double critLambda_a = k * shifted_temperature + b;



	// test if necessary to overwrite crit lambda 

	if ( eigD->ve[1] < 0 ){
		critLambda = critLambda_n;
	}else{
		critLambda = critLambda_a; 
	}



	V_FREE(strainRate); 

	return critLambda;





}
