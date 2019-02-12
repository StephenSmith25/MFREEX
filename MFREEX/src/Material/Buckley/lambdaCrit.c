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

double lambdaCrit(double critLambda_n, state_variables * state, VEC * para, double temperature, double dt){
	

	VEC * lambdaDot = state->lambdaDot;

	double critLambda = 1;
	int  index = 0;
	int dim = lambdaDot->max_dim;



	// double V1,V2,V3 = 0;

	// V1 = state->Vdot->me[1][1];
	// V2 = state->Vdot->me[2][2];

	// double maxSr = max(V1,V2);
	double maxSr = v_max(state->lambdaDot,&index);
	if ( maxSr == 0){
		maxSr = 0.01; 
	}




	MAT * d = state->dbar;
	double D1 = d->me[1][1];
	double D2 = d->me[2][2];


	//  D1 = d->me[0][0];
	// D2 = d->me[1][1];



	double theta = 0;
	if ( D1 >= D2)
	{
		if ( D1 == 0 )
		{
			theta = 0;
		}else{
			theta = D2/D1;
		}
	}else {
		if ( D2 == 0)
		{
			theta = 0;
		}else{
			theta = D1/D2;
		}
	}


	// find xi 
	double xi = ( 2.00 * theta + 1.00)/(theta +2.00);
	if ( xi > 1){
		xi = 1;
	}else if( xi < -1) {
		xi = -1;
	}else{
		// do nothing
	}


	// update critical network stretch 
	double C1 = para->ve[26];
	double C2 = para->ve[27];
	double beta = para->ve[28];
	double k =  para->ve[29];
	double b = para->ve[30];


	double shift_factor = pow(10,(C1*(maxSr-1)/(C2 + maxSr -1))*pow(beta,2-2*xi));
	double shifted_temperature = temperature*shift_factor;
	shifted_temperature = temperature;
	double critLambda_a = k * shifted_temperature + b;


	if ( (D2 < 0) && ( D1 < 0) ){
		// don't update the value
		critLambda = critLambda_n;
	}else{
		// update the value
		critLambda = critLambda_a; 
	}


	return critLambda;





}
// double D1 = 0;
// 	double D2 = 0;
// 	if ( IS_AXI == 1)
// 	{
// 		D1 = Dbar->me[1][1];
// 		D2 = Dbar->me[2][2];

// 		if ( D1 >= D2)
// 		{
// 			// correct order
// 		}else{
// 			double temp = D1;
// 			D1 = D2;
// 			D2 = temp;
// 		}


// 	}else{
// 		D1 = Dbar->me[0][0];
// 		D2 = Dbar->me[1][1];
// 		if ( D1 >= D2)
// 		{
// 			// correct order
// 		}else{
// 			double temp = D1;
// 			D1 = D2;
// 			D2 = temp;
// 		}
// 	}