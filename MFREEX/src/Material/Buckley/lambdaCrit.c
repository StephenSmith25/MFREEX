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

double lambdaCrit(double critLambda_n, state_Buckley * state, VEC * para, double temperature, double dt, int IS_AXI){
	

	VEC * lambdaDot = state->lambdaDot;

	double critLambda = 1;
	int  index = 0;
	double theta = 0;
	int dim = lambdaDot->max_dim;
	// find the strain rate
	// double V1 = 0;
	// double V2 = 0;
	// if ( IS_AXI){
	// 	V1 = state->Vdot->me[1][1]; 
	// 	V2 = state->Vdot->me[2][2];
	// }else{
	// 	V1 = state->Vdot->me[0][0];
	// 	V2 = state->Vdot->me[1][1];

	// }
	// double maxSr = max(V1,V2);

	//find the strain rate
	double V1 = 0;
	double V2 = 0;
	if ( IS_AXI == 1){
		V1 = state->Vdot->me[1][1]; 
		V2 = state->Vdot->me[2][2];
	}else{
		V1 = state->Vdot->me[0][0];
		V2 = state->Vdot->me[1][1];

	}
	double maxSr = max(V1,V2);

	//double maxSr = v_max(state->lambdaDot,&index);
	if ( maxSr == 0){
		maxSr = 0.01; 
	}

	// double D1 = state->eigValDBar->ve[2];
	// double D2 = state->eigValDBar->ve[1];
	//find the strain rate
	double D1 = 0;
	double D2 = 0;
	if ( IS_AXI == 1){
		D1 = state->d->me[1][1]; 
		D2 = state->d->me[2][2];
	}else{
		D1 = state->d->me[0][0];
		D2 = state->d->me[1][1];

	}


	// // find xi 
	// double xi = ( 2.00 * theta + 1.00)/(theta +2.00);
	// if ( xi > 1){
	// 	xi = 1;
	// }else if( xi < -1) {
	// 	xi = -1;
	// }else{
	// 	// do nothing
	// }
	double xi = 1.00;

	// update critical network stretch 
	double C1 = para->ve[0];
	double C2 = para->ve[1];
	double beta = para->ve[2];
	double k =  para->ve[3];
	double b = para->ve[4];


	double shift_factor = pow(10,(C1*(maxSr-1)/(C2 + maxSr -1))*pow(beta,2-2*xi));
	double shifted_temperature = temperature*shift_factor;

	double critLambda_a = k * shifted_temperature + b;


	if ( (D2 < 0) && ( D1 < 0) ){
		// don't update the value
		critLambda = critLambda_n;
	}else{
		// update the value
		critLambda = critLambda_a; 
	}

	critLambda = critLambda_a; 


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