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
	MAT * Dbar = state->dbar;


	// Find maximum strain rate 

	double V1 = state->Vdot->me[0][0];
	double V2 = state->Vdot->me[1][1];
	double V3 = state->Vdot->me[2][2];

	double maxSr = max(V1,V2);
	maxSr = max(maxSr,V3);

	//double maxSr = v_max(state->lambdaDot, &index);

	if ( maxSr == 0){
	 	maxSr = 0.01; 
	}



	double D1 = 0;
	double D2 = 0;
	
	D1 = Dbar->me[1][1];
	D2 = Dbar->me[2][2];

	if ( D1 >= D2)
	{
		// correct order
	}else{
		
		double temp = D1;
		D1 = D2;
		D2 = temp;
	}


	// Find deformation mode indicator 
	double xi = 0;
	double theta = 0;
	if (D1 > 0)
	{
		theta = D2/D1;
	}

	xi = (2*theta + 1.00)/(theta+2.00);

	if ( xi > 0)
	{
		xi = 1;
	}else if ( xi < -1)
	{
		xi = -1;
	}


	// // update critical network stretch 
	double C1 = para->ve[26];
	double C2 = para->ve[27];
	double beta = para->ve[28];
	double k =  para->ve[29];
	double b = para->ve[30];

	// Shifted temperature 
	double shift_factor = pow(10,(C1*(maxSr-1.00)/(C2 + maxSr -1.00))*pow(beta,2-2*xi)   );
	double shifted_temperature = temperature*shift_factor;
	critLambda = k * shifted_temperature + b;

	// return
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