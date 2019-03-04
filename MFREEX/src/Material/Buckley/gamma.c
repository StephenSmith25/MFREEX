/*
 * =====================================================================================
 *
 *       Filename:  gamma.c
 *
 *    Description:  Finds the viscosity of the dashpot used in the conformational part of
 *    		    the Buckley model. 
 *
 *        Version:  1.0
 *        Created:  10/05/18 17:56:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  STEPHEN SMITH 	 
 *   Organization:  QUEEN'S UNIVERSITY BELFAST 
 *
 * =====================================================================================
 */

#include "Material/Buckley/gamma.h"



double gammaV(state_variables * state, double maxLambdaN,double critLambda,
	VEC * para, double dt){




	double starT = para->ve[18];
	double refGamma = para->ve[19];
	double Cs = para->ve[20];
	double Tinf = para->ve[21];
	double shiftTSlope = para->ve[22];
	double shiftTIntercept = para->ve[23];
	double beta = para->ve[24];
	double gamma0 = 0;
	int indx = 0;
	double temperature = state->temperature;
	MAT * Dbar = state->dbar;


	double gamma_n = state->gamma;
	// initialise gamma 
	double gamma_n_1 = 0; 

	// Find maximum strain rate 

	double V1 = state->Vdot->me[0][0];
	double V2 = state->Vdot->me[1][1];
	double V3 = state->Vdot->me[2][2];

	double maxSr = max(V1,V2);
	maxSr = max(maxSr,V3);

	//double maxSr = v_max(state->lambdaDot, &indx);

	double log2sr = log2(maxSr);
	if ( maxSr < 0.01){
		log2sr = 0;
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

	}else if ( xi < 0)
	{
		xi = 0.001;
	}

	xi = 0;



	double shiftTemperature = temperature * 
	pow(10, ( ( shiftTSlope*log2sr ) / ( shiftTIntercept + log2sr))*pow(beta,2-2*xi)   );
	gamma0 = exp( Cs/(shiftTemperature - Tinf) - Cs/(starT - Tinf));
	gamma0 = gamma0*refGamma;


	if ( maxLambdaN >= critLambda ){
		gamma_n_1 = 1e50;
	}else{
		gamma_n_1 = gamma0/( 1 - maxLambdaN/critLambda);
	}


	
	state->gamma = gamma_n_1; 

	
	return gamma_n_1;
	


}

