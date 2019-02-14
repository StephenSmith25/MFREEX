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
	double expFactor = para->ve[24];
	double gamma0 = 0;
	VEC * lambdaDot = state->lambdaDot;
	MAT * Dbar = state->Dbar;


	int dim = Dbar->m;
	double temperature = state->temperature;

	double gamma_n = state->gamma;

	// initialise gamma 
	double gamma_n_1 = 0; 

	// find theta 
	int index = 0;


	// //double maxSr = state->Vdot->me[2][2];

	// double V1,V2,V3 = 0;
	// V1 = state->Vdot->me[1][1];
	// V2 = state->Vdot->me[2][2];

	// double maxSr = max(V1,V2);



	double maxSr = v_max(state->lambdaDot,&index);
	if ( maxSr <= 0.01){
		maxSr = 0.01; 
	}

	MAT * d = state->dbar;
	double D1 = d->me[1][1];
	double D2 = d->me[2][2];


	//D1 = d->me[0][0];
	//D2 = d->me[1][1];



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
	double xi = ( 2 * theta + 1)/(theta +2);
	if ( xi > 1){
		xi = 1;
	}else if( xi < 0 ) {
		xi = 0.001;
	}else{
		// do nothing
	}



	double log2sr = log(maxSr)/log(2);
	if ( maxSr < 0.01){
		log2sr = 0.01;
	}


	double shiftTemperature = temperature * pow(10, ( ( shiftTSlope*log2sr ) / ( shiftTIntercept + log2sr) )* pow(expFactor,2-2*xi));
	//shiftTemperature = temperature;

	gamma0 = exp( Cs/(shiftTemperature - Tinf) - Cs/(starT - Tinf));
	gamma0 = gamma0*refGamma;


	if ( maxLambdaN >= critLambda ){
		gamma_n_1 = 1e30;
	}else{
		gamma_n_1 = gamma0 * ( critLambda - 1) / ( critLambda - maxLambdaN);
		//( 1.000 - (maxLambdaN/critLambda)) ; 
		//gamma_n_1 = gamma0/( 1 - maxLambdaN/critLambda);
	}


	
	state->gamma = gamma_n_1; 

	
	return gamma_n_1;
	


}

