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



double gammaV(state_Buckley * state, double maxLambdaN,double critLambda,
	VEC * para, double dt, int IS_AXI){


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

	double gamma = state->gamma;

	// initialise gamma 
	double gamma_n_1 = 0; 

	// find theta 
	int index = 0;
	double theta = 0;

	// find the strain rate
	double V1 = 0;
	double V2 = 0;
	if ( IS_AXI){
		V1 = state->Vdot->me[1][1]; 
		V2 = state->Vdot->me[2][2];
	}else{
		V1 = state->Vdot->me[0][0];
		V2 = state->Vdot->me[1][1];

	}
	double maxSr = max(V1,V2);

	if ( maxSr == 0){
		maxSr = 0.0000001; 
	}


	double D1 = state->eigValDBar->ve[2];
	double D2 = state->eigValDBar->ve[1];

	// find theta, the ratio of in plane natural strain rate
	if ( D1 == 0){
		theta = 0;
	}else{
		theta = D2/D1;
	}

	// find theta, the ratio of in plane natural strain rate
	if ( D1 == 0){
		theta = 0;
	}else{
		theta = D2/D1;
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
	gamma0 = exp( Cs/(shiftTemperature - Tinf) - Cs/(starT - Tinf));
	gamma0 = gamma0*refGamma;


	if ( maxLambdaN >= critLambda ){
		gamma_n_1 = 1e30;
	}else{
		gamma_n_1 = gamma0/ ( 1.000 - (maxLambdaN/critLambda)) ; 
	}


	
	state->gamma = gamma_n_1; 

	gamma_n_1 = gamma_n_1 ;
	
	return gamma_n_1;
	


}

