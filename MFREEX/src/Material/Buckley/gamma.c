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

double gammaV(state_Buckley * state, double maxLambdaN,double critLambda,VEC * para, double dt){


	double starT = para->ve[18];
	double refGamma = para->ve[19];
	double Cs = para->ve[20];
	double Tinf = para->ve[21];
	double shiftTSlope = para->ve[22];
	double shiftTIntercept = para->ve[23];
	double expFactor = para->ve[24];
	double gamma0 = 0;
	VEC * lambdaDot = state->lambdaDot;
	VEC * eigD = state->eigValDBar;


	int dim = eigD->max_dim;
	double temperature = state->temperature;


	// initialise gamma 
	double gamma_n_1 = 0; 

	// find theta 
	int index = 0;
	double theta = 0;

	// find the strain rate
	double maxSr = v_max(lambdaDot,&index);

	if ( maxSr == 0){
		maxSr = 0.0000001; 
	}



	if ( eigD->ve[1] == 0){
		theta = 0;
	}else{
		theta = eigD->ve[2]/eigD->ve[1];
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
	if ( maxSr < 0.001){
		log2sr = 0;
	}


	double shiftTemperature = temperature * pow(10, ( ( shiftTSlope*log2sr ) / ( shiftTIntercept + log2sr) )* pow(expFactor,2-2*xi));
	gamma0 = exp( Cs/(shiftTemperature - Tinf) - Cs/(starT - Tinf));
	gamma0 = gamma0*refGamma;


	if ( maxLambdaN >= critLambda ){
		gamma_n_1 = 1e50;

	}else{
		gamma_n_1 = gamma0/ ( 1.000 - (maxLambdaN/critLambda)) ; 
	}

	
	
	gamma_n_1 = gamma_n_1 ;
	
	return gamma_n_1;
	


}

